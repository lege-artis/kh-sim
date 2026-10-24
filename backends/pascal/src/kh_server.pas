program kh_server;
{ KH-SIM Pascal backend - KH-007
  Framework : fphttpserver (Free Pascal fcl-web low-level)
  Port      : 8005
  Endpoints : POST /simulate  GET /health  GET /info
  Note      : Uses TFPHttpServer.OnRequest directly - no module registration.
              Server.Active := True blocks in non-threaded mode (accept loop).

  Error contract (shared/api/openapi.yaml, responses BadRequest/InternalError):
    400 {"error", "field"}    the request cannot be computed as given
    500 {"error", "backend"}  the computation itself failed
  Until 2026-09-26 every exception became a 400 carrying an FPU trap text, and
  a non-power-of-two grid wrote outside the FFT work array, corrupting the heap
  and leaving the process alive but no longer answering any request, /health
  included. See development record, finding 2026-09-16, B-3.

  Contract KH-VAL-2 (development record, modelling rules v0.1, ruling 2026-10-07)
  adds to KH-VAL-1, same rules as backends/cpp/src/kh_http_contract.hpp:
    K3 outside the model: domain_lx, domain_ly, reynolds_number < 0; dt < 0;
       steps < 0; linear_drag < 0; perturbation_mode < 0; band forcing with
       forcing_kmin > forcing_kmax                              -> 400 field
    K4 not representable (analytic field, amplitude <> 0): perturbation_mode
       0, or beyond the 2/3-rule band (3*mode > 2*(grid_nx div 2)) -> 400 field
    K4b (ruling 2026-10-09) a non-zero perturbation_amplitude sent with
       initial_omega (it would be silently ignored)           -> 400 field
    K4 forcing (amplitude <> 0, 2026-10-10): mode forcing_kx/ky beyond the
       2/3-rule band, or a band with 3*floor(kmax) > 2*(min(nx, ny) div 2) -> 400 field
    K3 a band with no integer mode kmin <= |k| <= kmax        -> 400 forcing_kmax
    K5 "warnings" in the response (KhWarnings).
    K7 (2026-10-24) a key the contract does not define, checked first; with
       several, the first in byte order; foreign names get a hint -> 400 field
  t_start (2026-10-10): the time of the initial field (forcing phase, t_final).
}
{$mode objfpc}{$H+}

uses
  SysUtils, Classes, Math,
  fphttpserver, httpdefs,
  fpjson, jsonparser,
  kh_physics;

const
  BACKEND_ID  = 'kh-pascal';             { the backend field of a 500 body, as kh-cpp, kh-rust ... }
  { One name for this server in /health, /info and every /simulate result, as the
    other four ("cpp-httplib", "rust-axum", ...). Until 2026-10-11 /health said
    "pascal" and /info and /simulate said "kh-pascal". }
  SERVICE_ID  = 'pascal-fphttpserver';
  SERVER_PORT = 8005;
  { Largest grid accepted, in cells. At 2048 x 2048 each complex work array is
    64 MB and the peak working set is roughly 1 GB. Every grid used by the test
    suites (up to 256 x 256) is far inside it.                                }
  MAX_CELLS   = 4194304;
  { Per-axis ceiling, so the cell product cannot overflow before it is checked. }
  MAX_AXIS    = 65536;

{ -- JSON helpers ----------------------------------------------------------- }
{ 2026-10-03 (finding v0.5, P-5): a key present with the wrong JSON type is a
  400 naming that key, as in the other four backends (contract KH-VAL-1). Until
  then fpjson converted silently ("32" became 32) or raised a conversion error
  reported against "body". A JSON null means the key is absent.               }
function JPresent(jo: TJSONObject; const key: string): TJSONData;
begin
  Result := jo.Find(key);
  if Assigned(Result) and (Result.JSONType = jtNull) then Result := nil;
end;

function JGetD(jo: TJSONObject; const key: string; def: Double): Double;
var jd: TJSONData;
begin
  jd := JPresent(jo, key);
  if not Assigned(jd) then Exit(def);
  if jd.JSONType <> jtNumber then
    raise EKhInvalid.CreateField(key, key + ' must be a number');
  Result := jd.AsFloat;
end;

function JGetI(jo: TJSONObject; const key: string; def: Integer): Integer;
var jd: TJSONData; f: Double;
begin
  jd := JPresent(jo, key);
  if not Assigned(jd) then Exit(def);
  if jd.JSONType <> jtNumber then
    raise EKhInvalid.CreateField(key, key + ' must be an integer');
  f := jd.AsFloat;
  if (f <> Int(f)) or (Abs(f) > 2147483647.0) then
    raise EKhInvalid.CreateField(key, key + ' must be an integer in the 32-bit range');
  Result := Trunc(f);
end;

function JGetS(jo: TJSONObject; const key: string; const def: string): string;
var jd: TJSONData;
begin
  jd := JPresent(jo, key);
  if not Assigned(jd) then Exit(def);
  if jd.JSONType <> jtString then
    raise EKhInvalid.CreateField(key, key + ' must be a string');
  Result := jd.AsString;
end;

procedure JAddArr(jo: TJSONObject; const key: string; const data: TRealArr);
var a: TJSONArray; k: Integer;
begin
  a := TJSONArray.Create;
  for k := 0 to Length(data)-1 do a.Add(data[k]);
  jo.Add(key, a);
end;

{ Error bodies are built with fpjson so that quotes, backslashes and control
  characters in exception messages are escaped correctly.                    }
function ErrorJSON(const Msg, KeyName, KeyValue: string): string;
var jo: TJSONObject;
begin
  jo := TJSONObject.Create;
  try
    jo.Add('error', Msg);
    jo.Add(KeyName, KeyValue);
    Result := jo.AsJSON;
  finally
    jo.Free;
  end;
end;

{ -- Request parsing -------------------------------------------------------- }
{ KH-VAL-2, K7 (step 2, 2026-10-24): a key the contract does not define is refused
  before any value is read; with several, the first in byte order is named. }
procedure CheckUnknownKeys(jo: TJSONObject);
const
  KNOWN: array[0..20] of string = (
    'grid_nx', 'grid_ny', 'domain_lx', 'domain_ly', 'dt', 'steps', 'reynolds_number',
    'velocity_shear', 'perturbation_amplitude', 'perturbation_mode', 'initial_omega',
    'init_omega', 'forcing_type', 'forcing_amplitude', 'forcing_kx', 'forcing_ky',
    'forcing_kmin', 'forcing_kmax', 'forcing_omega0', 'linear_drag', 't_start');
  FOREIGN_FROM: array[0..11] of string = (
    'nsteps', 'n_steps', 'num_steps', 'nx', 'ny', 'lx', 'ly', 're', 'Re', 'reynolds', 'amplitude', 'mode');
  FOREIGN_TO: array[0..11] of string = (
    'steps', 'steps', 'steps', 'grid_nx', 'grid_ny', 'domain_lx', 'domain_ly', 'reynolds_number', 'reynolds_number', 'reynolds_number', 'perturbation_amplitude', 'perturbation_mode');
var i, k: Integer; found, have: Boolean; u, msg: string;
begin
  have := False; u := '';
  for i := 0 to jo.Count - 1 do begin
    found := False;
    for k := Low(KNOWN) to High(KNOWN) do
      if jo.Names[i] = KNOWN[k] then begin found := True; Break; end;
    if (not found) and ((not have) or (jo.Names[i] < u)) then begin
      u := jo.Names[i]; have := True;
    end;
  end;
  if not have then Exit;
  msg := 'unknown key ''' + u + ''' (not part of the request contract)';
  for k := Low(FOREIGN_FROM) to High(FOREIGN_FROM) do
    if u = FOREIGN_FROM[k] then begin
      msg := 'unknown key ''' + u + ''' (did you mean ''' + FOREIGN_TO[k] + '''?)';
      Break;
    end;
  WriteLn(StdErr, 'KH-VAL-2 K7: unknown key ''', u, ''' refused');
  { StdErr is a buffered text file when it is a pipe (docker logs): flush it. }
  Flush(StdErr);
  raise EKhInvalid.CreateField(u, msg);
end;

function ParseSimRequest(const body: string): TSimRequest;
var
  jd:   TJSONData;
  jo:   TJSONObject;
  jarr: TJSONData;
  ja:   TJSONArray;
  k, n: Integer;
begin
  jd := GetJSON(body);
  if not (jd is TJSONObject) then begin
    jd.Free;
    raise EKhInvalid.CreateField('body', 'request body must be a JSON object');
  end;
  jo := TJSONObject(jd);
  try
    CheckUnknownKeys(jo);
    { Defaults = openapi.yaml = the oracle (128 x 64; was 64 x 32 until 2026-10-09). }
    Result.nx    := JGetI(jo, 'grid_nx',               128);
    Result.ny    := JGetI(jo, 'grid_ny',               64);
    Result.steps := JGetI(jo, 'steps',                 100);
    Result.mode  := JGetI(jo, 'perturbation_mode',     2);
    Result.lx    := JGetD(jo, 'domain_lx',             1.0);
    Result.ly    := JGetD(jo, 'domain_ly',             0.5);
    Result.dt    := JGetD(jo, 'dt',                    0.001);
    Result.re    := JGetD(jo, 'reynolds_number',       1000.0);
    Result.u0    := JGetD(jo, 'velocity_shear',        1.0);
    Result.amp   := JGetD(jo, 'perturbation_amplitude', 0.01);
    { Phase 7f forcing API (defaults reproduce pre-7f behaviour exactly) }
    Result.forcing_type      := JGetS(jo, 'forcing_type',      'none');
    Result.forcing_amplitude := JGetD(jo, 'forcing_amplitude', 0.0);
    Result.forcing_kx        := JGetI(jo, 'forcing_kx',        1);
    Result.forcing_ky        := JGetI(jo, 'forcing_ky',        1);
    Result.forcing_kmin      := JGetD(jo, 'forcing_kmin',      0.0);
    Result.forcing_kmax      := JGetD(jo, 'forcing_kmax',      0.0);
    Result.forcing_omega0    := JGetD(jo, 'forcing_omega0',    0.0);
    Result.linear_drag       := JGetD(jo, 'linear_drag',       0.0);
    Result.t_start           := JGetD(jo, 't_start',           0.0);
    { Optional custom IC (flattened row-major, length nx*ny).
      Its length is checked in ValidateRequest. Until 2026-09-26 an array of
      the wrong length was silently ignored and the canonical IC computed
      instead - a dropped parameter, the same class of defect as finding A. }
    Result.has_init_omega := False;
    SetLength(Result.init_omega, 0);
    { The canonical name is initial_omega (shared/physics/kh_physics.py and the
      C++ backend). This backend historically read only init_omega, so a client
      using the canonical name had its initial condition silently ignored.
      Both are accepted; sending both is ambiguous and refused.              }
    if Assigned(JPresent(jo, 'initial_omega')) and Assigned(JPresent(jo, 'init_omega')) then
      raise EKhInvalid.CreateField('initial_omega',
        'send initial_omega or its legacy alias init_omega, not both');
    jarr := JPresent(jo, 'initial_omega');
    if not Assigned(jarr) then jarr := JPresent(jo, 'init_omega');
    { A value that is not an array was ignored until 2026-10-03 (P-5). }
    if Assigned(jarr) and (jarr.JSONType <> jtArray) then
      raise EKhInvalid.CreateField('initial_omega', 'initial_omega must be an array of numbers');
    if Assigned(jarr) then begin
      ja := jarr as TJSONArray;
      n  := ja.Count;
      if n > MAX_CELLS then
        raise EKhInvalid.CreateField('initial_omega',
          Format('initial_omega has %d entries, more than the %d-cell limit', [n, MAX_CELLS]));
      SetLength(Result.init_omega, n);
      for k := 0 to n-1 do begin
        if ja.Items[k].JSONType <> jtNumber then
          raise EKhInvalid.CreateField('initial_omega', 'initial_omega must contain numbers only');
        Result.init_omega[k] := ja.Items[k].AsFloat;
      end;
      Result.has_init_omega := True;
      { K4b: an amplitude that was not sent cannot conflict with initial_omega. It has
        no effect on the computation either way, so 0 records "not sent". }
      if not Assigned(JPresent(jo, 'perturbation_amplitude')) then Result.amp := 0.0;
    end;
  finally
    jo.Free;
  end;
end;

{ -- Request validation ----------------------------------------------------- }
function IsPowerOfTwo(v: Integer): Boolean;
begin
  Result := (v >= 1) and ((v and (v - 1)) = 0);
end;

procedure CheckGridAxis(const field: string; v: Integer);
begin
  if (v < 1) or (v > MAX_AXIS) then
    raise EKhInvalid.CreateField(field,
      Format('%s = %d is out of range [1, %d]', [field, v, MAX_AXIS]));
  { The FFT is radix-2 Cooley-Tukey (kh_physics.FFT1D). }
  if not IsPowerOfTwo(v) then
    raise EKhInvalid.CreateField(field,
      Format('%s = %d is not a power of two; this backend uses a radix-2 FFT', [field, v]));
end;

{ Rejects exactly the inputs that cannot be computed: a grid the FFT cannot
  transform, a grid too large to allocate safely, and values the solver
  divides by. Nothing physically questionable-but-computable is refused here. }
procedure ValidateRequest(const req: TSimRequest);
var
  kfl, kr: Double;
  nmin, km, ix, iy: Integer;
  found: Boolean;
begin
  CheckGridAxis('grid_nx', req.nx);
  CheckGridAxis('grid_ny', req.ny);
  if Int64(req.nx) * Int64(req.ny) > MAX_CELLS then
    raise EKhInvalid.CreateField('grid_nx',
      Format('grid %d x %d = %d cells exceeds the %d-cell limit',
             [req.nx, req.ny, Int64(req.nx) * Int64(req.ny), MAX_CELLS]));
  { dx = lx/nx, dy = ly/ny and nu = u0/re are divisors in the solver. }
  if req.lx = 0.0 then
    raise EKhInvalid.CreateField('domain_lx', 'domain_lx must be non-zero');
  if req.ly = 0.0 then
    raise EKhInvalid.CreateField('domain_ly', 'domain_ly must be non-zero');
  if req.re = 0.0 then
    raise EKhInvalid.CreateField('reynolds_number', 'reynolds_number must be non-zero');
  if req.has_init_omega and (Length(req.init_omega) <> req.nx * req.ny) then
    raise EKhInvalid.CreateField('initial_omega',
      Format('initial_omega has %d entries; grid_nx * grid_ny = %d',
             [Length(req.init_omega), req.nx * req.ny]));
  { KH-VAL-2, K3: outside the model. }
  if req.lx < 0.0 then raise EKhInvalid.CreateField('domain_lx', 'domain_lx must be positive');
  if req.ly < 0.0 then raise EKhInvalid.CreateField('domain_ly', 'domain_ly must be positive');
  if req.re < 0.0 then
    raise EKhInvalid.CreateField('reynolds_number',
      'reynolds_number must be positive (a negative value is negative viscosity)');
  if req.dt < 0.0 then
    raise EKhInvalid.CreateField('dt', 'dt must be >= 0 (backward integration of a viscous flow is ill-posed)');
  if req.steps < 0 then raise EKhInvalid.CreateField('steps', 'steps must be >= 0');
  if req.linear_drag < 0.0 then
    raise EKhInvalid.CreateField('linear_drag', 'linear_drag must be >= 0 (negative drag injects energy)');
  if req.mode < 0 then raise EKhInvalid.CreateField('perturbation_mode', 'perturbation_mode must be >= 0');
  if (req.forcing_type = 'band') and (req.forcing_amplitude <> 0.0)
     and (req.forcing_kmin > req.forcing_kmax) then
    raise EKhInvalid.CreateField('forcing_kmin', 'forcing_kmin > forcing_kmax: the forcing band is empty');
  { KH-VAL-2, K4 (forcing): the forced modes must lie inside the 2/3-rule band. }
  if (req.forcing_type = 'mode') and (req.forcing_amplitude <> 0.0) then begin
    if 3 * Abs(Int64(req.forcing_kx)) > 2 * Int64(req.nx div 2) then
      raise EKhInvalid.CreateField('forcing_kx',
        Format('forcing_kx %d lies beyond the 2/3-rule band of grid_nx %d', [req.forcing_kx, req.nx]));
    if 3 * Abs(Int64(req.forcing_ky)) > 2 * Int64(req.ny div 2) then
      raise EKhInvalid.CreateField('forcing_ky',
        Format('forcing_ky %d lies beyond the 2/3-rule band of grid_ny %d', [req.forcing_ky, req.ny]));
  end;
  if (req.forcing_type = 'band') and (req.forcing_amplitude <> 0.0) then begin
    kfl := Floor(req.forcing_kmax);
    nmin := req.nx;  if req.ny < nmin then nmin := req.ny;
    if 3.0 * kfl > 2.0 * (nmin div 2) then
      raise EKhInvalid.CreateField('forcing_kmax',
        Format('forcing_kmax reaches beyond the 2/3-rule band of the %d x %d grid', [req.nx, req.ny]));
    { K3: the band must contain an integer mode (half-plane, kmin <= |k| <= kmax). }
    km := Trunc(kfl);
    found := False;
    for ix := 0 to km do begin
      for iy := -km to km do begin
        if (ix = 0) and (iy <= 0) then Continue;
        kr := Sqrt(Double(ix * ix + iy * iy));   { not 1.0*: FPC types that literal as Single }
        if (kr >= req.forcing_kmin) and (kr <= req.forcing_kmax) then begin found := True; Break; end;
      end;
      if found then Break;
    end;
    if not found then
      raise EKhInvalid.CreateField('forcing_kmax', 'the forcing band contains no integer mode: nothing would be forced');
  end;
  { KH-VAL-2, K4: the analytic perturbation must exist on the grid, inside the 2/3-rule band. }
  if (not req.has_init_omega) and (req.amp <> 0.0) then begin
    if req.mode = 0 then
      raise EKhInvalid.CreateField('perturbation_mode',
        'perturbation_mode 0 with a non-zero amplitude perturbs nothing');
    if 3 * Int64(req.mode) > 2 * Int64(req.nx div 2) then
      raise EKhInvalid.CreateField('perturbation_mode',
        Format('perturbation_mode %d lies beyond the 2/3-rule band of grid_nx %d (largest resolved mode %d)',
               [req.mode, req.nx, 2 * (req.nx div 2) div 3]));
  end;
  { KH-VAL-2, K4b: initial_omega replaces the analytic field; an amplitude sent with it would be ignored. }
  if req.has_init_omega and (req.amp <> 0.0) then
    raise EKhInvalid.CreateField('perturbation_amplitude',
      'perturbation_amplitude has no effect with initial_omega (send 0 or omit it)');
end;

{ KH-VAL-2, K5: computed, but possibly not resolved. Same conditions, order and
  prefixes in every backend and the oracle (kh_physics.request_warnings). }
function KhWarnings(const req: TSimRequest; cfl_peak: Double): TJSONArray;
var nu, dmin, diff: Double;
begin
  Result := TJSONArray.Create;
  nu   := Abs(req.u0) / req.re;
  dmin := req.lx / req.nx;
  if req.ly / req.ny < dmin then dmin := req.ly / req.ny;
  diff := 2.0 * nu * req.dt / (dmin * dmin);
  if cfl_peak > 1.0 then
    Result.Add(Format('cfl_advective: peak Courant number %.3g exceeds 1', [cfl_peak]));
  if diff > 1.0 then
    Result.Add(Format('cfl_diffusive: 2*nu*dt/min(dx,dy)^2 = %.3g exceeds 1', [diff]));
  if nu = 0.0 then
    Result.Add('inviscid: nu = 0, no physical dissipation; enstrophy accumulates at the 2/3 cutoff');
end;

{ -- Response building ------------------------------------------------------ }
{ B-2 (2026-09-26): the response states the configuration it computed, so a
  caller can prove the request landed. A response that does not say what it
  computed can only be checked by assumption - that is how finding A stayed
  hidden. steps_completed and t_final come from the solver, not the request. }
function BuildSimResponse(const req: TSimRequest; const res: TSimResult): string;
var jo, diag: TJSONObject;
begin
  jo := TJSONObject.Create;
  try
    diag := TJSONObject.Create;
    diag.Add('kinetic_energy',  res.ke);
    diag.Add('enstrophy',       res.enstrophy);
    diag.Add('max_vorticity',   res.max_vort);
    diag.Add('divergence_rms',  res.div_rms);
    diag.Add('cfl_peak',        res.cfl_peak);
    jo.Add('diagnostics', diag);
    jo.Add('backend',         SERVICE_ID);
    jo.Add('language',        'Pascal');
    jo.Add('grid_nx',         req.nx);
    jo.Add('grid_ny',         req.ny);
    jo.Add('steps_completed', res.steps_completed);
    jo.Add('t_final',         res.t_final);
    jo.Add('compute_time_ms', res.compute_ms);
    { P-6 (2026-10-03): the duplicate key compute_ms, kept on 2026-09-26 "for
      existing consumers", is gone. No consumer reads it from this response
      (log-service's compute_ms is its own ingest field), and it made the
      top-level keys differ from the other four backends, which
      tests/integration/test_cross_backend.py forbids.                      }
    JAddArr(jo, 'u_velocity', res.u);
    JAddArr(jo, 'v_velocity', res.v);
    JAddArr(jo, 'vorticity',  res.omega);
    JAddArr(jo, 'pressure',   res.psi);
    jo.Add('warnings', KhWarnings(req, res.cfl_peak));
    Result := jo.AsJSON;
  finally
    jo.Free;
  end;
end;

{ -- Handler object - owns routing; method bound to Server.OnRequest -------- }
type
  TKhHandler = class(TObject)
  private
    procedure DoSimulate(ARequest: TFPHTTPConnectionRequest;
                         AResponse: TFPHTTPConnectionResponse);
    procedure DoHealth  (AResponse: TFPHTTPConnectionResponse);
    procedure DoInfo    (AResponse: TFPHTTPConnectionResponse);
  public
    procedure HandleRequest(Sender: TObject;
                            var ARequest: TFPHTTPConnectionRequest;
                            var AResponse: TFPHTTPConnectionResponse);
  end;

procedure TKhHandler.DoSimulate(ARequest: TFPHTTPConnectionRequest;
                                AResponse: TFPHTTPConnectionResponse);
var req: TSimRequest; res: TSimResult;
begin
  { 1. parse - anything unparseable is the caller's error }
  try
    req := ParseSimRequest(ARequest.Content);
  except
    on E: EKhInvalid do begin
      AResponse.Code := 400;  AResponse.Content := ErrorJSON(E.Message, 'field', E.Field);  Exit;
    end;
    on E: Exception do begin
      AResponse.Code := 400;
      AResponse.Content := ErrorJSON('malformed request: ' + E.Message, 'field', 'body');
      Exit;
    end;
  end;

  { 2. validate - before any allocation or FFT touches the grid }
  try
    ValidateRequest(req);
  except
    on E: EKhInvalid do begin
      AResponse.Code := 400;  AResponse.Content := ErrorJSON(E.Message, 'field', E.Field);  Exit;
    end;
  end;

  { 3. compute - a request error raised inside the solver (forcing options)
       is still 400; any other failure is the computation's, so 500.     }
  try
    res := kh_physics.Simulate(req);
  except
    on E: EKhInvalid do begin
      AResponse.Code := 400;  AResponse.Content := ErrorJSON(E.Message, 'field', E.Field);  Exit;
    end;
    on E: Exception do begin
      AResponse.Code := 500;
      AResponse.Content := ErrorJSON(E.ClassName + ': ' + E.Message, 'backend', BACKEND_ID);
      Exit;
    end;
  end;

  AResponse.Code    := 200;
  AResponse.Content := BuildSimResponse(req, res);
end;

procedure TKhHandler.DoHealth(AResponse: TFPHTTPConnectionResponse);
begin
  AResponse.Code    := 200;
  AResponse.Content := '{"status":"ok","backend":"' + SERVICE_ID + '","port":8005}';
end;

{ B-1 (2026-09-26): the keys tests/integration/test_info.py requires, plus
  provenance. build_date is fixed at compile time, so a container whose image
  predates the current source is visible from outside - which is how finding A
  would have been caught in March rather than September.                     }
procedure TKhHandler.DoInfo(AResponse: TFPHTTPConnectionResponse);
var jo: TJSONObject;
begin
  jo := TJSONObject.Create;
  try
    jo.Add('backend',          SERVICE_ID);
    jo.Add('language',         'Pascal');
    jo.Add('language_version', {$I %FPCVERSION%});
    jo.Add('compiler',         'fpc ' + {$I %FPCVERSION%});
    jo.Add('framework',        'fphttpserver');
    jo.Add('fft_library',      'builtin-cooley-tukey-radix2');
    jo.Add('grid_constraint',  'power-of-two');
    jo.Add('openapi_spec',     'kh-sim/shared/api/openapi.yaml');
    jo.Add('build_date',       StringReplace({$I %DATE%}, '/', '-', [rfReplaceAll]));
    jo.Add('port',             SERVER_PORT);
    AResponse.Code    := 200;
    AResponse.Content := jo.AsJSON;
  finally
    jo.Free;
  end;
end;

procedure TKhHandler.HandleRequest(Sender: TObject;
                                   var ARequest: TFPHTTPConnectionRequest;
                                   var AResponse: TFPHTTPConnectionResponse);
var path: string; qi: Integer;
begin
  AResponse.ContentType := 'application/json';
  { Strip query string from URL }
  path := ARequest.URL;
  qi   := Pos('?', path);
  if qi > 0 then path := Copy(path, 1, qi-1);
  if path = '' then path := '/';

  if (path = '/simulate') and (ARequest.Method = 'POST') then
    DoSimulate(ARequest, AResponse)
  else if path = '/health' then
    DoHealth(AResponse)
  else if path = '/info' then
    DoInfo(AResponse)
  else begin
    AResponse.Code    := 404;
    AResponse.Content := '{"error":"not found"}';
  end;
end;

{ -- Entry point ------------------------------------------------------------ }
var
  Handler: TKhHandler;
  Server:  TFPHttpServer;

begin
  Handler := TKhHandler.Create;
  Server  := TFPHttpServer.Create(nil);
  try
    Server.Port      := SERVER_PORT;
    Server.Threaded  := False;          { single-threaded accept loop }
    Server.OnRequest := @Handler.HandleRequest;
    WriteLn('kh-sim Pascal backend listening on http://0.0.0.0:8005');
    WriteLn('Press Ctrl+C to stop.');
    Server.Active := True;              { blocks - runs accept loop }
  finally
    Server.Free;
    Handler.Free;
  end;
end.
