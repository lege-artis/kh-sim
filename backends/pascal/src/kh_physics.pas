unit kh_physics;
{ KH instability physics kernel - Pascal (Free Pascal 3.2+)
  Physics : 2D incompressible Navier-Stokes, vorticity-streamfunction form.
  Numerics : pseudo-spectral (Cooley-Tukey radix-2), RK4.
  Reference: kh-sim/shared/physics/KH-PHYSICS.md
  Canonical: kh-sim/shared/physics/kh_physics.py
}
{$mode objfpc}{$H+}

interface

uses Math, SysUtils, DateUtils;

{ -- Complex arithmetic ----------------------------------------------------- }
type
  TComplex    = record re, im: Double; end;
  TComplexArr = array of TComplex;
  TRealArr    = array of Double;

operator + (const a, b: TComplex): TComplex; inline;
operator - (const a, b: TComplex): TComplex; inline;
operator * (const a, b: TComplex): TComplex; inline;
operator * (const a: Double; const b: TComplex): TComplex; inline;
operator * (const a: TComplex; const b: Double): TComplex; inline;
operator / (const a: TComplex; const b: Double): TComplex; inline;

function CmplxOf(re, im: Double): TComplex; inline;
function CmplxFromReal(const r: TRealArr; n: Integer): TComplexArr;
function RealPart(const c: TComplexArr; n: Integer): TRealArr;

{ -- Simulation API --------------------------------------------------------- }
type
  TSimRequest = record
    nx, ny, steps, mode: Integer;
    lx, ly, dt, re, u0, amp: Double;
    { Optional custom IC: set has_init_omega=True and populate init_omega
      (flattened row-major, length nx*ny) to override canonical tanh IC.
      Added for phase 7d physics testbeds (TB-1a/1b/3a).                 }
    has_init_omega: Boolean;
    init_omega:     TRealArr;
    { Phase 7f forcing API (PHASE-7F-FORCING-API-DESIGN-v0.1.md S1 v0.1 + S7 v0.2).
      All optional; defaults reproduce pre-7f behaviour exactly (forcing_type=''
      or ''none'' or amplitude 0 => the RHS is bit-identical to the pre-7f path,
      so phase-6c golden vectors + phase-7d testbeds are unaffected).           }
    forcing_type:      String;   { ''none'' | ''mode'' | ''band'' (empty => none) }
    forcing_amplitude: Double;   { A }
    forcing_kx:        Integer;  { ''mode'' only: integer mode index along x }
    forcing_ky:        Integer;  { ''mode'' only: integer mode index along y }
    forcing_kmin:      Double;   { ''band'' only: inner integer-mode radius }
    forcing_kmax:      Double;   { ''band'' only: outer integer-mode radius }
    forcing_omega0:    Double;   { v0.2: phase-rotation rate scale (band only) }
    linear_drag:       Double;   { v0.2: Ekman-type linear drag alpha (-alpha*omega) }
    { Time of the initial field (2026-10-10, CP-7F-B/D harness): the forcing runs
      from t_start and t_final = t_start + steps*dt, so a continued run keeps its
      phase. Every caller sets it (a local record is not initialised). }
    t_start:           Double;
  end;

  TSimResult = record
    u, v, omega, psi: TRealArr;
    ke, enstrophy, max_vort, div_rms, cfl_peak, compute_ms: Double;
    { B-2 (2026-09-26): what the solver actually did, reported by the solver
      rather than echoed from the request. A 200 response implies every step
      ran, because a non-finite value raises an FPU exception first.        }
    steps_completed: Integer;
    t_final:         Double;
  end;

  { Raised for a request that cannot be computed as given. The HTTP layer maps
    it to 400 with a "field"; any other exception during Simulate maps to 500.
    See _audit/FINDING-PASCAL-HTTP-PARAMETER-DROP-2026-09-16.md, B-3.        }
  EKhInvalid = class(Exception)
  public
    Field: String;
    constructor CreateField(const AField, AMsg: String);
  end;

{ The radix-2 Cooley-Tukey transform the solver uses (n a power of two;
  inverse=True applies 1/n). Exported 2026-10-11 for CP-7F-E/E1, which runs the
  NR convlv example through it (tests/kh_cp7fe.pas); the solver is unchanged. }
procedure FFT1D(var a: TComplexArr; n: Integer; inverse: Boolean);
procedure DealiasField(var field_hat: TComplexArr; const kx, ky: TRealArr; nx, ny: Integer);
function CflDtBound(u_max, v_max, dx, dy, nu: Double): Double;
function Simulate(const req: TSimRequest): TSimResult;

implementation

constructor EKhInvalid.CreateField(const AField, AMsg: String);
begin
  inherited Create(AMsg);
  Field := AField;
end;

{ -- Operator implementations ----------------------------------------------- }
operator + (const a, b: TComplex): TComplex;
begin Result.re := a.re + b.re; Result.im := a.im + b.im; end;

operator - (const a, b: TComplex): TComplex;
begin Result.re := a.re - b.re; Result.im := a.im - b.im; end;

operator * (const a, b: TComplex): TComplex;
begin
  Result.re := a.re*b.re - a.im*b.im;
  Result.im := a.re*b.im + a.im*b.re;
end;

operator * (const a: Double; const b: TComplex): TComplex;
begin Result.re := a*b.re; Result.im := a*b.im; end;

operator * (const a: TComplex; const b: Double): TComplex;
begin Result.re := a.re*b; Result.im := a.im*b; end;

operator / (const a: TComplex; const b: Double): TComplex;
begin Result.re := a.re/b; Result.im := a.im/b; end;

function CmplxOf(re, im: Double): TComplex;
begin Result.re := re; Result.im := im; end;

function CmplxFromReal(const r: TRealArr; n: Integer): TComplexArr;
var i: Integer;
begin
  SetLength(Result, n);
  for i := 0 to n-1 do begin Result[i].re := r[i]; Result[i].im := 0.0; end;
end;

function RealPart(const c: TComplexArr; n: Integer): TRealArr;
var i: Integer;
begin
  SetLength(Result, n);
  for i := 0 to n-1 do Result[i] := c[i].re;
end;

{ -- 1D in-place Cooley-Tukey FFT (radix-2, DIT) --------------------------- }
{ Iterative DIT Cooley-Tukey FFT.
  Reference: Press, Flannery, Teukolsky & Vetterling (1986),
             Numerical Recipes in Pascal, S12.2, eq. 12.2.2
             (decimation-in-time butterfly recurrence, p.394 of the
             Pascal edition); routine FOUR1 is the textbook analogue.
  Sign convention: forward = exp(-2*pi*i*k*n/N) per S12.1 (standard);
             inverse applies 1/N normalisation. Phase-6c sign-convention
             fix is commit e77379b (conjugate-of-standard bug corrected).
  Complexity: O(N log N) -- log2(N) stages of N/2 butterflies, per S12.2;
             verified empirically by the phase-7d/7f timings.
  Doctrine 4.4 C3 citation header (Pascal canonical-reference backend). }
{ n must be a power of 2. inverse=True applies 1/n normalisation.             }
procedure FFT1D(var a: TComplexArr; n: Integer; inverse: Boolean);
var
  i, j, bit, len, k: Integer;
  ang: Double;
  wlen, w, u, v: TComplex;
begin
  if (n < 1) or ((n and (n - 1)) <> 0) then
    raise Exception.CreateFmt('FFT1D: length %d is not a power of two', [n]);
  { bit-reversal permutation }
  j := 0;
  for i := 1 to n-1 do begin
    bit := n shr 1;
    while (j and bit) <> 0 do begin
      j := j xor bit;
      bit := bit shr 1;
    end;
    j := j xor bit;
    if i < j then begin u := a[i]; a[i] := a[j]; a[j] := u; end;
  end;
  { butterfly stages }
  len := 2;
  while len <= n do begin
    ang := -2.0*Pi / len;
    if inverse then ang := -ang;
    wlen := CmplxOf(Cos(ang), Sin(ang));
    i := 0;
    while i < n do begin
      w := CmplxOf(1.0, 0.0);
      for k := 0 to (len shr 1)-1 do begin
        u := a[i+k];
        v := a[i+k+(len shr 1)] * w;
        a[i+k]              := u + v;
        a[i+k+(len shr 1)] := u - v;
        w := w * wlen;
      end;
      i := i + len;
    end;
    len := len shl 1;
  end;
  { normalise inverse }
  if inverse then
    for i := 0 to n-1 do a[i] := a[i] / n;
end;

{ -- 2D FFT: row-then-column (row-major a[i*ny+j]) ------------------------- }
procedure FFT2D(var a: TComplexArr; nx, ny: Integer; inverse: Boolean);
var
  i, j: Integer;
  row, col: TComplexArr;
begin
  SetLength(row, ny);
  SetLength(col, nx);
  { FFT along y-axis (rows of length ny) }
  for i := 0 to nx-1 do begin
    for j := 0 to ny-1 do row[j] := a[i*ny+j];
    FFT1D(row, ny, inverse);
    for j := 0 to ny-1 do a[i*ny+j] := row[j];
  end;
  { FFT along x-axis (columns of length nx) }
  for j := 0 to ny-1 do begin
    for i := 0 to nx-1 do col[i] := a[i*ny+j];
    FFT1D(col, nx, inverse);
    for i := 0 to nx-1 do a[i*ny+j] := col[i];
  end;
end;

{ -- Angular wavenumber meshgrid, row-major --------------------------------- }
{ kx[i*ny+j] = 2pi*fftfreq(i,nx)/dx,  ky[i*ny+j] = 2pi*fftfreq(j,ny)/dy     }
procedure MakeWavenumbers(nx, ny: Integer; dx, dy: Double;
  out kx, ky: TRealArr);
var i, j, fi, fj: Integer;
begin
  SetLength(kx, nx*ny);
  SetLength(ky, nx*ny);
  for i := 0 to nx-1 do begin
    if i <= nx div 2 then fi := i else fi := i - nx;
    for j := 0 to ny-1 do kx[i*ny+j] := 2.0*Pi*fi / (nx*dx);
  end;
  for j := 0 to ny-1 do begin
    if j <= ny div 2 then fj := j else fj := j - ny;
    for i := 0 to nx-1 do ky[i*ny+j] := 2.0*Pi*fj / (ny*dy);
  end;
end;

{ -- 2/3-rule dealiasing (KH-NM-4, Orszag 1971) ---------------------------- }
{ Zeros spectral modes where |kx| > 2/3*kx_max or |ky| > 2/3*ky_max.        }
{ -- Phase 7f forcing operator (PHASE-7F-FORCING-API-DESIGN-v0.1.md S3 + S7) -- }
{ Time-independent ''mode'' forcing and (optionally rotating-phase) ''band''
  forcing, plus an Ekman-type linear drag term. The deterministic, RNG-free
  phase construction makes the forcing field bit-near-identical across all six
  backends (S3: phi0(ix,iy)=2*pi*frac(ix*sqrt2+iy*sqrt3); S7:
  Om(ix,iy)=2*pi*frac(ix*sqrt11+iy*sqrt13)*omega0). When forcing is off the
  whole machinery is skipped (active=False) and the RHS path is bit-identical
  to the pre-7f kernel.                                                        }
type
  TForcing = record
    active:  Boolean;
    isband:  Boolean;          { True => rotating-band; False => static field }
    fstatic: TRealArr;         { static (''mode'' or v0.1 band) field, length n }
    cmode:   array of TRealArr; { band: per-mode C_m spatial field }
    smode:   array of TRealArr; { band: per-mode S_m spatial field }
    omegam:  TRealArr;          { band: per-mode rotation rate Om_m }
    nmodes:  Integer;
    n:       Integer;
  end;

function FracPos(x: Double): Double; inline;
begin
  Result := x - Floor(x);
  if Result < 0.0 then Result := Result + 1.0;
end;

{ Build the forcing operator once from the request. Returns active=False when
  forcing is off (type ''none''/'''' or amplitude 0) => pre-7f RHS unchanged.   }
function BuildForcing(const req: TSimRequest; nx, ny: Integer;
  lx, ly: Double): TForcing;
var
  ftype: String;
  amp, dx, dy, x, y, cx, cy, kr, scale, phi0, om, arg0: Double;
  kxm, kym, i, j, ix, iy, kmax_i, idx, m: Integer;
  modes_ix, modes_iy: array of Integer;
begin
  Result.active := False;
  Result.isband := False;
  Result.nmodes := 0;
  Result.n      := nx * ny;
  ftype := LowerCase(Trim(req.forcing_type));
  amp   := req.forcing_amplitude;
  if (ftype = '') or (ftype = 'none') or (amp = 0.0) then Exit;

  dx := lx / nx;  dy := ly / ny;

  if ftype = 'mode' then begin
    { F(x,y) = A * cos(2*pi*kx*x/lx) * cos(2*pi*ky*y/ly)  (steady) }
    kxm := req.forcing_kx;  kym := req.forcing_ky;
    SetLength(Result.fstatic, nx*ny);
    for i := 0 to nx-1 do begin
      x  := i * dx;
      cx := Cos(2.0*Pi*kxm*x/lx);
      for j := 0 to ny-1 do begin
        y  := j * dy;
        cy := Cos(2.0*Pi*kym*y/ly);
        Result.fstatic[i*ny+j] := amp * cx * cy;
      end;
    end;
    Result.active := True;
    Result.isband := False;
    Exit;
  end;

  if ftype = 'band' then begin
    { Isotropic shell of integer mode pairs (ix,iy) in the half-plane
      (ix >= 0; skip ix==0 and iy<=0) with kmin <= sqrt(ix^2+iy^2) <= kmax.
      Requires lx == ly (integer-mode-radius isotropy).                    }
    if Abs(lx - ly) > 1.0e-12 then
      raise EKhInvalid.CreateField('forcing_type', 'band forcing requires domain_lx == domain_ly');
    kmax_i := Trunc(req.forcing_kmax);
    SetLength(modes_ix, 0);  SetLength(modes_iy, 0);
    for ix := 0 to kmax_i do
      for iy := -kmax_i to kmax_i do begin
        if (ix = 0) and (iy <= 0) then Continue;
        kr := Sqrt(Double(ix*ix + iy*iy));
        if (kr >= req.forcing_kmin) and (kr <= req.forcing_kmax) then begin
          SetLength(modes_ix, Length(modes_ix)+1);
          SetLength(modes_iy, Length(modes_iy)+1);
          modes_ix[High(modes_ix)] := ix;
          modes_iy[High(modes_iy)] := iy;
        end;
      end;
    Result.nmodes := Length(modes_ix);
    if Result.nmodes = 0 then Exit;  { active stays False }
    scale := amp / Sqrt(Double(Result.nmodes));
    SetLength(Result.cmode,  Result.nmodes);
    SetLength(Result.smode,  Result.nmodes);
    SetLength(Result.omegam, Result.nmodes);
    for m := 0 to Result.nmodes-1 do begin
      ix := modes_ix[m];  iy := modes_iy[m];
      phi0 := 2.0*Pi*FracPos(ix*Sqrt(2.0) + iy*Sqrt(3.0));
      om   := 2.0*Pi*FracPos(ix*Sqrt(11.0) + iy*Sqrt(13.0)) * req.forcing_omega0;
      Result.omegam[m] := om;
      SetLength(Result.cmode[m], nx*ny);
      SetLength(Result.smode[m], nx*ny);
      for i := 0 to nx-1 do begin
        x := i * dx;
        for j := 0 to ny-1 do begin
          y := j * dy;
          arg0 := 2.0*Pi*(ix*x/lx + iy*y/ly) + phi0;
          idx  := i*ny+j;
          Result.cmode[m][idx] := scale * Cos(arg0);
          Result.smode[m][idx] := scale * Sin(arg0);
        end;
      end;
    end;
    Result.active := True;
    Result.isband := True;
    Exit;
  end;

  raise EKhInvalid.CreateField('forcing_type', 'unknown forcing_type: ' + req.forcing_type);
end;

{ Evaluate the forcing field at time t into F (length n). For ''mode''/v0.1-band
  the static field is returned (t-independent); for a rotating band:
  F(t) = sum_m [C_m*cos(Om_m t) - S_m*sin(Om_m t)].                            }
procedure EvalForcing(const fc: TForcing; t: Double; out F: TRealArr);
var idx, m: Integer; c, sn: Double;
begin
  SetLength(F, fc.n);
  if not fc.isband then begin
    for idx := 0 to fc.n-1 do F[idx] := fc.fstatic[idx];
    Exit;
  end;
  for idx := 0 to fc.n-1 do F[idx] := 0.0;
  for m := 0 to fc.nmodes-1 do begin
    c  := Cos(fc.omegam[m]*t);
    sn := Sin(fc.omegam[m]*t);
    for idx := 0 to fc.n-1 do
      F[idx] := F[idx] + fc.cmode[m][idx]*c - fc.smode[m][idx]*sn;
  end;
end;

procedure DealiasField(var field_hat: TComplexArr; const kx, ky: TRealArr; nx, ny: Integer);
var k, n: Integer;
    kx_max, ky_max, kx_cut, ky_cut: Double;
begin
  n := nx * ny;
  kx_max := 0.0;  ky_max := 0.0;
  for k := 0 to n-1 do begin
    if Abs(kx[k]) > kx_max then kx_max := Abs(kx[k]);
    if Abs(ky[k]) > ky_max then ky_max := Abs(ky[k]);
  end;
  kx_cut := (2.0 / 3.0) * kx_max;
  ky_cut := (2.0 / 3.0) * ky_max;
  for k := 0 to n-1 do
    if (Abs(kx[k]) > kx_cut) or (Abs(ky[k]) > ky_cut) then
      field_hat[k] := CmplxOf(0.0, 0.0);
end;

{ -- CFL stability bound (KH-NM-5) ----------------------------------------- }
{ Returns min of advective and diffusive CFL-limited time steps.              }
function CflDtBound(u_max, v_max, dx, dy, nu: Double): Double;
var dt_adv, dt_diff: Double;
begin
  dt_adv  := Min(dx / (u_max + 1e-30), dy / (v_max + 1e-30));
  { nu = 0 (inviscid) is a valid limit: no diffusive bound (modelling rules v0.1, D).
    Until 2026-10-07 this divided by zero and the FPU trap made the request a 500. }
  if nu = 0.0 then dt_diff := MaxDouble
  else dt_diff := Min(dx*dx / (2.0*nu), dy*dy / (2.0*nu));
  Result  := Min(dt_adv, dt_diff);
end;

{ -- Initial conditions: tanh shear + sinusoidal perturbation --------------- }
procedure InitConds(nx, ny: Integer; lx, ly, u0, delta, amp: Double;
  mode: Integer; out omega: TRealArr);
var i, j: Integer;
    dx, dy, x, y, z, pert: Double;
begin
  SetLength(omega, nx*ny);
  dx := lx / nx;  dy := ly / ny;
  for i := 0 to nx-1 do begin
    x    := i * dx;
    pert := amp * (2.0*Pi*mode/lx) * Cos(2.0*Pi*mode*x/lx);
    for j := 0 to ny-1 do begin
      y := j * dy;
      z := (y - ly*0.5) / delta;
      omega[i*ny+j] := (u0/delta) / (Cosh(z)*Cosh(z)) + pert;
    end;
  end;
end;

{ -- Spectral Poisson: psi_hat = omega_hat / k^2, zero mode = 0 ------------- }
procedure SolvePoisson(const oh: TComplexArr; const kx, ky: TRealArr;
  nx, ny: Integer; out ph: TComplexArr);
var i, j, idx: Integer; k2: Double;
begin
  SetLength(ph, nx*ny);
  for i := 0 to nx-1 do
    for j := 0 to ny-1 do begin
      idx := i*ny+j;
      if (i=0) and (j=0) then
        ph[idx] := CmplxOf(0.0, 0.0)
      else begin
        k2 := kx[idx]*kx[idx] + ky[idx]*ky[idx];
        ph[idx] := oh[idx] / k2;
      end;
    end;
end;

{ -- Velocity from streamfunction ------------------------------------------- }
{ u_hat = i*ky*psi_hat  ->  re=-ky*psi_hat.im, im=+ky*psi_hat.re                           }
{ v_hat = -i*kx*psi_hat ->  re=+kx*psi_hat.im, im=-kx*psi_hat.re                           }
procedure VelFromPsi(const ph: TComplexArr; const kx, ky: TRealArr;
  nx, ny: Integer; out u, v: TRealArr);
var i, j, idx, n: Integer;
    uh, vh: TComplexArr;
begin
  n := nx*ny;
  SetLength(uh, n); SetLength(vh, n);
  for i := 0 to nx-1 do
    for j := 0 to ny-1 do begin
      idx := i*ny+j;
      uh[idx] := CmplxOf(-ky[idx]*ph[idx].im,  ky[idx]*ph[idx].re);
      vh[idx] := CmplxOf( kx[idx]*ph[idx].im, -kx[idx]*ph[idx].re);
    end;
  FFT2D(uh, nx, ny, True);
  FFT2D(vh, nx, ny, True);
  u := RealPart(uh, n);
  v := RealPart(vh, n);
end;

{ -- Vorticity RHS: -u*domega/dx - v*domega/dy + nu*nabla^2omega -------------------------- }
{ Nonlinear term dealiased via 2/3-rule (KH-NM-4, Orszag 1971).             }
procedure VortRHS(const omega, uin, vin, kx, ky: TRealArr;
  nx, ny: Integer; nu: Double; drag: Double; use_force: Boolean;
  const force: TRealArr; out rhs: TRealArr);
var
  i, j, idx, n: Integer;
  k2: Double;
  oh, dxs, dys, laps, nl_hat: TComplexArr;
  dox, doy, lap, nl_phys: TRealArr;
begin
  n := nx*ny;
  oh := CmplxFromReal(omega, n);
  FFT2D(oh, nx, ny, False);
  SetLength(dxs, n); SetLength(dys, n); SetLength(laps, n);
  for i := 0 to nx-1 do
    for j := 0 to ny-1 do begin
      idx := i*ny+j;
      k2  := kx[idx]*kx[idx] + ky[idx]*ky[idx];
      { i*kx*omega_hat }  dxs[idx]  := CmplxOf(-kx[idx]*oh[idx].im,  kx[idx]*oh[idx].re);
      { i*ky*omega_hat }  dys[idx]  := CmplxOf(-ky[idx]*oh[idx].im,  ky[idx]*oh[idx].re);
      { -k^2*omega_hat  }  laps[idx] := oh[idx] * (-k2);
    end;
  FFT2D(dxs,  nx, ny, True);  dox := RealPart(dxs,  n);
  FFT2D(dys,  nx, ny, True);  doy := RealPart(dys,  n);
  FFT2D(laps, nx, ny, True);  lap := RealPart(laps, n);
  { Nonlinear product -> dealias -> IFFT (KH-NM-4) }
  SetLength(nl_hat, n);
  for idx := 0 to n-1 do
    nl_hat[idx] := CmplxOf(-uin[idx]*dox[idx] - vin[idx]*doy[idx], 0.0);
  FFT2D(nl_hat, nx, ny, False);
  DealiasField(nl_hat, kx, ky, nx, ny);
  FFT2D(nl_hat, nx, ny, True);
  nl_phys := RealPart(nl_hat, n);
  { Drag (-drag*omega) and forcing (+F) added per S2/S7; both no-ops when
    drag=0 and use_force=False => bit-identical to the pre-7f RHS.            }
  SetLength(rhs, n);
  for idx := 0 to n-1 do begin
    rhs[idx] := nl_phys[idx] + nu*lap[idx];
    if drag <> 0.0 then rhs[idx] := rhs[idx] - drag*omega[idx];
    if use_force   then rhs[idx] := rhs[idx] + force[idx];
  end;
end;

{ -- RK4 time step ---------------------------------------------------------- }
{ RK4 step with optional time-dependent forcing + linear drag (phase 7f).
  Forcing evaluated at RK4 stage times: k1@t, k2/k3@t+dt/2, k4@t+dt. When
  fc.active=False and drag=0 the four VortRHS calls take the pre-7f path
  (use_force=False) => bit-identical to the pre-7f kernel.                    }
procedure RK4Step(const omega, kx, ky: TRealArr;
  nx, ny: Integer; nu, dt, t, drag: Double; const fc: TForcing;
  out omega_out: TRealArr);
var
  n, idx: Integer;
  k1, k2, k3, k4, tmp, u, v: TRealArr;
  f1, f23, f4, fnull: TRealArr;
  oh, ph: TComplexArr;
  uf: Boolean;
begin
  n := nx*ny;
  SetLength(tmp, n);
  uf := fc.active;
  if uf then begin
    EvalForcing(fc, t,          f1);
    EvalForcing(fc, t + 0.5*dt, f23);
    EvalForcing(fc, t + dt,     f4);
  end else begin
    SetLength(fnull, 0); f1 := fnull; f23 := fnull; f4 := fnull;
  end;

  oh := CmplxFromReal(omega, n); FFT2D(oh, nx, ny, False);
  SolvePoisson(oh, kx, ky, nx, ny, ph);
  VelFromPsi(ph, kx, ky, nx, ny, u, v);
  VortRHS(omega, u, v, kx, ky, nx, ny, nu, drag, uf, f1, k1);

  for idx := 0 to n-1 do tmp[idx] := omega[idx] + 0.5*dt*k1[idx];
  oh := CmplxFromReal(tmp, n); FFT2D(oh, nx, ny, False);
  SolvePoisson(oh, kx, ky, nx, ny, ph);
  VelFromPsi(ph, kx, ky, nx, ny, u, v);
  VortRHS(tmp, u, v, kx, ky, nx, ny, nu, drag, uf, f23, k2);

  for idx := 0 to n-1 do tmp[idx] := omega[idx] + 0.5*dt*k2[idx];
  oh := CmplxFromReal(tmp, n); FFT2D(oh, nx, ny, False);
  SolvePoisson(oh, kx, ky, nx, ny, ph);
  VelFromPsi(ph, kx, ky, nx, ny, u, v);
  VortRHS(tmp, u, v, kx, ky, nx, ny, nu, drag, uf, f23, k3);

  for idx := 0 to n-1 do tmp[idx] := omega[idx] + dt*k3[idx];
  oh := CmplxFromReal(tmp, n); FFT2D(oh, nx, ny, False);
  SolvePoisson(oh, kx, ky, nx, ny, ph);
  VelFromPsi(ph, kx, ky, nx, ny, u, v);
  VortRHS(tmp, u, v, kx, ky, nx, ny, nu, drag, uf, f4, k4);

  SetLength(omega_out, n);
  for idx := 0 to n-1 do
    omega_out[idx] := omega[idx] + (dt/6.0)*(k1[idx] + 2.0*k2[idx] + 2.0*k3[idx] + k4[idx]);
end;

{ -- Main simulate entry point ---------------------------------------------- }
function Simulate(const req: TSimRequest): TSimResult;
var
  nx, ny, n, s, idx, k_cfl: Integer;
  lx, ly, dt, nu, delta, dx, dy: Double;
  omega, omega_next, kx, ky, u, v, psi_r: TRealArr;
  oh, ph, ph2, uh, vh, div_s: TComplexArr;
  ke_sum, ens_sum, max_w, div_sum, n_inv: Double;
  t_start: TDateTime;
  cfl_oh, cfl_ph: TComplexArr;
  cfl_u, cfl_v: TRealArr;
  cfl_umax, cfl_vmax, dt_bound: Double;
  fc: TForcing;        { phase 7f forcing operator (active=False => pre-7f path) }
  t_cur: Double;       { running simulation time for stage-time forcing eval }
begin
  nx := req.nx;  ny := req.ny;  n := nx*ny;
  lx := req.lx;  ly := req.ly;  dt := req.dt;
  { nu = |U| / Re: a reversed shear is the mirrored flow, not negative viscosity
    (modelling rules v0.1, K3). Identical for U >= 0. }
  nu    := Abs(req.u0) / req.re;
  delta := 0.05 * ly;
  dx    := lx / nx;
  dy    := ly / ny;

  MakeWavenumbers(nx, ny, dx, dy, kx, ky);
  if req.has_init_omega and (Length(req.init_omega) = n) then
    omega := req.init_omega
  else
    InitConds(nx, ny, lx, ly, req.u0, delta, req.amp, req.mode, omega);

  { Phase 7f: build the forcing operator once (active=False => pre-7f RHS). }
  fc := BuildForcing(req, nx, ny, lx, ly);

  { CFL guard (KH-NM-5): warn once if dt exceeds stability bound }
  cfl_oh := CmplxFromReal(omega, n);  FFT2D(cfl_oh, nx, ny, False);
  SolvePoisson(cfl_oh, kx, ky, nx, ny, cfl_ph);
  VelFromPsi(cfl_ph, kx, ky, nx, ny, cfl_u, cfl_v);
  cfl_umax := 0.0;  cfl_vmax := 0.0;
  for k_cfl := 0 to n-1 do begin
    if Abs(cfl_u[k_cfl]) > cfl_umax then cfl_umax := Abs(cfl_u[k_cfl]);
    if Abs(cfl_v[k_cfl]) > cfl_vmax then cfl_vmax := Abs(cfl_v[k_cfl]);
  end;
  dt_bound := CflDtBound(cfl_umax, cfl_vmax, dx, dy, nu);
  if dt > dt_bound then begin
    WriteLn(ErrOutput, 'KH-NM-5: dt=', dt:12:6, ' exceeds CFL bound ', dt_bound:12:6,
            '; simulation may be unstable');
    Flush(ErrOutput);   { buffered when redirected: the line arrives whole and on time }
  end;

  t_start := Now;
  Result.cfl_peak := 0.0;
  Result.steps_completed := 0;
  t_cur := req.t_start;
  for s := 1 to req.steps do begin
    { pre-step: recover velocity from current omega -> track CFL peak }
    cfl_oh := CmplxFromReal(omega, n);  FFT2D(cfl_oh, nx, ny, False);
    SolvePoisson(cfl_oh, kx, ky, nx, ny, cfl_ph);
    VelFromPsi(cfl_ph, kx, ky, nx, ny, cfl_u, cfl_v);
    cfl_umax := 0.0;  cfl_vmax := 0.0;
    for k_cfl := 0 to n-1 do begin
      if Abs(cfl_u[k_cfl]) > cfl_umax then cfl_umax := Abs(cfl_u[k_cfl]);
      if Abs(cfl_v[k_cfl]) > cfl_vmax then cfl_vmax := Abs(cfl_v[k_cfl]);
    end;
    dt_bound := dt * (cfl_umax / dx + cfl_vmax / dy);
    if dt_bound > Result.cfl_peak then Result.cfl_peak := dt_bound;
    RK4Step(omega, kx, ky, nx, ny, nu, dt, t_cur, req.linear_drag, fc, omega_next);
    omega := omega_next;
    t_cur := t_cur + dt;
    Inc(Result.steps_completed);
  end;
  Result.t_final := t_cur;
  Result.compute_ms := MilliSecondsBetween(Now, t_start);

  { final psi via ifft(solve_poisson(fft(omega))) }
  oh := CmplxFromReal(omega, n);  FFT2D(oh, nx, ny, False);
  SolvePoisson(oh, kx, ky, nx, ny, ph);
  SetLength(ph2, n);
  for idx := 0 to n-1 do ph2[idx] := ph[idx];
  FFT2D(ph2, nx, ny, True);
  psi_r := RealPart(ph2, n);

  { re-solve ph for velocity (ph was consumed by ifft above) }
  oh := CmplxFromReal(omega, n);  FFT2D(oh, nx, ny, False);
  SolvePoisson(oh, kx, ky, nx, ny, ph);
  VelFromPsi(ph, kx, ky, nx, ny, u, v);

  { diagnostics }
  n_inv := Double(1.0) / n;   { Double(): FPC types the literal 1.0 as Single }
  ke_sum := 0.0;  ens_sum := 0.0;  max_w := 0.0;
  for idx := 0 to n-1 do begin
    ke_sum  := ke_sum  + 0.5*(u[idx]*u[idx] + v[idx]*v[idx]);
    ens_sum := ens_sum + 0.5*omega[idx]*omega[idx];
    if Abs(omega[idx]) > max_w then max_w := Abs(omega[idx]);
  end;

  { divergence rms: div = i*kx*u_hat + i*ky*v_hat }
  uh := CmplxFromReal(u, n);  FFT2D(uh, nx, ny, False);
  vh := CmplxFromReal(v, n);  FFT2D(vh, nx, ny, False);
  SetLength(div_s, n);
  for idx := 0 to n-1 do
    div_s[idx] := CmplxOf(
      -kx[idx]*uh[idx].im - ky[idx]*vh[idx].im,
       kx[idx]*uh[idx].re + ky[idx]*vh[idx].re);
  FFT2D(div_s, nx, ny, True);
  div_sum := 0.0;
  for idx := 0 to n-1 do div_sum := div_sum + div_s[idx].re * div_s[idx].re;

  Result.u          := u;
  Result.v          := v;
  Result.omega      := omega;
  Result.psi        := psi_r;
  Result.ke         := ke_sum  * n_inv;
  Result.enstrophy  := ens_sum * n_inv;
  Result.max_vort   := max_w;
  Result.div_rms    := Sqrt(div_sum * n_inv);
  { cfl_peak already accumulated in time loop above }
end;

end.
