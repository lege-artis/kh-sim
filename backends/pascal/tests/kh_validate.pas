program kh_validate;
{ KH-007 Pascal validation - compares against the oracle golden vector V1
  (shared/golden-vectors/v1-canonical-reference.json, class A2), which replaced
  the retired shared/physics/kh_reference_output.json on 2026-10-08.
  The request is read from V1 "parameters"; the tolerance is V1's A2 band:
  |got - ref| <= max(scalar_abs * f, scalar_rel * f * |ref|), f = the
  "kh-pascal" scalar_factor override (1 when absent).
  Exit code 0 = all passed, 1 = any failure.
}
{$mode objfpc}{$H+}

uses
  SysUtils, Math, Classes, fpjson, jsonparser,
  kh_physics;

const
  BACKEND_ID = 'kh-pascal';
  DIV_TOL    = 1.0e-10; { absolute threshold for divergence_rms }

{ -- Golden vector loader --------------------------------------------------- }
type
  TRef = record
    ke, enstrophy, max_vort, div_rms: Double;
    nx, ny, steps, mode: Integer;
    lx, ly, dt, re, u0, amp: Double;
    s_abs, s_rel, factor: Double;
  end;

function LoadRef(const fname: string): TRef;
var
  sl:   TStringList;
  jo, p, diag, tols, ov: TJSONObject;
  jd:   TJSONData;
begin
  sl := TStringList.Create;
  try
    sl.LoadFromFile(fname);
    jo := GetJSON(sl.Text) as TJSONObject;
    try
      p    := jo.Objects['parameters'];
      diag := jo.Objects['diagnostics'];
      tols := jo.Objects['tolerances'];
      Result.nx    := p.Integers['grid_nx'];
      Result.ny    := p.Integers['grid_ny'];
      Result.steps := p.Integers['nsteps'];
      Result.mode  := p.Integers['perturbation_mode'];
      Result.lx    := p.Floats['domain_lx'];
      Result.ly    := p.Floats['domain_ly'];
      Result.dt    := p.Floats['dt'];
      Result.re    := p.Floats['reynolds_number'];
      Result.u0    := p.Floats['velocity_shear'];
      Result.amp   := p.Floats['perturbation_amplitude'];

      Result.ke        := diag.Floats['kinetic_energy'];
      Result.enstrophy := diag.Floats['enstrophy'];
      Result.max_vort  := diag.Floats['max_vorticity'];
      Result.div_rms   := diag.Floats['divergence_rms'];

      Result.factor := 1.0;
      jd := tols.FindPath('backend_overrides.' + BACKEND_ID);
      if Assigned(jd) and (jd.JSONType = jtObject) then begin
        ov := TJSONObject(jd);
        jd := ov.Find('scalar_factor');
        if Assigned(jd) then Result.factor := jd.AsFloat;
      end;
      Result.s_abs := tols.Floats['scalar_abs'] * Result.factor;
      Result.s_rel := tols.Floats['scalar_rel'] * Result.factor;
    finally
      jo.Free;
    end;
  finally
    sl.Free;
  end;
end;

{ -- Check helper ----------------------------------------------------------- }
var
  g_pass: Boolean = True;

procedure Check(const name: string; got, refval, s_abs, s_rel: Double);
var
  diff, band: Double;
begin
  diff := Abs(got - refval);
  band := Max(s_abs, s_rel * Abs(refval));
  Write(Format('%-15s : got=%.17g  ref=%.17g  |diff|=%.3e  band=%.3e',
               [name, got, refval, diff, band]));
  if not (diff <= band) then begin
    WriteLn('  FAIL');
    g_pass := False;
  end else
    WriteLn;
end;

procedure CheckDiv(got: Double);
begin
  Write(Format('%-15s : got=%.3e', ['divergence_rms', got]));
  if got > DIV_TOL then begin
    WriteLn(Format('  FAIL (>%.0e)', [DIV_TOL]));
    g_pass := False;
  end else
    WriteLn('  OK');
end;

{ -- Main ------------------------------------------------------------------- }
var
  ref: TRef;
  req: TSimRequest;
  res: TSimResult;
begin
  if ParamCount < 1 then begin
    WriteLn('Usage: kh-sim-pascal-validate <v1-canonical-reference.json>');
    Halt(1);
  end;

  ref := LoadRef(ParamStr(1));

  with req do begin
    nx := ref.nx;  ny := ref.ny;
    lx := ref.lx;  ly := ref.ly;
    dt := ref.dt;  steps := ref.steps;
    re := ref.re;  u0 := ref.u0;
    amp := ref.amp;  mode := ref.mode;
    { Explicit no-forcing defaults (phase 7f). A local TSimRequest record is
      NOT zero-initialised for unmanaged fields; without these, garbage in
      has_init_omega / linear_drag would corrupt the pre-7f RHS.             }
    has_init_omega := False;
    forcing_type      := 'none';
    forcing_amplitude := 0.0;
    forcing_kx := 1;  forcing_ky := 1;
    forcing_kmin := 0.0;  forcing_kmax := 0.0;
    forcing_omega0 := 0.0;  linear_drag := 0.0;  t_start := 0.0;
  end;

  WriteLn(Format('--- Pascal vs golden vector V1 (A2 band, factor %g) ---', [ref.factor]));
  res := Simulate(req);

  Check('kinetic_energy', res.ke,        ref.ke,        ref.s_abs, ref.s_rel);
  Check('enstrophy',      res.enstrophy, ref.enstrophy, ref.s_abs, ref.s_rel);
  Check('max_vorticity',  res.max_vort,  ref.max_vort,  ref.s_abs, ref.s_rel);
  CheckDiv(res.div_rms);
  WriteLn(Format('%-15s : %.3e ms', ['compute_ms', res.compute_ms]));
  WriteLn;

  if g_pass then
    WriteLn('All tests PASSED')
  else begin
    WriteLn('TESTS FAILED');
    Halt(1);
  end;
end.
