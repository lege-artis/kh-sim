program test_solver_guards;
{ kh-sim Pascal backend - native locks for the 2026-09-26 solver changes.
  Finding: kh-sim/_audit/FINDING-PASCAL-HTTP-PARAMETER-DROP-2026-09-16.md (v0.4)

  Calls kh_physics.Simulate directly, bypassing the HTTP validation in
  kh_server.pas, so it exercises the guards the HTTP tests cannot reach:

  G1/G2  FFT1D refuses a non-power-of-two length instead of writing outside
         its array. Over HTTP this is unreachable by design - the server
         rejects such grids first - so without this program the guard would be
         a lock that nothing can turn. Build it WITHOUT -Cr (as the Dockerfile
         does) to test the condition in which the guard matters: with range
         checks off, a missing guard corrupts the heap silently.
  G3     steps_completed and t_final are reported by the solver.
  G4     a negative step count reports what was done (zero steps), not what
         was asked.

  Build:  fpc -O2 -Fusrc -Fubuild -FUbuild -FEbuild tests/test_solver_guards.pas
  Run:    ./build/test_solver_guards      (exit 0 = all pass)
}
{$mode objfpc}{$H+}

uses SysUtils, kh_physics;

var
  failures: Integer = 0;

function BaseRequest(nx, ny, steps: Integer): TSimRequest;
begin
  Result.nx := nx;  Result.ny := ny;  Result.steps := steps;  Result.mode := 2;
  Result.lx := 1.0; Result.ly := 0.5; Result.dt := 0.001;
  Result.re := 1000.0; Result.u0 := 1.0; Result.amp := 0.01;
  Result.has_init_omega := False;
  SetLength(Result.init_omega, 0);
  Result.forcing_type := 'none';  Result.forcing_amplitude := 0.0;
  Result.forcing_kx := 1;  Result.forcing_ky := 1;
  Result.forcing_kmin := 0.0;  Result.forcing_kmax := 0.0;
  Result.forcing_omega0 := 0.0;  Result.linear_drag := 0.0;  Result.t_start := 0.0;
end;

procedure Check(const name: string; ok: Boolean; const detail: string);
begin
  if ok then WriteLn('PASS  ', name)
  else begin
    WriteLn('FAIL  ', name, '  -- ', detail);
    Inc(failures);
  end;
end;

procedure ExpectRefusal(const name: string; nx, ny: Integer);
var res: TSimResult; msg: string; raised: Boolean;
begin
  raised := False;  msg := '';
  try
    res := Simulate(BaseRequest(nx, ny, 2));
  except
    on E: Exception do begin raised := True; msg := E.Message; end;
  end;
  Check(name, raised and (Pos('not a power of two', msg) > 0),
        Format('grid %dx%d: raised=%s message="%s"', [nx, ny, BoolToStr(raised, True), msg]));
end;

var
  res: TSimResult;
begin
  WriteLn('kh-sim Pascal: solver guard locks');

  ExpectRefusal('G1 FFT1D refuses non-power-of-two nx (48x32)', 48, 32);
  ExpectRefusal('G2 FFT1D refuses non-power-of-two ny (64x24)', 64, 24);

  res := Simulate(BaseRequest(64, 32, 5));
  Check('G3a steps_completed reported by the solver', res.steps_completed = 5,
        Format('got %d, expected 5', [res.steps_completed]));
  Check('G3b t_final reported by the solver', Abs(res.t_final - 5 * 0.001) < 1e-15,
        Format('got %.17g, expected 0.005', [res.t_final]));

  res := Simulate(BaseRequest(64, 32, -3));
  Check('G4 negative steps reports zero steps done', (res.steps_completed = 0) and (res.t_final = 0.0),
        Format('steps_completed=%d t_final=%.17g', [res.steps_completed, res.t_final]));

  if failures = 0 then begin WriteLn('RESULT: PASS'); Halt(0); end
  else begin WriteLn(Format('RESULT: FAIL (%d)', [failures])); Halt(1); end;
end.
