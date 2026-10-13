program canonical_physics;
{ kh-sim phase 7f -- Pascal canonical-physics testbeds (CP-7F-A..E).

  Doctrine: Pascal is the canonical PEDAGOGICAL reference language
  (LEGE-ARTIS-LANGUAGE-DOCTRINE v0.1 sec.3.3) -- strict C1-C4.

  Where each testbed runs (ruling 2, 2026-10-07: one shared harness per testbed):
    CP-7F-A  tests/integration/test_canonical_physics_cp7fa.py  (MASTER s3.1 E9; moved 2026-10-12)
    CP-7F-B  tests/integration/test_canonical_physics_cp7fbd.py (s3.2 E7; moved 2026-10-10)
    CP-7F-C  tests/integration/test_canonical_physics_cp7fc.py  (s3.3 E6; moved 2026-10-07)
    CP-7F-D  tests/integration/test_canonical_physics_cp7fbd.py (s3.4 E7; moved 2026-10-10)
    CP-7F-E  HERE, the Pascal-binding strict-C1 anchor (s3.5 E8): tests/kh_cp7fe.pas --
             Numerical Recipes in Pascal s17.5 (d17r1 Dirichlet-Poisson, printed grid
             within 2e-2) and s12.4 (d12r6 convlv through the production FFT1D, 1e-6).
  The in-process CP-7F-A that stood here (12 seeded modes, a single 256^2 run each,
  sigma from the [10, 30] window, gate [0.180, 0.200]) used the double layer at
  delta = Ly/20, where the two layers couple (assessment 2026-09-15, s7); its text is
  in the git history of this file.

  Usage: canonical_physics [A|B|C|D|E|all] [golden-vector dir for E]
  ASCII-only source (KB-039 -- FPC is strict about high-bit chars). }
{$mode objfpc}{$H+}

uses
  SysUtils, kh_cp7fe;

var
  g_fail: Boolean = False;   { a failed CP-7F-E check -> exit 1 }

procedure RunCp7fA;
begin
  WriteLn('CP-7F-A: SKIP -- moved to tests/integration/test_canonical_physics_cp7fa.py (MASTER s3.1 E9)');
end;

procedure RunCp7fB;
begin
  WriteLn('CP-7F-B: SKIP -- moved to tests/integration/test_canonical_physics_cp7fbd.py (MASTER s3.2 E7)');
end;

procedure RunCp7fC;
begin
  WriteLn('CP-7F-C: SKIP -- moved to tests/integration/test_canonical_physics_cp7fc.py (MASTER s3.3 E6)');
end;

procedure RunCp7fD;
begin
  WriteLn('CP-7F-D: SKIP -- moved to tests/integration/test_canonical_physics_cp7fbd.py (MASTER s3.4 E7)');
end;

procedure RunCp7fEMode(const goldenDir: string);
begin
  if not RunCp7fE(goldenDir) then g_fail := True;
end;

var
  mode, goldenDir: string;
begin
  if ParamCount >= 1 then mode := ParamStr(1) else mode := 'all';
  goldenDir := '';
  if (mode = 'E') and (ParamCount >= 2) then goldenDir := ParamStr(2);

  WriteLn('=================================================');
  WriteLn(' kh-sim phase 7f: Pascal canonical-physics testbeds');
  WriteLn(' (A, B, C, D in the shared harnesses; E asserts here)');
  WriteLn('=================================================');

  if mode = 'all' then
  begin
    RunCp7fA; RunCp7fB; RunCp7fC; RunCp7fD; WriteLn;
    RunCp7fEMode(goldenDir); WriteLn;
  end
  else if mode = 'A' then RunCp7fA
  else if mode = 'B' then RunCp7fB
  else if mode = 'C' then RunCp7fC
  else if mode = 'D' then RunCp7fD
  else if mode = 'E' then RunCp7fEMode(goldenDir)
  else begin WriteLn('Unknown mode: ', mode); Halt(2); end;

  WriteLn('-------------------------------------------------');
  if g_fail then begin WriteLn('RESULT: FAIL'); Halt(1); end
  else begin WriteLn('RESULT: PASS'); Halt(0); end;
end.
