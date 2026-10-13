program test_canonical_physics;
{ kh-sim phase 7f -- Pascal canonical-physics testbeds (CP-7F-A..E).

  Doctrine status: Pascal is the canonical pedagogical reference language
  (LEGE-ARTIS-LANGUAGE-DOCTRINE v0.1 §3.3). Subject to §4.4 C1-C4 in full.
  CP-7F-E is the Pascal-binding STRICT-C1 anchor: a verbatim reproduction of a
  worked example from the direct historical lineage text --
  Press, Flannery, Teukolsky, Vetterling (1986), Numerical Recipes in Pascal,
  Chapter 19 (Partial Differential Equations). Pascal CP-7F-E green is the
  single strongest empirical event for promoting the doctrine v0.1 -> v1.0.

  This file is the SIBLING of phase-7d tests/test_physics_testbeds.pas.
  C3 discharge: every test procedure carries a citation header block.

  ----------------------------------------------------------------------------
  STATUS 2026-10-11: CP-7F-E ASSERTS (tests/kh_cp7fe.pas): E2 the d17r1
    SOR Dirichlet-Poisson grid (printed, 2e-2) by the NR s17.4 sine-transform
    method, E1 the d12r6 convlv channels (1e-6) through the production FFT1D.
    The loader stub that stood here read a draft schema (expected.value) the
    transcription of 2026-06-17 does not use; it is in git history.
    CP-7F-A, B, C and D run from the shared harnesses (A since 2026-10-12, MASTER
    E9: test_canonical_physics_cp7fa.py, green for the oracle and all five backends on
    the station that day).
  ----------------------------------------------------------------------------
}
{$mode objfpc}{$H+}

uses
  SysUtils, kh_cp7fe;

var
  g_pass: Boolean = True;
  g_skip: Boolean = False;

{ ===========================================================================
  CP-7F-E: Numerical Recipes in Pascal -- the Pascal-binding strict-C1 anchor.
  E2 (s17.5, d17r1 SOR Dirichlet-Poisson) and E1 (s12.4, d12r6 convlv through
  the production FFT1D), implemented once in tests/kh_cp7fe.pas and shared with
  canonical_physics.pas. ParamStr(1) may name the golden-vector directory;
  otherwise KH_GOLDEN_DIR or a walk up to shared/golden-vectors.
  =========================================================================== }
procedure test_cp7fe_nr1986(const goldenDir: string);
begin
  if not RunCp7fE(goldenDir) then g_pass := False;
end;

{ -- CP-7F-A, B, C, D: run from the shared harnesses (test_canonical_physics_cp7fa.py,
  ..._cp7fbd.py, ..._cp7fc.py). Pointers only.}
procedure test_cp7fa_kh_dispersion;      begin WriteLn('CP-7F-A: SKIP -- moved to tests/integration/test_canonical_physics_cp7fa.py (MASTER s3.1 E9)'); end;
procedure test_cp7fb_enstrophy_cascade;  begin WriteLn('CP-7F-B: SKIP -- moved to tests/integration/test_canonical_physics_cp7fbd.py (MASTER s3.2 E7)'); end;
procedure test_cp7fc_vortex_pair_merger; begin WriteLn('CP-7F-C: SKIP -- moved to tests/integration/test_canonical_physics_cp7fc.py (MASTER s3.3 E6)'); end;
procedure test_cp7fd_forced_taylor_green;begin WriteLn('CP-7F-D: SKIP -- moved to tests/integration/test_canonical_physics_cp7fbd.py (MASTER s3.4 E7)'); end;

{ ===========================================================================
  Main
  =========================================================================== }
var
  goldenDir: string;
begin
  { run-pascal-phase7f-tests.ps1 passes its -Mode here; only a directory is a
    golden-vector directory }
  goldenDir := '';
  if (ParamCount >= 1) and DirectoryExists(ParamStr(1)) then goldenDir := ParamStr(1);

  WriteLn('=================================================');
  WriteLn(' kh-sim phase 7f: Pascal canonical-physics testbeds');
  WriteLn(' (A, B, C, D in the shared harnesses; E asserts here)');
  WriteLn('=================================================');
  WriteLn;

  test_cp7fa_kh_dispersion;       WriteLn;
  test_cp7fb_enstrophy_cascade;   WriteLn;
  test_cp7fc_vortex_pair_merger;  WriteLn;
  test_cp7fd_forced_taylor_green; WriteLn;
  test_cp7fe_nr1986(goldenDir); WriteLn;

  WriteLn('-------------------------------------------------');
  if not g_pass then begin
    WriteLn('RESULT: FAIL'); Halt(1);
  end else if g_skip then begin
    WriteLn('RESULT: PASS (with SKIP/TODO -- not yet phase-7f-pascal-green)');
    Halt(0);
  end else begin
    WriteLn('RESULT: PASS'); Halt(0);
  end;
end.
