# run-pascal-phase7f-tests.ps1 -- native runner for the kh-sim Pascal phase-7f
# canonical-physics testbeds (CP-7F-A..E).
#
# Where each testbed runs (ruling 2, 2026-10-07: one shared harness per testbed):
#   CP-7F-A  tests/integration/test_canonical_physics_cp7fa.py  (MASTER s3.1 E9; moved 2026-10-12)
#   CP-7F-B  tests/integration/test_canonical_physics_cp7fbd.py (s3.2 E7; moved 2026-10-10)
#   CP-7F-C  tests/integration/test_canonical_physics_cp7fc.py  (s3.3 E6; moved 2026-10-07)
#   CP-7F-D  tests/integration/test_canonical_physics_cp7fbd.py (s3.4 E7; moved 2026-10-10)
#   CP-7F-E  here: the STRICT-C1 Pascal-binding anchor (tests/kh_cp7fe.pas, MASTER s3.5 E8):
#            E2 the d17r1 Dirichlet-Poisson grid, E1 the d12r6 convlv channels through
#            FFT1D, from shared/golden-vectors (found by walking up; KH_GOLDEN_DIR overrides).
# The test needs no running server (the in-process A and its server start are retired).
#
# Usage:
#   .\run-pascal-phase7f-tests.ps1            # build + run (A..D print where they run; E asserts)
#   .\run-pascal-phase7f-tests.ps1 -NoBuild   # reuse build\

param(
  [string]$Mode    = 'ALL',
  [string]$CaseArg = '',
  [switch]$NoBuild
)

$ErrorActionPreference = 'Stop'
Set-Location $PSScriptRoot

if (-not (Get-Command fpc -ErrorAction SilentlyContinue)) {
  Write-Error "fpc not found in PATH. Install: winget install FreePascal.FreePascal (then reopen shell)."
  exit 1
}
Write-Host ("FPC version: " + ((fpc -iV 2>&1) -join ''))

if (-not $NoBuild) {
  Write-Host "`n[build] canonical test ..."
  & "$PSScriptRoot\build.ps1" canonical
}

$testExe = Join-Path $PSScriptRoot 'build\test_canonical_physics.exe'
if (-not (Test-Path $testExe)) { Write-Error "test exe not found ($testExe); build first."; exit 1 }

Write-Host "[run] test_canonical_physics.exe`n"
& $testExe
$exitCode = $LASTEXITCODE
Write-Host "`n[done] test exit code: $exitCode"
exit $exitCode
