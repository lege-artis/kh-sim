# run-fortran-phase7f-tests.ps1
# Phase-7F canonical physics testbeds runner for ThinkPad (gfortran native).
# Builds test_canonical_physics from source, then runs all test cases.
#
# Usage:
#   .\run-fortran-phase7f-tests.ps1                   # full run (all cases)
#   .\run-fortran-phase7f-tests.ps1 -Mode A            # print SKIP: moved to the shared harness (E9)
#   .\run-fortran-phase7f-tests.ps1 -Mode B|C|D        # print SKIP: moved to the shared harnesses
#   .\run-fortran-phase7f-tests.ps1 -Mode E            # CP-7F-E: NR in Pascal E2 + E1 (MASTER E8)
#
# Acceptance gates:
#   CP-7F-B, CP-7F-D: moved 2026-10-10 to tests/integration/test_canonical_physics_cp7fbd.py
#             (MASTER s3.2/s3.4 E7; the HTTP server computes forcing since then)
#   CP-7F-A: moved 2026-10-12 to tests/integration/test_canonical_physics_cp7fa.py (MASTER s3.1 E9)
#   CP-7F-C: moved 2026-10-07 to tests/integration/test_canonical_physics_cp7fc.py (MASTER s3.3 E6)

param(
    [string]$Mode = "all",
    [string]$GFortran = "gfortran"
)

$ErrorActionPreference = "Stop"
$RepoRoot = Split-Path -Parent (Split-Path -Parent $PSScriptRoot)
$SrcDir   = Join-Path $RepoRoot "src"
$TestsDir = Join-Path $RepoRoot "tests"
$ModDir   = Join-Path $env:TEMP "kh_mods_phase7f"
$BinPath  = Join-Path $env:TEMP "test_canonical_physics.exe"
$CheckpointDir = $env:TEMP
$script:Failed = $false

# -- Build --
Write-Host "=== Building test_canonical_physics ==="
if (-not (Test-Path $ModDir)) { New-Item -ItemType Directory -Force -Path $ModDir | Out-Null }

$Sources = @(
    "kh_constants.f90", "kh_grid.f90", "kh_fft.f90",
    "kh_initial_conditions.f90", "kh_nonlinear.f90",
    "kh_poisson.f90", "kh_velocity.f90", "kh_diagnostics.f90",
    "kh_forcing.f90", "kh_physics.f90", "kh_io.f90"
) | ForEach-Object { Join-Path $SrcDir $_ }
$Sources += Join-Path $TestsDir "test_canonical_physics.f90"

$CompileArgs = @("-std=f2008", "-O2", "-J$ModDir", "-I$ModDir") + $Sources + @("-o", $BinPath)
& $GFortran @CompileArgs
if ($LASTEXITCODE -ne 0) { Write-Error "Compilation failed"; exit 1 }
Write-Host "  Binary: $BinPath"

# -- Helper: run a test-binary invocation --
function Run-Test {
    param([string[]]$Args)
    & $BinPath @Args
    if ($LASTEXITCODE -ne 0 -and $LASTEXITCODE -ne 1) {
        Write-Warning "test_canonical_physics returned exit=$LASTEXITCODE"
    }
}

# -- CP-7F-B, CP-7F-D: moved to tests/integration/test_canonical_physics_cp7fbd.py (E7) --
function Run-CPBD {
    Write-Host ""
    Run-Test @("B")
    Run-Test @("D")
}

# -- CP-7F-A: moved 2026-10-12 to tests/integration/test_canonical_physics_cp7fa.py (E9) --
function Run-CPA {
    Write-Host ""
    Run-Test @("A")
}

# -- CP-7F-C: moved to tests/integration/test_canonical_physics_cp7fc.py (E6) --
function Run-CPC {
    Write-Host ""
    Run-Test @("C")
}

# -- CP-7F-E: NR in Pascal s17.5 (E2) + s12.4 convlv (E1), MASTER E8 --
function Run-CPE {
    Write-Host ""
    $env:KH_FORTRAN_TESTS_DIR = $TestsDir
    try {
        & $BinPath "E"
        if ($LASTEXITCODE -ne 0) { $script:Failed = $true }
    } finally { Remove-Item Env:KH_FORTRAN_TESTS_DIR -ErrorAction SilentlyContinue }
}

# -- Dispatch --
switch ($Mode.ToLower()) {
    "all" {
        Run-CPA
        Run-CPBD
        Run-CPC
        Run-CPE
    }
    "e"       { Run-CPE }
    "b"       { Run-CPBD }
    "d"       { Run-CPBD }
    "a"       { Run-CPA }
    "c"       { Run-CPC }
    default   { Write-Error "Unknown mode '$Mode'. Use: all, A, B, C, D, E" }
}
if ($script:Failed) { exit 1 }
