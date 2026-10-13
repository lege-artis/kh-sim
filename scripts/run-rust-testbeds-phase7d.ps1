# run-rust-testbeds-phase7d.ps1
# Phase-7d physics testbed runner for the Rust backend.
#
# Usage (from kh-sim repo root):
#   cd /path/to/kh-sim
#   .\scripts\run-rust-testbeds-phase7d.ps1
#
# Prerequisites:
#   - Python 3 with: pip install pytest requests numpy
#   - Pre-built Rust server binary OR Rust toolchain (cargo) installed
#
# The script:
#   1. Starts kh-sim-rust.exe on port 8001 (pre-built or cargo-run)
#   2. Waits for /health to respond
#   3. Runs pytest on test_physics_testbeds_rust.py with -v --tb=short
#   4. Kills the server
#   5. Exits with the pytest exit code

$ErrorActionPreference = "Stop"
$KhSimRoot = Split-Path -Parent $PSScriptRoot

Write-Host "=== kh-sim phase-7d Rust testbeds ===" -ForegroundColor Cyan
Write-Host "Monorepo root: $KhSimRoot"

# ── Step 1: start the server ──────────────────────────────────────────────────

$ServerExe  = Join-Path $KhSimRoot "backends\rust\target\release\kh-sim-rust.exe"
$ServerPort = 8001
$ServerUrl  = "http://localhost:$ServerPort/health"

if (Test-Path $ServerExe) {
    Write-Host "Starting pre-built server: $ServerExe" -ForegroundColor Green
    $ServerProc = Start-Process -FilePath $ServerExe -PassThru -WindowStyle Hidden
} else {
    Write-Host "Pre-built binary not found; using 'cargo run --release'" -ForegroundColor Yellow
    $RustDir = Join-Path $KhSimRoot "backends\rust"
    $ServerProc = Start-Process -FilePath "cargo" -ArgumentList "run","--release" `
                                -WorkingDirectory $RustDir -PassThru -WindowStyle Hidden
}

Write-Host "Server PID: $($ServerProc.Id)"

# ── Step 2: wait for server health ───────────────────────────────────────────

Write-Host "Waiting for server on port $ServerPort ..."
$maxWait = 30   # seconds
$waited  = 0
$ready   = $false
while ($waited -lt $maxWait) {
    Start-Sleep -Seconds 1
    $waited++
    try {
        $r = Invoke-WebRequest -Uri $ServerUrl -UseBasicParsing -TimeoutSec 1 -ErrorAction Stop
        if ($r.StatusCode -eq 200) { $ready = $true; break }
    } catch { }
}

if (-not $ready) {
    Write-Host "ERROR: server did not start within ${maxWait}s" -ForegroundColor Red
    if (-not $ServerProc.HasExited) { Stop-Process -Id $ServerProc.Id -Force }
    exit 1
}
Write-Host "Server ready." -ForegroundColor Green

# ── Step 3: run pytest ────────────────────────────────────────────────────────

$TestFile = Join-Path $KhSimRoot "tests\integration\test_physics_testbeds_rust.py"
Write-Host "Running: pytest $TestFile -v --tb=short" -ForegroundColor Cyan

$env:KH_RUST_REQUIRE_SERVER = "1"
& python -m pytest $TestFile -v --tb=short
$PytestExit = $LASTEXITCODE

# ── Step 4: stop the server ───────────────────────────────────────────────────

Write-Host "Stopping server (PID $($ServerProc.Id))..."
if (-not $ServerProc.HasExited) {
    Stop-Process -Id $ServerProc.Id -Force
}

# ── Step 5: report ────────────────────────────────────────────────────────────

if ($PytestExit -eq 0) {
    Write-Host "=== ALL TESTS PASSED ===" -ForegroundColor Green
} else {
    Write-Host "=== TESTS FAILED (exit $PytestExit) ===" -ForegroundColor Red
}
exit $PytestExit
