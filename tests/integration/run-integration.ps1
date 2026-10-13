# run-integration.ps1 -- KH-018 integration test runner (ThinkPad LDE)
# Task: KH-018
#
# PREREQUISITES
#   LDE stack running:  .\_config\Start-LocalEnv.ps1 -Action up
#   Python + pip:       python --version >= 3.10
#
# USAGE
#   .\kh-sim\tests\integration\run-integration.ps1           # all tests
#   .\kh-sim\tests\integration\run-integration.ps1 -Subset health   # health only
#   .\kh-sim\tests\integration\run-integration.ps1 -Subset simulate # physics only
#   .\kh-sim\tests\integration\run-integration.ps1 -Subset log      # log service
#   .\kh-sim\tests\integration\run-integration.ps1 -Verbose         # -v output
#
# EXIT CODES
#   0 = all tests passed
#   1 = test failures or setup error
#
# NOTE: ASCII-only characters throughout -- PS 5.1 parse errors on em dashes.

[CmdletBinding()]
param(
    [ValidateSet("all","health","info","simulate","log","cross")]
    [string]$Subset = "all",

    [switch]$VerboseOutput,   # Pass -v to pytest for per-test output
    [switch]$SkipInstall       # Skip pip install (use if deps already installed)
)

Set-StrictMode -Version Latest
$ErrorActionPreference = "Stop"

$ScriptDir   = Split-Path -Parent $MyInvocation.MyCommand.Path
$ProjectRoot = Split-Path -Parent (Split-Path -Parent $ScriptDir)
$ReqFile     = Join-Path $ScriptDir "requirements.txt"

function Write-Header([string]$Text) {
    $line = "-" * 60
    Write-Host ""
    Write-Host $line -ForegroundColor DarkCyan
    Write-Host "  $Text" -ForegroundColor Cyan
    Write-Host $line -ForegroundColor DarkCyan
}
function Write-OK([string]$Text)   { Write-Host "  [OK]   $Text" -ForegroundColor Green  }
function Write-FAIL([string]$Text) { Write-Host "  [FAIL] $Text" -ForegroundColor Red    }
function Write-INFO([string]$Text) { Write-Host "  [INFO] $Text" -ForegroundColor Yellow }

Write-Header "KH-018 Integration Tests -- ThinkPad LDE"

# ── Quick LDE health gate ─────────────────────────────────────────────────────
Write-Header "Pre-flight: LDE health check"
$endpoints = @(
    @{ name = "kh-rust";       url = "http://localhost:8001/health" },
    @{ name = "kh-scala";      url = "http://localhost:8002/health" },
    @{ name = "kh-cpp";        url = "http://localhost:8003/health" },
    @{ name = "kh-fortran";    url = "http://localhost:8004/health" },
    @{ name = "kh-pascal";     url = "http://localhost:8005/health" },
    @{ name = "kh-log-service";url = "http://localhost:8006/health" }
)
$reachable = 0
foreach ($ep in $endpoints) {
    try {
        $r = Invoke-WebRequest -Uri $ep.url -UseBasicParsing -TimeoutSec 3 -ErrorAction Stop
        if ($r.StatusCode -lt 400 -or $r.StatusCode -eq 503) {
            Write-OK  "$($ep.name) --> $($ep.url)"
            $reachable++
        } else {
            Write-FAIL "$($ep.name) --> HTTP $($r.StatusCode)"
        }
    } catch {
        Write-FAIL "$($ep.name) --> unreachable (start LDE stack first)"
    }
}
Write-INFO "$reachable / $($endpoints.Count) services reachable"
if ($reachable -lt 2) {
    Write-FAIL "Fewer than 2 services reachable -- aborting. Start LDE stack:"
    Write-INFO '  .\_config\Start-LocalEnv.ps1 -Action up'
    exit 1
}

# ── Install test dependencies ─────────────────────────────────────────────────
if (-not $SkipInstall) {
    Write-Header "Installing test dependencies"
    pip install -r $ReqFile --quiet --break-system-packages
    if ($LASTEXITCODE -ne 0) {
        Write-FAIL "pip install failed"
        exit 1
    }
    Write-OK "Dependencies installed."
}

# ── Build pytest arguments ────────────────────────────────────────────────────
$pytestArgs = @($ScriptDir)

switch ($Subset) {
    "health"   { $pytestArgs += @("-k", "health") }
    "info"     { $pytestArgs += @("-k", "info") }
    "simulate" { $pytestArgs += @("-k", "simulate") }
    "log"      { $pytestArgs += @("-k", "log") }
    "cross"    { $pytestArgs += @("-k", "cross") }
    default    { }   # all
}

if ($VerboseOutput) {
    $pytestArgs += "-v"
} else {
    $pytestArgs += "-v"   # always -v for integration suites -- per-test pass/fail is informative
}

$pytestArgs += "--tb=short"

# ── Run ───────────────────────────────────────────────────────────────────────
Write-Header "pytest $($pytestArgs -join ' ')"
Set-Location $ProjectRoot

python -m pytest @pytestArgs
$exitCode = $LASTEXITCODE

Write-Host ""
if ($exitCode -eq 0) {
    Write-OK  "KH-018 integration suite: ALL PASSED"
} else {
    Write-FAIL "KH-018 integration suite: FAILURES (exit $exitCode)"
    Write-INFO "Run with -Subset <module> to isolate failing group"
    Write-INFO "Detailed output: python -m pytest $ScriptDir -v --tb=long"
}
exit $exitCode
