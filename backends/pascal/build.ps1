# build.ps1 — KH-SIM Pascal backend (KH-007)
# Requirements: Free Pascal Compiler (fpc) in PATH
#
# Install FPC on Windows (if not present):
#   winget install FreePascal.FreePascal
# Or via MSYS2 (preferred for consistency with other backends):
#   pacman -S mingw-w64-x86_64-fpc
#
# Usage:
#   .\build.ps1             # build everything (validate + server)
#   .\build.ps1 validate    # build validation binary only
#   .\build.ps1 server      # build HTTP server only
#   .\build.ps1 canonical   # build phase-7f canonical-physics test only
#   .\build.ps1 clean       # remove build artifacts

param(
  [ValidateSet('all', 'validate', 'server', 'canonical', 'clean')]
  [string]$Target = 'all'
)

$ErrorActionPreference = 'Stop'
Set-Location $PSScriptRoot

# ── Check FPC ─────────────────────────────────────────────────────────────────
if (-not (Get-Command fpc -ErrorAction SilentlyContinue)) {
  Write-Error @"
fpc not found in PATH.
Install options:
  winget install FreePascal.FreePascal
  pacman -S mingw-w64-x86_64-fpc   (MSYS2 MinGW64 shell)
After install, open a new shell so PATH is refreshed.
"@
  exit 1
}

$fpcVer = (fpc -iV 2>&1) -join ''
Write-Host "FPC version: $fpcVer"

# ── FPC flags ─────────────────────────────────────────────────────────────────
# -O2          optimise
# -Fusrc       search src/ for units
# -FUbuild     write compiled units (.ppu/.o) to build/
# -FEbuild     write executables to build/
# -Fubuild     also search build/ for pre-compiled units (shared across targets)
$FLAGS = @('-O2', '-Fusrc', '-Futests', '-Fubuild', '-FUbuild', '-FEbuild')

if (-not (Test-Path 'build')) { New-Item -ItemType Directory 'build' | Out-Null }

function Invoke-FPC([string]$src, [string]$label) {
  Write-Host "`n  [$label] fpc $src"
  & fpc @FLAGS $src
  if ($LASTEXITCODE -ne 0) {
    throw "fpc compilation failed for $src (exit $LASTEXITCODE)"
  }
  $exe = [IO.Path]::GetFileNameWithoutExtension($src)
  Write-Host "  -> build\$exe.exe"
}

switch ($Target) {
  'validate' {
    Invoke-FPC 'tests\kh_validate.pas' 'validate'
  }
  'server' {
    Invoke-FPC 'src\kh_server.pas' 'server'
  }
  'all' {
    # validate first: physics-only, no fcl-web dependency — fast smoke test
    Invoke-FPC 'tests\kh_validate.pas' 'validate'
    Invoke-FPC 'src\kh_server.pas'     'server'
    Invoke-FPC 'tests\test_canonical_physics.pas' 'canonical'
  }
  'canonical' {
    # Phase 7f canonical-physics testbed (CP-7F-A..E). Requires the running
    # HTTP server (.\build.ps1 server then run build\kh_server.exe).
    Invoke-FPC 'tests\test_canonical_physics.pas' 'canonical'
  }
  'clean' {
    Get-ChildItem 'build' -Include '*.o','*.ppu','*.exe' -Recurse |
      Remove-Item -Force
    Write-Host 'build/ cleaned'
  }
}

Write-Host "`nDone."
