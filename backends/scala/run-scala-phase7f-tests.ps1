# run-scala-phase7f-tests.ps1
# Phase-7F canonical physics testbeds, Scala backend.
#
# Every CP-7F testbed runs once, for the oracle and all five backends, from a shared
# harness in tests/integration (ruling 2, 2026-10-07):
#   CP-7F-A: test_canonical_physics_cp7fa.py   (MASTER s3.1 E9; moved 2026-10-12)
#   CP-7F-B: test_canonical_physics_cp7fbd.py  (MASTER s3.2 E7; moved 2026-10-10)
#   CP-7F-C: test_canonical_physics_cp7fc.py   (MASTER s3.3 E6; moved 2026-10-07)
#   CP-7F-D: test_canonical_physics_cp7fbd.py  (MASTER s3.4 E7; moved 2026-10-10)
#   CP-7F-E: the Pascal-binding anchor (MASTER s3.5 E8): the oracle, Pascal and Fortran
# khsim.CanonicalPhysicsTest only prints these pointers; it needs no running server.
#
# Usage:  .\run-scala-phase7f-tests.ps1      (the former -Mode / -KeepServer are ignored)

param(
    [string]$Mode = "native",
    [switch]$KeepServer
)

$ErrorActionPreference = "Stop"
Set-Location $PSScriptRoot
Write-Host "=== sbt testOnly khsim.CanonicalPhysicsTest (pointers to the shared harnesses) ==="
sbt "testOnly khsim.CanonicalPhysicsTest"
exit $LASTEXITCODE
