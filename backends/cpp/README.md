# KH-SIM — C++ Backend

**Task:** KH-005 | **Framework:** cpp-httplib 0.18 | **Port:** 8003 | **FFT:** Cooley-Tukey (built-in)

## Architecture

```
CMakeLists.txt              CMake 3.20+, FetchContent (cpp-httplib, nlohmann/json)
src/
  main.cpp                  cpp-httplib server on :8003, CORS, JSON error handling
  models.hpp                SimulationRequest / SimulationResult (nlohmann/json codecs)
  physics.hpp               simulate() declaration
  physics.cpp               Cooley-Tukey FFT, IC, Poisson, velocity recovery, RK4
tests/
  validation_test.cpp       standalone binary vs golden vector V1 (A2 band)
```

## Endpoint contract

Implements `kh-sim/shared/api/openapi.yaml`

| Endpoint  | Method | Description                                |
|-----------|--------|--------------------------------------------|
| /simulate | POST   | Run KH simulation, return SimulationResult |
| /health   | GET    | Liveness probe                             |
| /info     | GET    | Backend language/framework metadata        |

## Physics kernel

Exact port of `kh-sim/shared/physics/kh_physics.py`:
- IC: tanh shear + sinusoidal perturbation
- FFT: in-place Cooley-Tukey radix-2 (std::complex<double>, C++20)
- Poisson: psi_hat = omega_hat / k^2, zero mean
- Velocity: u = i*ky*psi_hat, v = -i*kx*psi_hat
- Time: RK4, 4 stages per step

No external FFT library. Requires grid dimensions to be powers of 2.

## Build and run

Prerequisites: CMake 3.20+, C++20 compiler (MSVC 19.29+ or GCC 11+), Git.
All HTTP and JSON dependencies are downloaded automatically at configure time.

```powershell
cd kh-sim\backends\cpp

# Configure (downloads cpp-httplib + nlohmann/json via FetchContent)
cmake -B build -DCMAKE_BUILD_TYPE=Release

# Build server + validation binary
cmake --build build --config Release

# Run server (binds 0.0.0.0:8003)
.\build\Release\kh-sim-cpp.exe
```

## Validation test

```powershell
# Prerequisite
python kh-sim\shared\physics\kh_physics.py

# Build and run validation binary directly
cmake --build build --config Release
.\build\Release\kh-sim-cpp-validate.exe ..\..\shared\golden-vectors\v1-canonical-reference.json

# Or via CTest
ctest --test-dir build -C Release -V
```

Expected output:
```
--- C++ vs Python reference ---
kinetic_energy : got=0.112810  ref=0.112810  err=0.00%
enstrophy      : got=43.705142  ref=43.705142  err=0.00%
max_vorticity  : got=31.907572  ref=31.907572  err=0.00%
divergence_rms : got=~1e-14    ref=1.20e-14
All tests PASSED
```

## Smoke test (server running)

```powershell
curl http://localhost:8003/health

$body = '{"grid_nx":64,"grid_ny":32,"steps":100,"reynolds_number":1000}'
Invoke-RestMethod -Method POST -Uri http://localhost:8003/simulate `
    -ContentType "application/json" -Body $body | ConvertTo-Json -Depth 5

```

## Status

v0.1.0 (2026-10), validated against the Python reference oracle on every layer of the suite:

- golden vectors v1–v9 of `shared/golden-vectors/`, element by element (tolerance factor 1×);
- the physics testbeds (conservation, decay, linear growth) of `tests/integration/`;
- the canonical textbook testbeds CP-7F-A…D, run once for all six implementations from the shared harnesses in `tests/integration/test_canonical_physics_cp7f*.py`;
- request validation and response schema per `shared/api/openapi.yaml` (`tests/integration/test_conformance.py`).

Container: `docker build -t kh-sim-cpp .` then `docker run --rm -p 8003:8003 kh-sim-cpp`.
