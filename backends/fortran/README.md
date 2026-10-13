# KH-SIM — Fortran Backend

**Task:** KH-006 | **Port:** 8004 | **FFT:** Cooley-Tukey (native Fortran 2008)

## Architecture

```
CMakeLists.txt              Mixed-language: Fortran + C++ (gfortran + g++, Ninja)
src/
  kh_physics.f90            Fortran module: FFT, IC, Poisson, RK4, iso_c_binding export
  kh_shim.hpp / kh_shim.cpp C++ wrapper: calls kh_simulate_c(), packages SimResult
  main.cpp                  cpp-httplib HTTP server on :8004
tests/
  validation_test.cpp       standalone binary vs golden vector V1 (A2 band)
```

## Language architecture

Physics is entirely in Fortran (`kh_physics.f90`). The `kh_simulate_c` subroutine
is exported via `bind(c)` / `iso_c_binding`, making it callable from C/C++ without
any name-mangling issues. The C++ shim (`kh_shim.cpp`) calls this symbol directly,
handles memory allocation, and packages the result for the HTTP layer. cpp-httplib
and nlohmann/json handle HTTP and serialisation (same as KH-005 C++).

## Build and run

Prerequisites: CMake 3.20+, gfortran 11+, g++ 11+, Ninja (optional; any CMake generator works).

```powershell
cd kh-sim\backends\fortran

cmake -B build -DCMAKE_BUILD_TYPE=Release -G Ninja `
      -DCMAKE_Fortran_COMPILER=gfortran `
      -DCMAKE_CXX_COMPILER=g++

cmake --build build

.\build\kh-sim-fortran.exe
```

First configure downloads cpp-httplib + nlohmann/json (same cache as KH-005 if shared).

## Validation

```powershell
.\build\kh-sim-fortran-validate.exe ..\..\shared\golden-vectors\v1-canonical-reference.json
```

Expected:
```
--- Fortran vs Python reference ---
kinetic_energy : got=0.112810  ref=0.112810  err=0.00%
enstrophy      : got=43.705142  ref=43.705142  err=0.00%
divergence_rms : got=~1e-14    ref=1.20e-14
All tests PASSED
```

## Status

v0.1.0 (2026-10), validated against the Python reference oracle on every layer of the suite:

- golden vectors v1–v9 of `shared/golden-vectors/`, element by element (tolerance factor 1×);
- the physics testbeds (conservation, decay, linear growth) of `tests/integration/`;
- the canonical textbook testbeds CP-7F-A…D, run once for all six implementations from the shared harnesses in `tests/integration/test_canonical_physics_cp7f*.py`; CP-7F-E (the Pascal-binding anchor) is asserted natively in the Fortran test suite as well.
- request validation and response schema per `shared/api/openapi.yaml` (`tests/integration/test_conformance.py`).

Container: `docker build -t kh-sim-fortran .` then `docker run --rm -p 8004:8004 kh-sim-fortran`.
