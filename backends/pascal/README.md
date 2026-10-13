# KH-SIM — Pascal Backend

**Framework:** fphttpserver (Free Pascal fcl-web) | **Port:** 8005 | **FFT:** radix-2 Cooley–Tukey, hand-rolled | **Compiler:** FPC 3.2.2

The pedagogical reference among the six implementations: the solver reads like the textbook it is bound to. Its numerics reproduce worked examples of *Numerical Recipes in Pascal* (Press et al., 1986): the Poisson problem of §17.5 and the convolution of §12.4 (canonical testbed CP-7F-E), so the code can be followed with the book open.

## Architecture

```
src/
  kh_physics.pas      the solver: FFT1D (radix-2), 2D transforms, Poisson solve,
                      2/3-rule de-aliasing, RK4, diagnostics, K5 warnings
  kh_forcing.pas      steady and rotating-phase forcing ("mode", "band")
  kh_server.pas       HTTP layer: POST /simulate, GET /health, GET /info;
                      request validation KH-VAL-1/KH-VAL-2 (400 {"error","field"}),
                      500 {"error","backend"} when the computation fails
tests/
  kh_validate.pas            golden vector v1 against the Python oracle (factor 10x)
  test_canonical_physics.pas canonical testbeds: A-D print where they run (shared
                             harnesses), E asserts NR §17.5 and §12.4 (kh_cp7fe.pas)
  test_solver_guards.pas     solver guards the HTTP layer cannot reach (FFT length,
                             reported steps and time)
```

The physics is a line-by-line implementation of `shared/canonical-equations/` (KH-GE, KH-NM); the request and response follow `shared/api/openapi.yaml`, identically to the other four HTTP backends.

## Build and run

POSIX (Makefile):

```sh
make server           # build/kh_server
make validate         # build/kh_validate
make canonical        # build/test_canonical_physics
./build/kh_server     # listens on :8005
```

Windows (PowerShell): `.\build.ps1 [all|server|validate|canonical|clean]`.

Container: `docker build -t kh-sim-pascal .` then `docker run --rm -p 8005:8005 kh-sim-pascal`.

## Tests

```sh
./build/kh_validate ../../shared/golden-vectors/v1-canonical-reference.json
./build/test_canonical_physics        # RESULT: PASS (CP-7F-E: PASS (E2 + E1))
mkdir -p build/guards && fpc -O2 -Fusrc -FUbuild/guards -FEbuild/guards tests/test_solver_guards.pas && ./build/guards/test_solver_guards
```

The HTTP-level suites (golden vectors v1–v9, physics testbeds, CP-7F-A…D, conformance) run from `tests/integration/` against a running server.

## Status

v0.1.0 (2026-10): golden vectors v1–v9 element by element against the Python oracle (tolerance factor 10×: manual FFT loops in FPC `double` accumulate rounding differently); the physics testbeds; CP-7F-A…D through the shared harnesses; CP-7F-E natively, the strict textbook-binding test of this backend; conformance to `shared/api/openapi.yaml`.
