# KH-SIM — Scala Backend

**Task:** KH-004 | **Framework:** http4s 0.23 (Ember) | **Port:** 8002 | **FFT:** JTransforms 3.1

## Architecture

```
build.sbt                           sbt build (Scala 3.3.4, http4s 0.23, JTransforms 3.1)
project/
  build.properties                  sbt 1.10.6
  plugins.sbt                       sbt-native-packager
src/main/scala/khsim/
  Main.scala                        IOApp.Simple: EmberServer on :8002 + CORS
  Models.scala                      SimulationRequest / SimulationResult (circe codecs)
  Routes.scala                      HttpRoutes[IO]: POST /simulate, GET /health, GET /info
  physics/
    Fft2D.scala                     DoubleFFT_2D wrapper + Wavenumbers.make()
    Solver.scala                    IC, Poisson, velocity recovery, RK4 time loop
src/test/scala/khsim/
  ValidationTest.scala              ScalaTest vs golden vector V1 (A2 band)
```

## Endpoint contract

Implements `kh-sim/shared/api/openapi.yaml`

| Endpoint  | Method | Description                                |
|-----------|--------|--------------------------------------------|
| /simulate | POST   | Run KH simulation, return SimulationResult |
| /health   | GET    | Liveness probe                             |
| /info     | GET    | Backend language/framework metadata        |

## Physics kernel

Exact port of `kh-sim/shared/physics/kh_physics.py` -- same equations, same algorithm:
- IC: tanh shear layer + sinusoidal perturbation (analytic vorticity)
- Poisson: psi_hat = omega_hat / k^2, zero mean mode
- Velocity: u = i*ky*psi_hat, v = -i*kx*psi_hat (spectral differentiation)
- Time integrator: RK4, 4 stages per step
- FFT: JTransforms DoubleFFT_2D (full complex, both axes)

## Build and run

Prerequisites: JDK 17+ and sbt 1.10.x in PATH.

```powershell
cd kh-sim\backends\scala

# Download dependencies + compile (first run ~2-3 min)
sbt compile

# Run validation tests
sbt test

# Run server (binds 0.0.0.0:8002)
sbt run

# Production binary distribution
sbt stage
.\target\universal\stage\bin\kh-sim-scala.bat
```

## Validation test

```powershell
# Prerequisite: generate reference output
python kh-sim\shared\physics\kh_physics.py

sbt test
```

Expected output:
```
--- Scala vs Python reference ---
kinetic_energy : got=0.112810  ref=0.112810  err=~0%
enstrophy      : got=43.7051   ref=43.7051   err=~0%
max_vorticity  : got=31.9076   ref=31.9076   err=~0%
divergence_rms : got=~1e-14    ref=1.20e-14
```

## Smoke test (server running)

```powershell
curl http://localhost:8002/health

$body = '{"grid_nx":64,"grid_ny":32,"steps":100,"reynolds_number":1000}'
Invoke-RestMethod -Method POST -Uri http://localhost:8002/simulate `
    -ContentType "application/json" -Body $body | ConvertTo-Json -Depth 5

```

## Status

v0.1.0 (2026-10), validated against the Python reference oracle on every layer of the suite:

- golden vectors v1–v9 of `shared/golden-vectors/`, element by element (tolerance factor 2×);
- the physics testbeds (conservation, decay, linear growth) of `tests/integration/`;
- the canonical textbook testbeds CP-7F-A…D, run once for all six implementations from the shared harnesses in `tests/integration/test_canonical_physics_cp7f*.py`;
- request validation and response schema per `shared/api/openapi.yaml` (`tests/integration/test_conformance.py`).

Container: `docker build -t kh-sim-scala .` then `docker run --rm -p 8002:8002 kh-sim-scala`.
