# KH-SIM — Rust Backend

**Task:** KH-003 | **Framework:** axum 0.7 | **Port:** 8001 | **FFT:** rustfft 6

## Architecture

```
src/
  lib.rs          -- library root (re-exported by integration tests)
  main.rs         -- axum server, binds 0.0.0.0:8001
  models.rs       -- SimulationRequest / SimulationResult (mirrors openapi.yaml)
  handlers.rs     -- POST /simulate, GET /health, GET /info
  physics/
    mod.rs        -- pub use solver::simulate
    fft2d.rs      -- 2D full-complex FFT (Fft2DPlans, make_wavenumbers)
    solver.rs     -- IC, Poisson solve, velocity recovery, RK4 time loop
tests/
  validation_test.rs -- acceptance test vs golden vector V1 (A2 band)
```

## Endpoint contract

Implements `kh-sim/shared/api/openapi.yaml`

| Endpoint  | Method | Description                                   |
|-----------|--------|-----------------------------------------------|
| /simulate | POST   | Run KH simulation, return SimulationResult    |
| /health   | GET    | Liveness probe                                |
| /info     | GET    | Backend language/framework metadata           |

## Physics kernel

Exact port of `kh-sim/shared/physics/kh_physics.py`:

- Initial conditions: tanh shear layer + sinusoidal perturbation (analytic vorticity)
- Poisson solver: spectral (psi_hat = omega_hat / k^2), zero mean mode
- Velocity recovery: u = i*ky*psi_hat, v = -i*kx*psi_hat (spectral differentiation)
- Time integrator: RK4 (4 stages per step)
- FFT: 2D full-complex DFT via rustfft, axis-0 fftfreq(nx), axis-1 fftfreq(ny)

**Note on kh-instability-sim.zip (uploaded 2026-03-28):**
The zip contains a browser-based D2Q9 Lattice Boltzmann simulation (JavaScript/Canvas,
single HTML file -- a different physics formulation and architecture). It has been
preserved at `kh-sim/reference/kh-instability-sim.zip`. The Rust backend implements
the canonical Navier-Stokes vorticity-streamfunction method per `KH-PHYSICS.md`.

## Build and run

```powershell
# Prerequisites: Rust stable toolchain (rustup), cargo in PATH
cd kh-sim\backends\rust

# Debug build (fast compile, slower run)
cargo build

# Release build (optimised -- use for benchmarking)
cargo build --release

# Run (binds 0.0.0.0:8001)
cargo run --release
```

Environment:
- `RUST_LOG=debug`  -- verbose request/response tracing
- `RUST_LOG=info`   -- default, operational logs only

## Validation test

```powershell
# Prerequisite: generate reference output first
cd kh-sim\shared\physics
python kh_physics.py

# Run acceptance tests
cd ..\..\backends\rust
cargo test --test validation_test -- --nocapture
```

Expected output:
```
--- Rust vs Python reference ---
kinetic_energy : got=0.112810  ref=0.112810  err=~0%
enstrophy      : got=43.7051   ref=43.7051   err=~0%
max_vorticity  : got=31.9076   ref=31.9076   err=~0%
divergence_rms : got=~1e-14    ref=1.20e-14
```

Acceptance criteria (KH-PHYSICS.md S7):
- kinetic_energy within +-5% of reference
- enstrophy within +-5% of reference
- divergence_rms < 1e-10

## Smoke test (server running)

```powershell
# Health check
curl http://localhost:8001/health

# Simulation request (small grid, fast)
$body = '{"grid_nx":64,"grid_ny":32,"steps":100,"reynolds_number":1000}'
Invoke-RestMethod -Method POST -Uri http://localhost:8001/simulate `
    -ContentType "application/json" -Body $body | ConvertTo-Json -Depth 5

```

## Status

v0.1.0 (2026-10), validated against the Python reference oracle on every layer of the suite:

- golden vectors v1–v9 of `shared/golden-vectors/`, element by element (tolerance factor 1×);
- the physics testbeds (conservation, decay, linear growth) of `tests/integration/`;
- the canonical textbook testbeds CP-7F-A…D, run once for all six implementations from the shared harnesses in `tests/integration/test_canonical_physics_cp7f*.py`;
- request validation and response schema per `shared/api/openapi.yaml` (`tests/integration/test_conformance.py`).

Container: `docker build -t kh-sim-rust .` then `docker run --rm -p 8001:8001 kh-sim-rust`.
