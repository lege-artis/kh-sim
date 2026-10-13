# lege-artis/kh-sim — engineer tier — 00: Quick start (EN)

> **Audience.** Engineer, scientist, data scientist, hobbyist. You want to run kh-sim, read what it returns, and know whether the result can be trusted. You need a working mental model and runnable code, not the derivations.
>
> **Reading time.** About 10 minutes for this chapter, about 40 minutes for the engineer tier (00–02).
>
> **License.** CC-BY-SA-4.0 (text); the code snippets are Apache-2.0.
>
> **Companion.** The canonical tier, [`shared/canonical-equations/`](../../../shared/canonical-equations/), states the same scope as equations with their derivations.

---

## §1. What this library does

In one sentence: **kh-sim integrates a two-dimensional incompressible flow forward in time, starting from a shear layer, and returns the flow field and its diagnostics at the end time.**

| You give it | kh-sim does | You get back |
|---|---|---|
| A request: grid, box size, time step, number of steps, Reynolds number, the ripple that seeds the shear layer (or your own initial vorticity field) | Integrates the 2D Navier–Stokes equations in vorticity–stream-function form: pseudo-spectral in space, classical RK4 in time, 2/3-rule de-aliasing | The vorticity, both velocity components and the stream function on the grid at `t_final`, plus five diagnostics and any resolution warnings |
| The same request, as JSON, to any of five HTTP backends | The same computation in Rust, Scala, C++, Fortran or Pascal | The same response, field for field, within a published tolerance |
| The field a previous run returned, with `t_start` | Continues that run | The next segment of the same trajectory |

The Python module [`shared/physics/kh_physics.py`](../../../shared/physics/kh_physics.py) is the **reference oracle**: the five other implementations are validated against it, element by element.

## §2. The intuition

Two layers of fluid slide past each other, in kh-sim's convention the upper one to the left and the lower one to the right. The interface between them is a thin band in which the velocity changes sign; that band is the **shear layer**. Give the interface a small ripple and the shear amplifies it: the ripple grows exponentially, steepens, and rolls up into a row of vortices — the **Kelvin–Helmholtz billows** familiar from cloud tops and stratified water. At high enough Reynolds number the billows go on to interact and merge.

kh-sim computes exactly this, in a box that is periodic in both directions. Chapter [01](01-what-kh-sim-computes.md) says what is being computed and which questions the solver can answer; chapter [02](02-numerical-methods.md) says how.

## §3. Hello-World example (Python reference)

Save the snippet as `hello_kh.py` in the repository root (the directory that holds `shared/`) and run `python hello_kh.py`. It needs Python 3.10+ and NumPy, and runs in about two seconds.

<!-- hello-world -->
```python
from shared.physics.kh_physics import SimulationRequest, simulate

# The canonical shear layer: 128 x 64 grid on a 1.0 x 0.5 box, mode-2 ripple.
req = SimulationRequest(grid_nx=128, grid_ny=64, dt=0.001, steps=500,
                        reynolds_number=1000.0, perturbation_amplitude=0.01,
                        perturbation_mode=2)
res = simulate(req)

d = res.diagnostics
print(f"t_final        = {res.t_final:.3f}")
print(f"kinetic_energy = {d.kinetic_energy:.6f}")
print(f"enstrophy      = {d.enstrophy:.4f}")
print(f"max_vorticity  = {d.max_vorticity:.4f}")
print(f"cfl_peak       = {d.cfl_peak:.4f}")
print(f"divergence ok  = {d.divergence_rms < 1e-10}")
print(f"warnings       = {res.warnings}")
print(f"field length   = {len(res.vorticity)}  (= 128 * 64, row-major [ix*ny + iy])")
```

**Expected output** (the test suite runs this snippet and compares every line):

<!-- hello-world-output -->
```text
t_final        = 0.500
kinetic_energy = 0.091226
enstrophy      = 29.2264
max_vorticity  = 21.4639
cfl_peak       = 0.0994
divergence ok  = True
warnings       = []
field length   = 8192  (= 128 * 64, row-major [ix*ny + iy])
```

Three lines are worth checking by eye. `divergence ok = True`: the velocity is divergence-free to rounding error, as it must be for incompressible flow. `cfl_peak = 0.0994`: the time step is ten times below the advective stability limit. `warnings = []`: nothing in the request was computed but possibly not resolved (§4).

If your numbers differ in the fourth significant digit, suspect the request before the solver: a different `dt`, `steps` or grid changes every line.

## §4. Reading the output

`simulate()` returns a `SimulationResult`; an HTTP backend returns the same object as JSON (schema: [`shared/api/openapi.yaml`](../../../shared/api/openapi.yaml)).

| Field | Meaning |
|---|---|
| `vorticity` | ω(x, y) at `t_final`, `grid_nx * grid_ny` values |
| `u_velocity`, `v_velocity` | the velocity components, same layout |
| `pressure` | **the stream function ψ** — the field keeps its historical name; kh-sim never computes a pressure |
| `t_final` | `t_start + steps * dt` |
| `diagnostics.kinetic_energy` | ½ · mean(u² + v²) |
| `diagnostics.enstrophy` | ½ · mean(ω²) |
| `diagnostics.max_vorticity` | max \|ω\| |
| `diagnostics.divergence_rms` | RMS of ∂u/∂x + ∂v/∂y; about 1e-14 when healthy |
| `diagnostics.cfl_peak` | the largest Courant number dt · (max\|u\|/dx + max\|v\|/dy) seen during the run |
| `warnings` | computed, but possibly not resolved: `cfl_advective:` (Courant number above 1), `cfl_diffusive:` (2νdt/min(dx,dy)² above 1), `inviscid:` (ν = 0) |

Every field is flattened **row-major with x as the outer index**: the value at grid point (ix, iy) sits at `ix * grid_ny + iy`. With NumPy, `np.array(res.vorticity).reshape(grid_nx, grid_ny)` gives an array indexed `[ix, iy]`; transpose it for `imshow` with x horizontal.

To take snapshots of one run, continue it: pass the returned `vorticity` as `initial_omega`, the returned `t_final` as `t_start`, and `perturbation_amplitude=0` in the next request. A continued run equals the unsplit run. The HTTP backends refuse a non-zero `perturbation_amplitude` sent together with `initial_omega`, because the field replaces the analytic initial condition and the ripple would be silently ignored. The figure scripts of this tier continue their runs exactly this way (chapter [01](01-what-kh-sim-computes.md)).

## §5. Which backend to use when

All six implementations compute the same thing. They differ in what they are good for.

| Situation | Use | Why |
|---|---|---|
| Learning the algorithm, small grids, scripting experiments | **Python reference** | Shortest path from equation to code; it is the oracle the others are checked against |
| Larger grids, many runs, a service | **Rust** (rustfft) or **Scala** (JTransforms) | Library-backed FFTs, compiled speed |
| Seeing every arithmetic step of the FFT | **C++**, **Fortran** or **Pascal** | Hand-rolled radix-2 Cooley–Tukey FFT; grid sizes must be powers of two |
| Reading the numerics next to a textbook | **Pascal** | Bound to the Pascal edition of *Numerical Recipes* (1986); see chapter [02](02-numerical-methods.md) §4 |

Each HTTP backend is built from its own `Dockerfile` and serves `POST /simulate`, `GET /health` and `GET /info` on a fixed port:

| Backend | Directory | Port |
|---|---|---|
| Rust | [`backends/rust/`](../../../backends/rust/) | 8001 |
| Scala | [`backends/scala/`](../../../backends/scala/) | 8002 |
| C++ | [`backends/cpp/`](../../../backends/cpp/) | 8003 |
| Fortran | [`backends/fortran/`](../../../backends/fortran/) | 8004 |
| Pascal | [`backends/pascal/`](../../../backends/pascal/) | 8005 |

For example, Rust, from the repository root:

```sh
docker build -t kh-sim-rust backends/rust
docker run --rm -p 8001:8001 kh-sim-rust
curl -s -X POST http://localhost:8001/simulate -H "Content-Type: application/json" \
     -d '{"grid_nx":128,"grid_ny":64,"dt":0.001,"steps":500,"reynolds_number":1000,
          "perturbation_amplitude":0.01,"perturbation_mode":2}'
```

The request fields are the fields of `SimulationRequest`; `GET /info` reports the backend's grid constraint and FFT library. Native builds (CMake, Cargo, sbt, FPC) are described in each backend's `README.md`.

## §6. The one habit that matters

Compare fields **element-wise**, not by aggregate. A Kelvin–Helmholtz run's total energy does not change when an FFT sign convention is flipped, so an energy check passes a solver whose field is the mirror image of the right one. Every cross-backend defect found in kh-sim's history was caught by an element-wise comparison and would have survived a comparison of norms. The test architecture (chapter [02](02-numerical-methods.md) §5) is built on that lesson; if you compare two runs of your own, do the same: `np.max(np.abs(a - b))`, not `abs(a.sum() - b.sum())`.

## §7. Where to next

- **What is being computed:** [01 — What kh-sim actually computes](01-what-kh-sim-computes.md) — the equations in engineer's terms, the billows, what you can and cannot ask of the solver, how to read the diagnostics.
- **How it is computed:** [02 — Numerical methods](02-numerical-methods.md) — pseudo-spectral derivatives, RK4, de-aliasing, and why there are six implementations.
- **The mathematics:** [`shared/canonical-equations/kh-governing-equations.md`](../../../shared/canonical-equations/kh-governing-equations.md) (KH-GE).
- **The story:** the shad tier, [*Just Shad's Guide to the Kelvin–Helmholtz Galaxy*](../../shad/00-prologue.md).

---

*Next:* [01 — What kh-sim actually computes](01-what-kh-sim-computes.md)
