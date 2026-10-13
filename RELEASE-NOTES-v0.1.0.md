# lege-artis/kh-sim — v0.1.0

**Tag:** `v0.1.0`
**Status:** First public release. Six independent implementations of one solver, validated element by element against each other and against textbook physics.

---

## What this is

A canonical reference implementation of two-dimensional incompressible Navier–Stokes flow in vorticity–stream-function form, specialised to the Kelvin–Helmholtz shear instability: pseudo-spectral in space, classical RK4 in time, 2/3-rule de-aliasing, doubly periodic box.

The same equations — stated once in `shared/canonical-equations/` (KH-GE, KH-NM, KH-LS) — are implemented six times:

| Implementation | FFT | Role |
|---|---|---|
| Python / NumPy (`shared/physics/kh_physics.py`) | NumPy | the reference oracle |
| Fortran (`backends/fortran`) | hand-rolled radix-2 | the canonical compiled implementation |
| Pascal (`backends/pascal`) | hand-rolled radix-2 | the pedagogical one, bound to *Numerical Recipes in Pascal* (1986) |
| C++ (`backends/cpp`) | hand-rolled radix-2 | every arithmetic step visible |
| Rust (`backends/rust`) | rustfft | independent library cross-check, HTTP service |
| Scala (`backends/scala`) | JTransforms | independent library cross-check, HTTP service |

The five compiled implementations serve one HTTP contract (`shared/api/openapi.yaml`: `POST /simulate`, `GET /health`, `GET /info`), with identical request validation and resolution warnings.

## Validation

| Layer | What | Result |
|---|---|---|
| Golden vectors | Nine reference runs `shared/golden-vectors/v1…v9`, diagnostics and field elements, against the oracle within a √N·ε rounding bound | All six; tolerance factor 1× (Python, Rust, C++, Fortran), 2× (Scala), 10× (Pascal) |
| Contract | Request validation, response schema, `/health`, `/info` (`tests/integration/test_conformance.py`) | All five HTTP backends |
| Physics testbeds | Conservation (total vorticity, inviscid energy and enstrophy), viscous energy decay, incompressibility, Taylor–Green decay | All six |
| Canonical textbook testbeds | A: dispersion σ(k) of the tanh layer (Michalke 1964) · B: forced 2D-turbulence spectrum · C: co-rotating vortex merger (Melander, Zabusky & McWilliams 1988) · D: forced Taylor–Green fixed point · E: *Numerical Recipes in Pascal* §17.5 and §12.4 | A: σ_∞ = 0.18818 / 0.18745 at kδ = 0.4 / 0.5 against Michalke's 0.1881 / 0.1875, all six identical to the printed digit · B: slopes −4.12 … −4.16 (finite-Re coherent-vortex regime) · C: merge at 2.5, 3.0, 3.3 radii, none at 3.5 and 4.0, all six · D: within 1e-6 at five Reynolds numbers · E: Python, Pascal, Fortran reproduce the printed results |
| Documentation | Every code snippet and example script of the docs is run, and its output compared with what the chapter prints (`tests/docs/`) | 9/9 |

The textbook testbeds A–D run once, from shared harnesses in `tests/integration/test_canonical_physics_cp7f*.py`, against all six implementations through the same request API users call.

## What's in this release

- **Code** (Apache 2.0): the Python oracle and the five backends, each with its build (CMake, Cargo, sbt, Make/FPC) and Dockerfile; tools for golden-vector generation (`tools/`).
- **Tests**: `tests/integration/` (golden vectors, contract, physics and canonical testbeds; the HTTP suites skip when no backend is running), `tests/docs/`, and each backend's native tests.
- **Documentation** (CC-BY-SA-4.0), three tiers:
  - canonical — `shared/canonical-equations/` (KH-GE v0.1.1, KH-NM, KH-LS; Markdown + LaTeX);
  - engineer — `docs/engineer/en/` 00 Quick start, 01 What kh-sim computes, 02 Numerical methods;
  - shad — `docs/shad/` 00 Prologue, 01 The billows, with `examples/shad/c1-shear-layers/`.
- **CI**: six GitHub Actions workflows, one per implementation.
- **Governance**: `LICENSE` (Apache 2.0), `LICENSE-DOCS` (CC-BY-SA-4.0), `NOTICE`, `TRADEMARK.md`, `CONTRIBUTING.md`, `CODE_OF_CONDUCT.md`, `SECURITY.md`.

## Known limits of v0.1.0

- **Model:** two-dimensional, incompressible, one density, doubly periodic. The analytic initial condition fixes the layer width at L_y/20; a single layer in a periodic box is realised with a uniform background vorticity of the opposite sign, which lowers the growth rate of the default box to about three quarters of the unbounded layer's (engineer tier 01 §2, §5). The textbook comparison A uses a box built for it.
- **ETDRK4** is specified (KH-NM-3b) and present in the Fortran library with its unit tests, but no request uses it.
- **CI** runs each implementation's own build and tests; the HTTP suites of `tests/integration/` need running backends and skip in CI. The full cross-implementation run (all six, including the hours-long CP-7F-A and B) is made on a reference workstation before each release.
- **Development record:** some code comments cite design and audit documents of the development history (`_audit/…`, `_handoffs/…`); those documents are not part of this repository.

## Next (v0.1.x)

Shad tier chapters 2–5 (vorticity and stream function, reading the growth rate, why six languages, shear layers in the wild); translations of the documentation; the asymptotic k⁻³ regime of testbed B at higher Reynolds number.

## Sibling projects

- [lege-artis/fourier](https://github.com/lege-artis/fourier) — the DFT/FFT canonical reference; same multi-implementation discipline.
- [lege-artis/jwst-l2-coupled-dynamics](https://github.com/lege-artis/jwst-l2-coupled-dynamics) — two-rigid-body coupled dynamics.
