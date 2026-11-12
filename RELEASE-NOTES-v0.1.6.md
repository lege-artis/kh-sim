# lege-artis/kh-sim — v0.1.6

**Tag:** `v0.1.6`
**Type:** documentation. Three new chapters of the shad tier, *Just Shad's Guide to the Kelvin–Helmholtz Galaxy*, each with a runnable script. No solver, backend or contract change.

## The chapters

| Chapter | What it does | What it shows |
|---|---|---|
| [02 — The vortex and the stream](docs/shad/02-the-vortex-and-the-stream.md) | One Lamb–Oseen vortex left to decay; two co-rotating vortices orbiting each other (128², Re = 1000) | The peak vorticity follows Γ/(π(a² + 4νt)) to every printed digit for ten time units. The pair orbits at 1.282 against 1.389 for point vortices in the plane: the box holds no net spin, which removes Γ/L² (prediction 1.256, measured 2 % above it) |
| [03 — Reading the growth rate](docs/shad/03-reading-the-growth-rate.md) | Eleven ripple sizes kδ = 0.1 … 1.1 in a box of two mirror layers (128 × 256, Re = 40 000), each seeded with the box's own unstable mode | Up to kδ = 0.8 the measured growth rates lie 0.001–0.002 below the box's inviscid modes; at kδ = 0.4 the rate is linear in 1/Re and extrapolates to 0.1883, the box's own number, against Michalke's 0.1881. Near the cut-off (kδ = 0.9) the grid, not the physics, sets the error |
| [04 — Why six languages](docs/shad/04-why-six-languages.md) | Chapter 1's billows run with two bugs kh-sim really had, put back into the reference: the conjugate FFT convention (once in C++, Fortran, Pascal) and the missing 2/3 rule (once in Rust) | Under the conjugate FFT every summary and the vorticity are bit-identical and the velocity is exactly reversed. The missing 2/3 rule moves the summaries by 4e-6 to 8e-3 at t = 1.2, and by less than 1e-3 at t = 0.3 or at golden-vector V3 conditions. Only element-wise comparison against the reference, at a rounding tolerance, sees both |

The prologue, chapter 1 and the tier's README link the new chapters; the engineer tier's chapter 02 links chapter 4.

## The scripts

`examples/shad/c2-vortex-and-stream/main.py`, `c3-growth-rate/main.py`, `c4-six-languages/main.py` — Apache-2.0, ASCII-only, computed with the Python reference `shared/physics/kh_physics.py`. Each writes its figures to `docs/shad/figures/` (or `--outdir`). Runtimes on a reference machine: about 50 s, 2 min 15 s and 16 s. Chapter 3 also needs SciPy (already a project dependency) for the Rayleigh eigenvalue problem.

## Tests

`tests/docs/test_doc_tiers.py` now covers chapters 02–04: their links, licence lines and figures (12 figures, each at most 300 KB), the ASCII rule for the new scripts, and one slow test per chapter that runs its script and compares the printout with the chapter's block line by line. 12 tests, about 4.5 minutes; `-m "not slow"` skips the script runs.
