# lege-artis/kh-sim — engineer tier — 02: Numerical methods (EN)

> **Audience.** You know what kh-sim computes (chapter [01](01-what-kh-sim-computes.md)) and want to know how: what turns the vorticity equation into arithmetic, where the method can fail, and how the test suite makes sure it does not.
>
> **Reading time.** About 15 minutes.
>
> **License.** CC-BY-SA-4.0.
>
> **Companion.** Canonical tier: [`kh-numerical-methods.md`](../../../shared/canonical-equations/kh-numerical-methods.md) (KH-NM). This chapter cites its equation labels.

---

## §1. The big picture

Three decisions define the solver, and each has a testable consequence.

**Space: pseudo-spectral.** The box is doubly periodic and the fields are smooth, which is exactly where Fourier methods are at their best. Every spatial derivative is computed by transforming the field to Fourier space, multiplying by the wavenumber, and transforming back. Derivatives are then exact for every resolved mode — no finite-difference truncation error — and the Poisson equation for the stream function becomes one division per mode. The nonlinear advection term is the "pseudo" part: it is multiplied out on the grid, not in Fourier space.

**Time: classical fourth-order Runge–Kutta (RK4).** The right-hand side — advection plus viscous diffusion — is advanced with RK4, an explicit method: fourth-order accurate, four evaluations of the right-hand side per step, and the same few lines in every language. That last property matters when six implementations must agree element by element.

**Aliasing: the 2/3-rule.** The product in the advection term creates wavenumbers the grid cannot hold; the top third of the spectrum is zeroed after every evaluation so that they cannot fold back onto the resolved modes.

One step of the solver, from ω to the next ω, is therefore: FFT, divide by k² (ψ̂), multiply by i k (û, v̂, ∂ω/∂x, ∂ω/∂y), inverse FFTs, multiply on the grid, FFT, cut the top third, add the viscous term −ν k² ω̂ — four times, once per RK4 stage.

## §2. Spectral derivatives via the FFT

**Why the FFT.** A derivative in Fourier space costs one multiplication per mode; the transforms cost O(N log N) for N grid points with the fast Fourier transform, against O(N²) for a direct transform. A 256 × 128 grid has 32 768 points: the FFT is the difference between a step that takes milliseconds and one that takes seconds.

**The wavenumbers (KH-NM-1).** For a grid of N_x × N_y points on L_x × L_y:

```text
k_x = (2π / L_x) · {0, 1, …, N_x/2 − 1, −N_x/2, …, −1}      full transform along x
k_y = (2π / L_y) · {0, 1, …, N_y/2}                          real transform along y
```

The field is real, so its transform is conjugate-symmetric and only the non-negative k_y half is stored (`rfft2`). Then

```text
∂/∂x → i k_x ·      ∂/∂y → i k_y ·      ∇² → −(k_x² + k_y²) ·
ψ̂ = ω̂ / (k_x² + k_y²)   for (k_x, k_y) ≠ (0, 0),     ψ̂(0, 0) = 0     (KH-NM-2)
```

**The convention.** The forward transform is `X[k] = Σ x[n] exp(−2πi k n / N)`, unnormalised; the inverse carries the 1/N. All six implementations use this convention, and the test suite checks it field by field, because a transform with the opposite sign is the complex conjugate of the right one: it differentiates with the wrong sign, and the energy and enstrophy totals do not notice. kh-sim has met this defect in practice: it hid behind matching energy totals until the fields were compared element by element. That history is why every comparison in the suite is element-wise (§5).

The FFT algorithms are standard: Press et al., *Numerical Recipes* (3rd ed., 2007), §12.1–12.2 for the radix-2 Cooley–Tukey transform that the C++, Fortran and Pascal implementations write out by hand; Python, Rust and Scala call NumPy's, rustfft's and JTransforms' transforms.

**The 2/3-rule (KH-NM-4; Orszag, 1971).** The advection term u ∂ω/∂x + v ∂ω/∂y multiplies two fields. If each holds wavenumbers up to k_max, their product holds wavenumbers up to 2 k_max, and on a periodic grid the part above k_max does not disappear: it **aliases** onto lower wavenumbers as spurious content. In a roll-up or a turbulent run that content piles up at the smallest resolved scales and can make the run blow up. The remedy is to zero, after every evaluation of the product, every mode with |k_x| above two thirds of the largest k_x or k_y above two thirds of the largest k_y. The aliased products then land only on the zeroed band. The price is a third of the resolution in each direction; the return is a clean, stable cascade.

Two consequences you will meet in requests: a ripple or a forcing mode beyond the 2/3 band would be cut at the first step, so the HTTP backends refuse it (`perturbation_mode` with 3·m > 2·(N_x/2)); and the effective resolution of a 256² run is about 85 modes per direction, not 128.

## §3. RK4 time stepping and its limits

Write the vorticity equation as ∂ω/∂t = f(ω), with f the advection plus the viscous term, the velocity recovered from ω at every stage. One RK4 step (KH-NM-3a; Press et al. 2007, §17.1):

```text
k1 = f(ω_n)
k2 = f(ω_n + Δt/2 · k1)
k3 = f(ω_n + Δt/2 · k2)
k4 = f(ω_n + Δt · k3)
ω_{n+1} = ω_n + Δt/6 · (k1 + 2 k2 + 2 k3 + k4)
```

RK4 is explicit, so the step size is bounded (KH-NM-5):

| Bound | Condition | What it limits |
|---|---|---|
| Advective (CFL) | Δt · (max\|u\|/Δx + max\|v\|/Δy) below about 1 | The flow must not cross more than about a cell per step |
| Diffusive | 2 ν Δt / min(Δx, Δy)² below about 1 | Viscosity must not smooth faster than the step resolves |

kh-sim does **not** adapt Δt. It computes what was asked, measures the peak Courant number during the run (`cfl_peak`), and returns a warning — `cfl_advective:` or `cfl_diffusive:` — when a bound is exceeded. The decision to trust or discard such a run is yours. For the canonical layer the advective bound is the active one: the Hello-World of chapter [00](00-quick-start.md) runs at a peak Courant number of 0.10, the 256 × 128 roll-up of chapter [01](01-what-kh-sim-computes.md) at about 0.5.

RK4 is the canonical integrator of every implementation and of every request. The canonical tier also specifies an exponential integrator, ETDRK4 (KH-NM-3b; Cox and Matthews, 2002), which treats the viscous term exactly and relaxes the diffusive bound for stiff, high-resolution runs. It exists as an alternate scheme in the Fortran library, with its own unit tests (linear decay, convergence order); no request uses it, and no cross-implementation result depends on it.

## §4. Why six implementations

kh-sim implements the same equations six times, on purpose. The short version:

- **The Python reference is the oracle.** It is the shortest path from equation to code; every other implementation is compared against its fields.
- **C++, Fortran and Pascal write the FFT out by hand.** Every arithmetic step of the transform is visible and reviewable. Fortran is the canonical compiled implementation: the most extensively documented and tested of the three. Pascal is bound to a textbook: its numerics reproduce worked examples of the Pascal edition of *Numerical Recipes* (1986) — the Poisson solution of §17.5 and the convolution of §12.4 — so that a reader can follow the code with the book open.
- **Rust and Scala call third-party FFT libraries** (rustfft, JTransforms) that share no code with the hand-rolled transforms or with NumPy. When they agree with the oracle element by element, the agreement is evidence about the equations, not about one particular FFT. They also serve the same API as typed, compiled services.

If six independent implementations, three hand-rolled and three library-backed, produce the same field to within rounding, a shared mistake would have to be a mistake in the equations, and the equations are checked against textbook physics (§5). That is the claim "canonical reference" makes, and the reason the claim is testable.

## §5. Validation discipline

The suite checks the solver in layers, each with a number to hit and a source that could have shown it wrong.

| Layer | What is compared | Against what |
|---|---|---|
| Unit | FFT round trip, spectral derivatives, Poisson solve, de-aliasing, CFL bound | Analytic fields, at machine precision |
| Golden vectors | Nine reference runs (`shared/golden-vectors/v1…v9`): diagnostics and field elements | The oracle, element by element, within a rounding bound that grows as √N · machine ε (Higham's accumulation bound); the bound's per-backend factor (1× to 10×) is stored in the vector, never in the test |
| Contract | Request validation, response schema, `/health`, `/info` | [`shared/api/openapi.yaml`](../../../shared/api/openapi.yaml), identically for the five HTTP backends |
| Physics testbeds | Conservation, viscous energy decay (dE/dt = −2νZ with the mean mode removed), Taylor–Green decay | The analytic identities of the equations |
| Textbook physics | Five canonical testbeds, run once for all six implementations | Published results (below) |

The textbook layer, as of v0.1.0:

| Testbed | Physics | Result, all six implementations |
|---|---|---|
| A — dispersion σ(k) | Linear growth of a tanh layer (Michalke, 1964) | σ_∞ = 0.18818 / 0.18745 at kδ = 0.4 / 0.5, against 0.1881 / 0.1875 (chapter [01](01-what-kh-sim-computes.md) §6) |
| B — forced 2D turbulence | Energy spectrum slope in the forward cascade (Kraichnan, 1967; McWilliams, 1984) | Slopes −4.12 … −4.16: the finite-Re coherent-vortex regime, steeper than the asymptotic −3, which 256² at Re = 5000 cannot reach |
| C — vortex-pair merger | Two co-rotating patches merge below a critical separation of about 3.3 radii (Melander, Zabusky and McWilliams, 1988) | Merge at 2.5, 3.0, 3.3 radii; no merger at 3.5 and 4.0; identical classification |
| D — forced Taylor–Green | A forced single mode is a fixed point of the RK4 map | Every sample within 1e-6 of the analytic amplitude, at five Reynolds numbers |
| E — *Numerical Recipes* (1986) | Poisson problem of §17.5, convolution of §12.4 | Python, Pascal and Fortran reproduce the printed results (the Poisson grid to its printed precision, the convolution to 2e-16) |

All of it runs through the same request API you use, so a backend cannot pass by a path users do not take.

## §6. Where to next

- **The specification:** [KH-NM](../../../shared/canonical-equations/kh-numerical-methods.md) — every method above with its equations and the per-implementation code mapping.
- **Running it:** chapter [00](00-quick-start.md); the backends' own `README.md` files under [`backends/`](../../../backends/).
- **The same story for a general reader:** the shad tier, [chapter 1 — The billows](../../shad/01-the-billows.md), with a script that measures a growth rate and renders the billows.

### References

- Cox, S. M. and Matthews, P. C. (2002). Exponential time differencing for stiff systems. *J. Comput. Phys.* 176, 430–455.
- Higham, N. J. (2002). *Accuracy and Stability of Numerical Algorithms*, 2nd ed. SIAM.
- Kraichnan, R. H. (1967). Inertial ranges in two-dimensional turbulence. *Phys. Fluids* 10, 1417–1423.
- McWilliams, J. C. (1984). The emergence of isolated coherent vortices in turbulent flow. *J. Fluid Mech.* 146, 21–43.
- Melander, M. V., Zabusky, N. J. and McWilliams, J. C. (1988). Symmetric vortex merger in two dimensions. *J. Fluid Mech.* 195, 303–340.
- Michalke, A. (1964). On the inviscid instability of the hyperbolic-tangent velocity profile. *J. Fluid Mech.* 19, 543–556.
- Orszag, S. A. (1971). On the elimination of aliasing in finite-difference schemes by filtering high-wavenumber components. *J. Atmos. Sci.* 28, 1074.
- Press, W. H., Flannery, B. P., Teukolsky, S. A. and Vetterling, W. T. (1986). *Numerical Recipes in Pascal*. Cambridge University Press.
- Press, W. H., Teukolsky, S. A., Vetterling, W. T. and Flannery, B. P. (2007). *Numerical Recipes*, 3rd ed. Cambridge University Press.

---

*Previous:* [01 — What kh-sim actually computes](01-what-kh-sim-computes.md) · *Back to:* [00 — Quick start](00-quick-start.md)
