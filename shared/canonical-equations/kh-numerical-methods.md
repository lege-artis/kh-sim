# kh-sim canonical equation — numerical methods (KH-NM)

**Status:** canonical (TeX-grade)
**Version:** 0.1 — first cut, 2026-05-12
**Companion to:** `kh-numerical-methods.tex`
**Cited by:** all five backends; specifically `kh-sim/backends/fortran/tests/test_num_*.f90`
which traces to this document for the TC-NUM-KH-001..008 specifications.

This document specifies the **numerical methods** the kh-sim solver uses to discretise
and integrate the governing equations from `kh-governing-equations.md`. Three pieces:
spatial discretisation (pseudo-spectral), time integration (RK4 canonical, ETDRK4
alternate), and aliasing control (2/3-rule).

---

## KH-NM-1: pseudo-spectral spatial discretisation

The doubly periodic domain Ω = [0, Lx] × [0, Ly] is discretised on a uniform grid
of Nx × Ny points. Spatial derivatives are computed exactly in Fourier space via
the 2D discrete Fourier transform.

For a real-valued field f on the grid, the 2D real-input FFT (rfft2) maps
f(x, y) → f̂(kx, ky) where the wavenumber arrays are:

$$
k_x = \frac{2\pi}{L_x}\,\{0, 1, 2, \ldots, N_x/2-1, -N_x/2, \ldots, -1\}
$$
$$
k_y = \frac{2\pi}{L_y}\,\{0, 1, 2, \ldots, N_y/2\}
$$

(rfft2 exploits the conjugate symmetry of the real input by storing only the
non-negative-ky half; the kx axis carries the full FFT range.)

Spatial derivatives become multiplications in Fourier space:

$$
\widehat{\partial_x f}(k_x, k_y) = i k_x\,\hat f(k_x, k_y),
\qquad
\widehat{\partial_y f}(k_x, k_y) = i k_y\,\hat f(k_x, k_y).
$$

The Laplacian becomes:

$$
\widehat{\nabla^2 f}(k_x, k_y) = -(k_x^2 + k_y^2)\,\hat f(k_x, k_y).
\quad\text{(KH-NM-1)}
$$

## KH-NM-2: spectral Poisson solver

The streamfunction Poisson equation KH-GE-2 ψ_xx + ψ_yy = -ω is diagonal in
Fourier space (per KH-NM-1):

$$
-(k_x^2 + k_y^2)\,\hat\psi(k_x, k_y) = -\hat\omega(k_x, k_y),
$$

which gives the immediate solution:

$$
\hat\psi(k_x, k_y) = \frac{\hat\omega(k_x, k_y)}{k_x^2 + k_y^2}
\qquad\text{for}\quad (k_x, k_y) \neq (0, 0),
$$
$$
\hat\psi(0, 0) = 0
\qquad\text{(fixes mean streamfunction to zero).}
\quad\text{(KH-NM-2)}
$$

The zero-mode is set to zero because the Poisson equation determines ψ only
up to an additive constant; the choice ψ̂(0,0) = 0 is conventional. This
corresponds to TC-NUM-KH-002 (spectral Poisson solver correctness).

## KH-NM-3: time integration

### KH-NM-3a: canonical Runge-Kutta 4 (RK4)

The vorticity-transport equation KH-GE-1 is written abstractly as

$$
\frac{\partial \omega}{\partial t} = f(\omega) = -u\,\frac{\partial \omega}{\partial x} - v\,\frac{\partial \omega}{\partial y} + \nu\,\nabla^2 \omega
$$

where u, v are recovered from ω via KH-NM-2 + KH-GE-3. The classical 4th-order
Runge-Kutta scheme advances ω from time tⁿ to tⁿ⁺¹ = tⁿ + Δt via:

$$
\begin{aligned}
k_1 &= f(\omega^n),\\
k_2 &= f(\omega^n + \tfrac{1}{2}\Delta t\,k_1),\\
k_3 &= f(\omega^n + \tfrac{1}{2}\Delta t\,k_2),\\
k_4 &= f(\omega^n + \Delta t\,k_3),\\
\omega^{n+1} &= \omega^n + \tfrac{\Delta t}{6}\,(k_1 + 2k_2 + 2k_3 + k_4).
\end{aligned}
\quad\text{(KH-NM-3a)}
$$

This is the canonical scheme used by `kh_physics.py` (Python reference) and
`backends/rust/src/physics/solver.rs`. Per Hairer-Wanner standard treatment.

### KH-NM-3b: alternate ETDRK4 (Cox & Matthews 2002)

When the viscous diffusion term ν∇²ω is stiff (high Re or under-resolved grids),
the **exponential time-differencing Runge-Kutta 4** scheme handles the linear
diffusion term **exactly** by matrix-exponential integration, with RK4-like
stages for the non-linear advection term. For the linear-plus-nonlinear split

$$
\frac{\partial \omega}{\partial t} = L\,\omega + N(\omega)
$$

where L = ν∇² (linear, diagonal in Fourier space: $L = -\nu(k_x^2 + k_y^2)$) and
N(ω) = -u∂ω/∂x - v∂ω/∂y (non-linear advection), the ETDRK4 update is:

$$
\omega^{n+1} = e^{L\Delta t}\omega^n + \Delta t\,\big[\,\phi_1\,N_1 + 2\phi_2(N_2 + N_3) + \phi_3\,N_4\,\big]
\quad\text{(KH-NM-3b)}
$$

where φ_i are contour-integral-evaluated φ-functions (Kassam & Trefethen 2005)
and N_i are non-linear-term evaluations at intermediate stages. The non-linear
stages use RK4-like sub-step structure.

This is the scheme used by `backends/fortran/src/kh_etdrk4.f90`. Per
AUDIT-kh-sim §6.1, the choice between KH-NM-3a and KH-NM-3b is **flagged for
disambiguation in phase 5**. Until then, both schemes co-exist with
documented divergence behaviour.

## KH-NM-4: 2/3-rule dealiasing (Orszag 1971)

The non-linear advection term f(ω) involves products of fields (e.g. u·∂ω/∂x),
which in Fourier space corresponds to convolution. On a finite grid this
convolution produces aliasing of high-wavenumber content.

The 2/3-rule is the standard mitigation: zero out the highest one-third of
Fourier modes after each non-linear-term evaluation. For wavenumber arrays of
size Nx × (Ny/2+1), this means setting:

$$
\hat g(k_x, k_y) := 0 \quad\text{if}\quad |k_x| > \frac{N_x}{3} \;\;\text{or}\;\; |k_y| > \frac{N_y}{3}.
\quad\text{(KH-NM-4)}
$$

After dealiasing, the surviving Fourier modes contain only the aliasing-free
representation of the product. This is TC-NUM-KH-005.

**Audit note (AUDIT-kh-sim §6.3):** the Python reference kernel `kh_physics.py`
may or may not implement KH-NM-4. If it does not, there is a numerical drift
between Python ref + Rust backend (no dealiasing) and Fortran backend (has
dealiasing). To be confirmed and aligned in phase 5.

## KH-NM-5: CFL stability constraint

For the explicit RK4 time integrator (KH-NM-3a) and a periodic Fourier-spectral
spatial discretisation, the maximum stable timestep is bounded by:

$$
\Delta t \leq \min\!\left( \frac{\Delta x}{|u|_{\max}}, \frac{\Delta y}{|v|_{\max}} \right)
\qquad\text{(advective CFL)},
$$
$$
\Delta t \leq \min\!\left( \frac{(\Delta x)^2}{2\nu}, \frac{(\Delta y)^2}{2\nu} \right)
\qquad\text{(diffusive CFL)}.
\quad\text{(KH-NM-5)}
$$

The kh-sim default Δt = 0.001 is conservative for the default Re = 1000 on a
128 × 64 grid (Δx ≈ 7.8e-3, |u|_max ≈ 1 → advective bound 7.8e-3; ν = 1e-3 →
diffusive bound ≈ 3e-2). This is TC-NUM-KH-006.

The ETDRK4 scheme (KH-NM-3b) relaxes the diffusive CFL by treating the linear
term exactly; the advective bound remains.

---

## Test-case traceability

| Spec | KH-NM section | Fortran test |
|---|---|---|
| TC-NUM-KH-001 | KH-NM-1 (2D FFT round-trip) | `test_num_001_fft_roundtrip.f90` |
| TC-NUM-KH-002 | KH-NM-2 (spectral Poisson) | `test_num_002_poisson.f90` |
| TC-NUM-KH-003 | KH-NM-3b (ETDRK4 on linear scalar ODE) | `test_num_003_etdrk4_linear.f90` |
| TC-NUM-KH-004 | KH-NM-3b (ETDRK4 temporal convergence order) | `test_num_004_convergence.f90` |
| TC-NUM-KH-005 | KH-NM-4 (2/3-rule dealiasing) | `test_num_005_dealias.f90` |
| TC-NUM-KH-006 | KH-NM-5 (CFL stability across Re) | `test_num_006_cfl.f90` |
| TC-NUM-KH-007 | KH-NM-3a + KH-NM-2 (energy conservation, inviscid limit) | `test_num_007_energy.f90` |
| TC-NUM-KH-008 | KH-NM-1..5 composed (reference run validation) | `test_num_008_reference.f90` |

---

## Code mapping per backend

| Equation | Python ref | Rust | Scala | C++ | Fortran | Pascal |
|---|---|---|---|---|---|---|
| KH-NM-1 (FFT + wavenumbers) | `shared/physics/kh_physics.py:69-77` | `backends/rust/src/physics/fft2d.rs:107-145` | `backends/scala/src/main/scala/khsim/physics/Fft2D.scala:50-80` | `backends/cpp/src/physics.cpp:103-113` (1D helper) + `backends/cpp/src/physics.cpp:274-281` (meshgrid in `simulate`) | `backends/fortran/src/kh_fft.f90:114-166` (wavenumber array building via `kh_grid.f90`) | `backends/pascal/src/kh_physics.pas:154-168` |
| KH-NM-2 (spectral Poisson) | `shared/physics/kh_physics.py:109-121` | `backends/rust/src/physics/solver.rs:52-77` | `backends/scala/src/main/scala/khsim/physics/Solver.scala:39-64` | `backends/cpp/src/physics.cpp:138-155` | `backends/fortran/src/kh_poisson.f90:42-59` | `backends/pascal/src/kh_physics.pas:190-205` |
| KH-NM-3a (RK4) | `shared/physics/kh_physics.py:153-169` | `backends/rust/src/physics/solver.rs:159-187` | `backends/scala/src/main/scala/khsim/physics/Solver.scala:148-172` | `backends/cpp/src/physics.cpp:208-226` | `backends/fortran/src/kh_physics.f90:129-175` (`rk4_step`, HTTP-shim path) | `backends/pascal/src/kh_physics.pas:259-295` |
| KH-NM-3b (ETDRK4) | (not used) | (not used) | (not used) | (not used) | `backends/fortran/src/kh_etdrk4.f90:128-219` (alternate-scheme solver via `kh_solver.f90`) | (not used) |
| KH-NM-4 (2/3-rule dealiasing) | `shared/physics/kh_physics.py:137-145` (`dealias_two_thirds`); applied in `vorticity_rhs:167` | `backends/rust/src/physics/solver.rs:109-119` (`dealias_two_thirds`); applied at line 175 | `backends/scala/src/main/scala/khsim/physics/Solver.scala:147-164` (`dealiasTwoThirds`); applied at line 135 | `backends/cpp/src/physics.cpp:177-188` (`dealias_two_thirds`); applied at line 223 | `backends/fortran/src/kh_nonlinear.f90:111-122` (`kh_dealias`); applied in `kh_physics.f90:120` (HTTP-shim RK4 path) | `backends/pascal/src/kh_physics.pas:174-189` (`DealiasField`); applied at line 290 |
| KH-NM-5 (CFL guard) | `shared/physics/kh_physics.py:147-151` (`cfl_dt_bound`); guard at lines 220-227 | `backends/rust/src/physics/solver.rs:122-126` (`cfl_dt_bound`); guard at lines 284-294 | `backends/scala/src/main/scala/khsim/physics/Solver.scala:167-170` (`cflDtBound`); guard at lines 284-292 | `backends/cpp/src/physics.cpp:191-195` (`cfl_dt_bound`); guard at lines 327-340 | `backends/fortran/src/kh_physics.f90:247-259` (inline CFL guard in `kh_simulate_c`) | `backends/pascal/src/kh_physics.pas:193-199` (`CflDtBound`); guard at lines 360-371 |

---

## References

See `shared/reference-bibliography/refs.bib`:
- Boyd (2001) chs. 1–2 — Fourier spectral methods
- Trefethen (2000) chs. 3–4 — spectral methods with 2D Poisson via FFT
- Hairer, Nørsett, Wanner (1993) — classical RK4
- Cox & Matthews (2002) — ETDRK4 seminal paper
- Kassam & Trefethen (2005) — ETDRK4 contour-integral coefficient stabilisation
- Orszag (1971) — 2/3-rule dealiasing
