# kh-sim canonical equation — governing equations (KH-GE)

**Status:** canonical (TeX-grade)
**Version:** 0.1.1 — 2026-10-12 (KH-GE-4: sign of u corrected); 0.1 first cut 2026-05-12
**Companion to:** `kh-governing-equations.tex`
**Cited by:** all five backends (Python reference, Rust, Scala, C++, Fortran, Pascal)

This is the canonical statement of what kh-sim simulates. Every backend's code is
a line-by-line implementation of the equations stated here (Equations KH-GE-1 through
KH-GE-5). If a backend computes anything different, either the backend is wrong or
this document is incomplete — never both.

---

## KH-GE-1: vorticity-transport equation

For two-dimensional incompressible Navier-Stokes flow in vorticity-streamfunction
form, the scalar vorticity ω evolves under advection by the velocity field
(u, v) and viscous diffusion at kinematic viscosity ν:

$$
\frac{\partial \omega}{\partial t}
+ u\,\frac{\partial \omega}{\partial x}
+ v\,\frac{\partial \omega}{\partial y}
= \nu\left( \frac{\partial^2 \omega}{\partial x^2} + \frac{\partial^2 \omega}{\partial y^2} \right)
\quad\text{(KH-GE-1)}
$$

where the scalar vorticity is

$$
\omega = \frac{\partial v}{\partial x} - \frac{\partial u}{\partial y}
\quad\text{(z-component, 2D)}
$$

and the kinematic viscosity is set by the Reynolds number Re and the characteristic
shear velocity U₀:

$$
\nu = \frac{U_0}{\mathrm{Re}}.
$$

## KH-GE-2: streamfunction Poisson equation

The streamfunction ψ recovers the velocity field from the vorticity. It satisfies
the Poisson equation:

$$
\frac{\partial^2 \psi}{\partial x^2} + \frac{\partial^2 \psi}{\partial y^2}
= -\omega
\quad\text{(KH-GE-2)}
$$

## KH-GE-3: velocity recovery

$$
u = \frac{\partial \psi}{\partial y},
\qquad
v = -\frac{\partial \psi}{\partial x}.
\quad\text{(KH-GE-3)}
$$

Continuity ∂u/∂x + ∂v/∂y = 0 is satisfied identically by construction.

## KH-GE-4: initial conditions

Velocity field — tanh shear layer with sinusoidal perturbation:

$$
u(x, y, 0) = -\,U_0\,\tanh\!\left( \frac{y - L_y/2}{\delta} \right),
$$
$$
v(x, y, 0) = A\,\sin\!\left( \frac{2\pi k x}{L_x} \right).
$$

Initial vorticity (analytically from the velocity IC):

$$
\omega(x, y, 0) = A\,\frac{2\pi k}{L_x}\,\cos\!\left( \frac{2\pi k x}{L_x} \right)
+ \frac{U_0}{\delta}\,\operatorname{sech}^2\!\left( \frac{y - L_y/2}{\delta} \right).
\quad\text{(KH-GE-4)}
$$

The tanh shear gives the sech² vorticity sheet; the v-perturbation gives the
cosine vorticity perturbation. For U₀ > 0 the upper stream moves in −x and the
lower one in +x.

*Corrected 2026-10-12 (v0.1.1): the sign of u. The vorticity above is the one every
backend has always computed and every golden vector holds; with
ω = ∂v/∂x − ∂u/∂y it belongs to u = −U₀ tanh, not to +U₀ tanh as v0.1 wrote. The
two are mirror images; the linear stability (KH-LS) is invariant under the
reflection.*

**Default parameters** (matching kh-sim's `kh_physics.py`):

| Symbol | Meaning | Default |
|---|---|---|
| Lx | domain length | 1.0 |
| Ly | domain height | 0.5 |
| U₀ | shear velocity | 1.0 |
| Re | Reynolds number | 1000 |
| A | perturbation amplitude | 0.01 |
| k | perturbation mode | 2 |
| δ | shear-layer thickness | 0.05·Ly |

## KH-GE-5: boundary conditions

The domain is **doubly periodic**:

$$
\omega(0, y, t) = \omega(L_x, y, t),
\qquad
\omega(x, 0, t) = \omega(x, L_y, t),
\quad\text{(KH-GE-5)}
$$

and equivalently for ψ. Doubly periodic BCs are exact for the Fourier basis
(pseudo-spectral discretisation, KH-NM document).

---

## Equation-to-code mapping per backend

Every backend's solver header MUST reference this document and provide a
line-by-line mapping from each KH-GE-N equation to specific code lines.

| Code path | KH-GE-1 | KH-GE-2 | KH-GE-3 | KH-GE-4 | KH-GE-5 |
|---|---|---|---|---|---|
| `shared/physics/kh_physics.py` | `vorticity_rhs()` lines 136-150 | `solve_poisson()` lines 109-121 | `velocity_from_psi_hat()` lines 124-133 | `initial_conditions()` lines 80-106 | implicit (periodic FFT basis) |
| `backends/rust/src/physics/solver.rs` | `vorticity_rhs()` lines 109-151 | `solve_poisson()` lines 52-77 | `velocity_from_psi_hat()` lines 80-104 | `initial_conditions()` lines 16-49 | implicit |
| `backends/cpp/src/physics.cpp` | `vorticity_rhs()` lines 174-198 | `solve_poisson()` lines 138-155 | `velocity_from_psi()` lines 159-171 | `initial_conditions()` lines 116-135 | implicit |
| `backends/fortran/src/kh_physics.f90` | `vorticity_rhs()` lines 203-233 (RK4 shim); nonlinear term in `kh_nonlinear.f90:53-97` (ETDRK4 path) | `kh_poisson.f90:42-59` | `kh_velocity.f90:45-73` | `initial_conditions()` lines 139-157 | implicit |
| `backends/pascal/src/kh_physics.pas` | `VortRHS()` lines 230-256 | `SolvePoisson()` lines 190-205 | `VelFromPsi()` lines 210-227 | `InitConds()` lines 171-187 | implicit |
| `backends/scala/.../physics/Solver.scala` | `vorticityRhs()` lines 92-134 | `solvePoisson()` lines 39-64 | `velocityFromPsiHat()` lines 67-89 | `initialConditions()` lines 17-36 | implicit |

(TBD entries land in phase 5 of the discipline lift, after the canonical-scheme
disambiguation in AUDIT-kh-sim §6.1.)

---

## References

See `shared/reference-bibliography/refs.bib`:
- Drazin & Reid (2004) ch. 4 — KH instability theory
- Schmid & Henningson (2001) ch. 3 — free shear layers, modern numerical-friendly treatment
- Kelvin (1871), Helmholtz (1868) — historical
