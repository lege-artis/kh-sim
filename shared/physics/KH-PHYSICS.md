# Kelvin-Helmholtz Instability — Physics Reference
**Version:** 1.0.0 | **Date:** 2026-03-28
**Kernel:** `kh_physics.py` (Python/NumPy reference)
**Status:** Canonical — all backends must validate against this

---

## 1. Physical Problem

The Kelvin-Helmholtz (KH) instability occurs at the interface between two fluid layers
moving at different velocities. Small perturbations at the shear layer amplify into
characteristic rolling vortex structures (KH billows).

**Canonical setup for kh-sim:**
- 2D incompressible flow
- Doubly periodic domain: x ∈ [0, Lx], y ∈ [0, Ly]
- Shear layer centred at y = Ly/2
- Reynolds number Re = U₀·Lx/ν controls diffusion rate

---

## 2. Governing Equations

### 2D Incompressible Navier-Stokes (vorticity-streamfunction form)

The vorticity-streamfunction formulation eliminates pressure and enforces
incompressibility exactly, making it the canonical form for 2D flow.

**Vorticity transport:**
```
∂ω/∂t + u·∂ω/∂x + v·∂ω/∂y = ν(∂²ω/∂x² + ∂²ω/∂y²)

where: ω = ∂v/∂x - ∂u/∂y  (scalar vorticity, z-component)
       ν = U₀ / Re          (kinematic viscosity)
```

**Streamfunction Poisson equation:**
```
∂²ψ/∂x² + ∂²ψ/∂y² = -ω
```

**Velocity recovery from streamfunction:**
```
u =  ∂ψ/∂y
v = -∂ψ/∂x
```

The continuity equation ∂u/∂x + ∂v/∂y = 0 is satisfied identically by construction.

---

## 3. Initial Conditions

### Velocity field (shear layer + perturbation)
```
u(x, y, 0) = -U₀ · tanh((y - Ly/2) / δ)     (upper stream in -x; sign corrected 2026-10-12)

v(x, y, 0) = A · sin(2π·k·x / Lx)

where: δ = 0.05·Ly  (shear layer thickness)
       A = perturbation_amplitude  (e.g. 0.01)
       k = perturbation_mode       (e.g. 2)
```

### Initial vorticity (analytically from IC velocity field)
```
ω(x, y, 0) = ∂v/∂x - ∂u/∂y
            = A·(2πk/Lx)·cos(2πkx/Lx)
            + U₀/δ · sech²((y - Ly/2)/δ)

Note: tanh shear gives sech² vorticity sheet;
      v-perturbation gives cosine vorticity perturbation.
```

**Default parameters:**
| Parameter              | Symbol | Default  |
|------------------------|--------|----------|
| Domain length          | Lx     | 1.0      |
| Domain height          | Ly     | 0.5      |
| Shear velocity         | U₀     | 1.0      |
| Reynolds number        | Re     | 1000     |
| Perturbation amplitude | A      | 0.01     |
| Perturbation mode      | k      | 2        |
| Shear layer thickness  | δ      | 0.05·Ly  |

---

## 4. Numerical Method

### Spatial discretisation: Pseudo-spectral (Fourier)
- Spatial derivatives computed exactly via 2D FFT
- Doubly periodic BC enforced by Fourier basis
- Grid: Nx × Ny uniform points (default 128 × 64)

### Poisson solver (spectral)
```
In Fourier space: ψ̂(kx, ky) = ω̂(kx, ky) / (kx² + ky²)

Zero-mode: ψ̂(0,0) = 0  (fixes mean streamfunction = 0)
```

### Time integration: 4th-order Runge-Kutta (RK4)
```
k₁ = f(ωⁿ)
k₂ = f(ωⁿ + ½·dt·k₁)
k₃ = f(ωⁿ + ½·dt·k₂)
k₄ = f(ωⁿ + dt·k₃)
ωⁿ⁺¹ = ωⁿ + (dt/6)(k₁ + 2k₂ + 2k₃ + k₄)

where f(ω) = -u∂ω/∂x - v∂ω/∂y + ν∇²ω
```

### Stability criterion (CFL)
```
dt ≤ min(dx/|U|_max, dy/|V|_max)   (advection)
dt ≤ min(dx², dy²) / (2ν)          (diffusion)
```
Default dt = 0.001 is conservative for Re=1000 on a 128×64 grid.

---

## 5. Output Fields (REST API — SimulationResult)

| Field       | Shape    | Description                              |
|-------------|----------|------------------------------------------|
| u_velocity  | Nx × Ny  | x-velocity field (row-major flattened)   |
| v_velocity  | Nx × Ny  | y-velocity field (row-major flattened)   |
| vorticity   | Nx × Ny  | ω = ∂v/∂x - ∂u/∂y (row-major flattened) |
| pressure    | Nx × Ny  | ψ (streamfunction, pressure surrogate)   |

**Note on pressure field:** The REST API specifies `pressure` in the field snapshot.
In the vorticity-streamfunction formulation, pressure is not the primary variable.
We expose ψ (streamfunction) as the `pressure` field — it provides equivalent
visualisation information and encodes the flow topology. For the true pressure field,
backends may optionally solve the pressure Poisson equation:
∇²p = -2(∂u/∂x · ∂v/∂y - ∂u/∂y · ∂v/∂x)

### Diagnostics
| Diagnostic    | Formula                        |
|---------------|--------------------------------|
| kinetic_energy| ½ · mean(u² + v²)             |
| enstrophy     | ½ · mean(ω²)                  |
| max_vorticity | max(|ω|)                      |
| divergence_rms| rms(∂u/∂x + ∂v/∂y) ≈ 0 (spectral) |

---

## 6. Backend Adaptation Notes

### All backends must:
1. Accept `SimulationRequest` JSON (see `shared/api/openapi.yaml`)
2. Implement the same IC formula (tanh shear + sinusoidal perturbation)
3. Return fields in **row-major order** (x varies first, i.e. field[ix*ny + iy])
4. Validate against the golden vector `shared/golden-vectors/v1-canonical-reference.json` (tolerance: its A2 band, §7)

### Language-specific FFT libraries
| Backend | FFT library              |
|---------|--------------------------|
| Python  | numpy.fft (reference)    |
| Rust    | rustfft                  |
| Scala   | Breeze or JTransforms    |
| C++     | FFTW3                    |
| Fortran | fftpack5 or FFTW3        |
| Pascal  | Custom DFT (small grids) or FFTW3 binding |

### Simplification for early backends
Backends may use finite-difference + SOR Poisson solver instead of FFT
for initial implementation. Validate diagnostics within 5% of reference output.
Once correctness is confirmed, optimise to spectral method.

---

## 7. Reference Validation

The reference is the oracle golden vector **V1**, `shared/golden-vectors/v1-canonical-reference.json`
(class A2, `_audit/KH-SIM-MODELLING-RULES-v0.1.md` §4): the request (`parameters`), the expected
diagnostics and field summaries, and the tolerance (`tolerances`) in one file, generated by
`kh_physics.py` at a recorded scheme version (post-5c, 2/3-rule dealiasing). It replaced
`shared/physics/kh_reference_output.json`, retired 2026-10-08: a pre-5c copy that differed from the
oracle by 8.5 × 10⁻¹⁰ and was graded at ±5 %, which hid the difference.

**Reference diagnostics (V1)** — t=0.1 (100 steps, dt=0.001), 64×32 grid, Re=1000:
```json
{
  "t_final": 0.1,
  "diagnostics": {
    "kinetic_energy":  0.11280980368298389,
    "enstrophy":       43.70514370056882,
    "max_vorticity":   31.908804110189003,
    "divergence_rms":  1.148e-14,
    "cfl_peak":        0.049731045202224225
  }
}
```
Each backend's diagnostics — natively (`validation_test.*`, `kh_validate.pas`) and over HTTP
(`tests/integration/test_simulate.py`, `test_golden_vectors_cross_backend.py`) — must agree with V1
within its A2 band: |got − ref| ≤ max(scalar_abs·f, scalar_rel·f·|ref|), with f the backend's
`backend_overrides` factor (Scala 2, Pascal 10, others 1).
divergence_rms < 1e-10 is the incompressibility acceptance criterion.
