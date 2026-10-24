# kh-sim golden vectors

Reference outputs from the Python oracle (`shared/physics/kh_physics.py`) for
cross-backend validation. All backends (Rust, Scala, C++, Fortran, Pascal) must
reproduce these diagnostics and field summaries within the specified tolerances.

Per phase 6a design: development record, phase-6a golden-vector design v0.1.

---

## Vector matrix

9 vectors spanning the active parameter space (Re sweep, grid refinement,
long-time evolution, alternate modes, dealiasing edge, zero-perturbation base flow).

| ID | File | Grid | T_end | Re | amp | mode | Purpose |
|----|------|------|-------|----|-----|------|---------|
| V1 | [v1-canonical-reference.json](v1-canonical-reference.json) | 64×32 | 0.1 | 1000 | 0.01 | 2 | **Canonical reference.** Phase-5 baseline. |
| V2 | [v2-low-re-diffusive.json](v2-low-re-diffusive.json) | 64×32 | 0.1 | 100 | 0.01 | 2 | **Low-Re diffusive.** ν=0.01; strong viscous damping; KE decreases monotonically. |
| V3 | [v3-high-re-nonlinear.json](v3-high-re-nonlinear.json) | 64×32 | 0.1 | 10000 | 0.01 | 2 | **High-Re nearly-inviscid.** ν=0.0001; nonlinear cascade; dealiasing stress test. |
| V4 | [v4-double-resolution.json](v4-double-resolution.json) | 128×64 | 0.1 | 1000 | 0.01 | 2 | **Double-resolution refinement.** 4× spatial points; asymptotic convergence check. |
| V5 | [v5-half-resolution-fast.json](v5-half-resolution-fast.json) | 32×16 | 0.1 | 1000 | 0.01 | 2 | **Half-resolution fast.** 1/4 spatial points; suitable for fast CI unit-test runs. |
| V6 | [v6-long-time-evolution.json](v6-long-time-evolution.json) | 64×32 | 0.5 | 1000 | 0.01 | 2 | **Long-time evolution.** T_end=0.5; nonlinear saturation phase; vortex rollup. |
| V7 | [v7-alt-mode-k1.json](v7-alt-mode-k1.json) | 64×32 | 0.1 | 1000 | 0.01 | 1 | **Alternate mode k=1.** Different unstable wavenumber; checks scheme symmetry. |
| V8 | [v8-dealiasing-edge-k4.json](v8-dealiasing-edge-k4.json) | 64×32 | 0.1 | 1000 | 0.01 | 4 | **Dealiasing-edge mode k=4.** Energy near 2/3 cutoff; dealiasing precision test. |
| V9 | [v9-zero-perturbation.json](v9-zero-perturbation.json) | 64×32 | 0.1 | 1000 | 0 | 2 | **Zero-perturbation base-flow.** amp=0; quasi-steady viscous decay; no KH growth. |

### Oracle diagnostics (generation run)

| ID | KE | Enstrophy | max_vort | div_rms | cfl_peak |
|----|-----|-----------|---------|---------|---------|
| V1 | 0.112810 | 43.7051 | 31.909 | 1.15e-14 | 0.0497 |
| V2 | 0.073757 | 22.5734 | 16.120 | 1.25e-14 | 0.0497 |
| V3 | 0.119914 | 52.0398 | 38.948 | 1.11e-14 | 0.0502 |
| V4 | 0.112807 | 43.7029 | 31.901 | 2.41e-14 | 0.0497 |
| V5 | 0.116437 | 46.2215 | 33.057 | 4.44e-15 | 0.0124 |
| V6 | 0.091228 | 29.2271 | 21.466 | 9.83e-15 | 0.0497 |
| V7 | 0.112810 | 43.7024 | 31.827 | 1.23e-14 | 0.0497 |
| V8 | 0.112808 | 43.7151 | 32.078 | 1.19e-14 | 0.0497 |
| V9 | 0.112786 | 43.7015 | 31.761 | 0.00e+00 | 0.0491 |

All `div_rms` < 1.5×10⁻¹⁴ — well below the 10⁻¹⁰ incompressibility acceptance criterion.

---

## JSON schema

Each vector file is validated against [`schema.json`](schema.json)
(JSON Schema Draft 2020-12). The schema enforces:

- `metadata` — vector identity, provenance, oracle commit, environment
- `parameters` — exact simulation parameters (grid, domain, dt, nsteps, Re, amp, mode)
- `diagnostics` — 6 scalar values: KE, enstrophy, max_vort, div_rms, cfl_peak, t_final
- `field_summaries` — 4 fields × 6 statistics each: u, v, ω (vorticity), ψ (pressure/streamfunction)
- `tolerances` — default Higham §4.2 hybrid bounds + per-backend override factors

Validate any vector manually:
```
python -c "
import json, jsonschema
schema = json.load(open('shared/golden-vectors/schema.json'))
doc    = json.load(open('shared/golden-vectors/v1-canonical-reference.json'))
jsonschema.validate(instance=doc, schema=schema)
print('PASS')
"
```

---

## Regeneration

Vectors are **immutable ground truth**. Regenerate only when:
- `shared/physics/kh_physics.py` changes numerically (not pure refactoring)
- Canonical equations spec (KH-GE-*, KH-NM-*) is amended
- 2/3-rule dealiasing or CFL guard semantics change
- A new vector is added to the matrix

To regenerate all 9 vectors:
```
python tools/generate_kh_golden_vectors.py --out shared/golden-vectors/ --force
```

To regenerate a subset:
```
python tools/generate_kh_golden_vectors.py --only V1,V4 --force
```

The generator script requires NumPy and (optionally) jsonschema:
```
pip install numpy scipy jsonschema
```

Exit codes: `0` = success, `1` = oracle error, `2` = schema validation failed, `3` = output not writable.

---

## Comparison rule — Higham §4.2 hybrid

For each scalar diagnostic:
```
|backend_value - oracle_value| ≤ max(scalar_abs, scalar_rel × |oracle_value|)
```

For field summary statistics (min/max/mean per field):
```
|backend - oracle| ≤ max(field_summary_abs, field_summary_rel × |oracle|)
```

For first_10 / last_10 element-wise comparisons, an additional
sqrt(N) × machine_eps × |oracle| term accounts for accumulated rounding
(Higham §4.2 floating-point accumulation error growth), where N = nx × ny.

---

## Per-backend tolerance overrides

Default tolerances apply to Python ref (self-comparison), Rust, C++, and Fortran.
Some backends have wider rounding profiles:

| Backend | Class | Override factor | Rationale |
|---------|-------|-----------------|-----------|
| `kh-rust` | — | 1× | f64 throughout; tight bounds achievable |
| `kh-cpp` | — | 1× | f64 throughout |
| `kh-fortran` | — | 1× | gfortran f64 canonical |
| `kh-scala` | scalar + field | **2×** | JTransforms FFT butterfly order differs; accumulated rounding ~2× wider |
| `kh-pascal` | scalar + field | **10×** | FPC `double` + manual FFT loops; wider rounding profile |

Override factors multiply the default abs and rel thresholds stored in each
vector's `tolerances.backend_overrides` block. Example: Pascal scalar abs =
1×10⁻¹² × 10 = 1×10⁻¹¹.

If a backend repeatedly fails even with its override factor, that signals a
**port bug**, not tolerance miscalibration. Open a phase 6c follow-up issue.

---

## Open questions (phase 6a §10)

- **OQ-GOLD-1:** Pascal `double` precision on Apple Silicon — may tighten or widen
  the 10× override; measure during phase 6c Pascal-port test.
- **OQ-GOLD-2:** V6 long-time field-element tolerance — accumulated rounding grows
  sqrt(500)/sqrt(100) ≈ 2.24× relative to V1; relax field_element bounds if backends
  show systematic drift during phase 6c.
- **OQ-GOLD-3:** V9 v-velocity near-zero — field_rel × 0 = 0 → pure abs comparison;
  handled correctly by hybrid rule (max(1e-10, 1e-7 × 0) = 1e-10 abs threshold).
- **OQ-GOLD-4:** Cross-backend HTTP tests require all 5 backends running; intended for
  local bring-up on the reference workstation initially; Docker compose deferred to phase 9.
