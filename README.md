# lege-artis / kh-sim

Canonical reference implementation of the 2D incompressible Navier-Stokes equations,
specialised to the Kelvin-Helmholtz shear instability: vorticity-streamfunction form,
pseudo-spectral in space, RK4 in time, 2/3-rule de-aliasing. One specification
(`shared/canonical-equations/`), six independent implementations — the Python oracle,
Fortran, Pascal, C++, Rust and Scala — validated element by element against each other
and against textbook physics.

---

## Status

`v0.1.1` — CI fix on the first public release; see [`RELEASE-NOTES-v0.1.0.md`](RELEASE-NOTES-v0.1.0.md)
for what is validated and how, and [`RELEASE-NOTES-v0.1.1.md`](RELEASE-NOTES-v0.1.1.md) for the fix.

---

## Documentation

Three tiers; each chapter links one tier deeper.

| Tier | For | Start here |
|------|-----|------------|
| Shad | the general reader: the physics as a story, with runnable scripts | [`docs/shad/00-prologue.md`](docs/shad/00-prologue.md) |
| Engineer | running kh-sim, reading its output, knowing what it computes and how | [`docs/engineer/en/00-quick-start.md`](docs/engineer/en/00-quick-start.md) |
| Canonical | the equations, with derivations and the code mapping | [`shared/canonical-equations/`](shared/canonical-equations/) |

---

## Plan

| Version | Scope |
|---------|-------|
| v0.1.0 | Six implementations; golden vectors, contract, physics and canonical textbook testbeds; three documentation tiers; license stack, community files, CI |
| v0.1.x | Shad-tier chapters 2–5; documentation translations; the asymptotic k⁻³ regime of the turbulence testbed at higher Reynolds number |
| v0.3.0 | Lyapunov estimator + invariant-measure tests |
| v0.5.0+ | Stretch: GPU bridge, REST microservice form |

---

## Backends

**Python/NumPy** is the oracle: it defines the canonical vorticity field at each golden
checkpoint and all other backends are validated against it element-wise.
**Rust** (rustfft) and **Scala** (JTransforms) call third-party FFT libraries that share
no code with NumPy or the hand-rolled transforms: their agreement with the oracle is an
independent cross-check of the equations, and they serve the same HTTP contract as typed
services.
**C++** carries a hand-rolled iterative Cooley-Tukey FFT, making every arithmetic step
visible and reviewable — the same design philosophy as the fourier sibling project.
**Fortran** is the canonical hand-rolled reference per the lege-artis language doctrine:
best-referential, best-documented, best-tested among the systems languages here.
**Pascal** is the canonical pedagogical reference: its FFT loop is readable by anyone
who learned numerical methods from a textbook, and it binds directly to the
Numerical Recipes 1986 historical anchor.

---

## License

- Code: Apache License 2.0 — see `LICENSE`
- Documentation: CC-BY-SA-4.0 — see `LICENSE-DOCS`
- Names: `MIM2000` / `Improwave` / `Petr Yamyang` not licensed for derivative use — see `TRADEMARK.md`

---

## Central references

- Orszag (1971) — 2/3-rule spectral dealiasing
- Press, Teukolsky, Vetterling, Flannery — *Numerical Recipes* (2007 3rd ed, general spectral
  methods + Cooley-Tukey FFT reference; 1986 Pascal ed is the direct Pascal-backend anchor)

---

## See also

- Sibling repo: [`lege-artis/fourier`](https://github.com/lege-artis/fourier) — canonical
  DFT/FFT reference implementation; same multi-backend, multi-language discipline this repo
  extends to the full 2D Navier-Stokes + Kelvin-Helmholtz benchmark context.
