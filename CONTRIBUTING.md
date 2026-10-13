# Contributing to kh-sim

Thank you for your interest in contributing. kh-sim is a polyglot reference
implementation of 2D incompressible Navier-Stokes with the Kelvin-Helmholtz
instability benchmark. Contributions in any of the six backend languages
(Python, Rust, Scala, C++, Fortran, Pascal), the shared physics layer, or
the validation infrastructure are all welcome.

Every contribution, regardless of size, is expected to pass the validation
protocol before it is merged. This keeps the cross-backend golden-vector
matrix meaningful as a quality signal.

---

## Ground Rules

- All contributions must pass the **validation protocol** before merging
  (see [Validation Protocol](#validation-protocol) below)
- Code style must match the conventions already present in each backend
  directory -- follow the existing file structure and naming patterns
- One feature or fix per pull request -- keep diffs reviewable
- Commit messages use the imperative mood: `Add spectral Poisson solver`,
  `Fix RK4 time step`, `Port dealiasing kernel to Pascal`
- New contributions must not introduce project-brand identifiers that are
  not part of kh-sim itself; the public-mirror sanitisation rule requires
  that source, documentation, comments, and fixture files stay free of
  unrelated internal or personal-brand tokens (see SECURITY.md for the
  maintainer contact channel)

---

## Validation Protocol

Every backend change must produce diagnostics within tolerance of the
Python/NumPy reference oracle at t=0.1 (100 steps, 64x32 grid, Re=1000).
The canonical anchor point is the wave-1 golden-vector matrix.

| Diagnostic     | Reference   | Tolerance |
|----------------|-------------|-----------|
| kinetic_energy | 0.112810    | +/-5%     |
| enstrophy      | 43.705142   | +/-5%     |
| divergence_rms | 1.20e-14    | < 1e-10   |

For element-wise field-level validation at machine precision, run the
golden-vector matrix tests:

```bash
cd shared/physics && python kh_physics.py
cd tests && python validate_all_backends.py
```

The `shared/physics/golden/` directory contains the canonical golden-vector
JSON files. Each backend is expected to reproduce the field values within
the cross-backend bit-identity tolerance documented in those files.

---

## Local Development

```bash
# Python reference oracle
cd backends/python-numpy && pip install -e ../../[dev] && pytest

# Rust
cd backends/rust && cargo test

# Scala
cd backends/scala && sbt test

# C++
cd backends/cpp && cmake -B build && cmake --build build && ctest --test-dir build

# Fortran
cd backends/fortran && cmake -B build && cmake --build build && ctest --test-dir build

# Pascal
cd backends/pascal && fpc src/kh_server.pas -o build/kh_server
```

---

## Pull Request Checklist

- [ ] Validation protocol passes locally (diagnostics within tolerance)
- [ ] Golden-vector matrix test passes for the affected backend
- [ ] Sanitisation check passes: no project-brand identifiers from outside
  kh-sim have been introduced in source, docs, or comments
- [ ] New code carries SPDX or equivalent license header where the language
  convention supports it (e.g. `// SPDX-License-Identifier: Apache-2.0`)
- [ ] One feature or one fix per PR
- [ ] Commit messages are imperative
- [ ] PR description includes a one-sentence summary of the change

For security-sensitive findings, see SECURITY.md before opening a public issue.
