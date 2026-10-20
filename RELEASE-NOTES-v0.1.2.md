# lege-artis/kh-sim — v0.1.2

**Tag:** `v0.1.2`
**Type:** audits and comments. The solver's behaviour is unchanged; the source changes are comment lines only.

## Citation audit (doctrine C3)

The routines of the hand-rolled backends that implement a textbook method carry a citation header: reference (author, year, section, equation), sign or convention, and complexity. `_audit/textbook-citations-allowlist.yaml` lists them; `scripts/audit_citations.py` (wrapper `scripts/audit-citations.sh`) checks the sources against it. The workflow `citations-c3.yml` runs the audit on every push that touches the backends.

**What the audit checks**

- Each cited routine is defined in its file, and the comment block above it names the reference, section, equation (and page where given), the convention and the complexity.
- Every public Fortran and Pascal routine is in the allowlist. It is either cited, pending (anchor still to be chosen), or exempt with a reason. A new public routine without an entry fails the audit.
- `--self-test` proves that the audit fails on a removed reference line, a removed complexity line, a stale entry and an unlisted public routine.

**State in this release**

- 7 cited routines:
  - FFT — Fortran, Pascal, C++ (Press et al., Numerical Recipes, §12.2, eq. 12.2.2);
  - classical RK4 — Fortran (Numerical Recipes 3rd ed., §17.1, eq. 17.1.3) and Pascal (Numerical Recipes in Pascal, §15.1, eq. 15.1.3);
  - 2/3-rule de-aliasing — Fortran, Pascal (Orszag 1971).
- 8 pending:
  - spectral Poisson solve, spectral derivative and advection term in Fortran and Pascal;
  - the ETDRK4 alternate scheme.
- 30 exempt: drivers, set-up, I/O, forcing construction, helpers.
- Rust (rustfft) and Scala (JTransforms) are library-backed and exempt from the doctrine criteria C1–C4.

## Complexity audit (doctrine C4)

`backends/pascal/tests/test_complexity_c4.pas` and `backends/fortran/tests/test_complexity_c4.f90` time the hand-rolled FFTs at N = 2^6 … 2^16. Each size is run five times and the median taken. A least-squares fit gives the slope of ln t over ln(N log₂ N); the gate is [0.9, 1.1].

| Routine | Slope | Result |
|---|---|---|
| Pascal `FFT1D` | 0.98 | PASS |
| Fortran `fft1d_inplace` (1D) | 0.94 | PASS |

These are the first measurements, on a Linux x86-64 container. For the Fortran 2D routines the slope is printed for information only; no complexity claim is made on them.

The workflow `complexity-c4.yml` runs both audits nightly and on demand. Run time depends on the runner, so one failed night is re-run before it is read as a regression.

## Source changes

Comment lines only:
- citation headers on the routines above (`kh_fft.f90`, `kh_physics.f90`, `kh_grid.f90`, `kh_physics.pas`);
- the complexity line on the C++ `fft1d` (`physics.cpp`).
