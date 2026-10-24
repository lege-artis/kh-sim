# lege-artis/kh-sim — v0.1.3

**Tag:** `v0.1.3`
**Type:** comments, audits and housekeeping. The solver's behaviour is unchanged. With one exception, the source changes are comment and docstring lines only: the banner string of `tests/integration/run-integration.ps1`.

## Citation audit (doctrine C3): every routine cited

v0.1.2 left eight routines pending. They now carry citation headers. The allowlist moves to schema 0.3.0 and lists:
- 18 hand-rolled routines, cited and checked;
- 2 library-backed FFTs;
- 30 exempt routines;
- no pending routines.

| Method | Routines | Reference |
|---|---|---|
| FFT | Fortran `fft1d_inplace`, Pascal `FFT1D`, C++ `fft1d` | Press et al., *Numerical Recipes*, §12.2, eq. 12.2.2 |
| Classical RK4 | Fortran `rk4_step`, Pascal `RK4Step` | *Numerical Recipes* 3rd ed. §17.1, eq. 17.1.3; *in Pascal* §15.1, eq. 15.1.3 |
| 2/3-rule de-aliasing | Fortran `kh_grid_dealias_mask`, Pascal `DealiasField` | Orszag (1971), J. Atmos. Sci. 28(6), p. 1074 |
| Spectral Poisson solve | Fortran `kh_poisson_solve`, `solve_poisson`; Pascal `SolvePoisson` | Trefethen (2000), *Spectral Methods in MATLAB*, ch. 3 |
| Spectral derivative (velocity) | Fortran `kh_velocity_from_psi`, `velocity_from_psi`; Pascal `VelFromPsi` | Trefethen (2000), ch. 3 |
| Advection term | Fortran `kh_nonlinear_rhs`, `vorticity_rhs`; Pascal `VortRHS` | Trefethen (2000), ch. 3; Orszag (1971) |
| ETDRK4 (alternate scheme) | Fortran `kh_etdrk4_precompute`, `kh_etdrk4_step` | Cox & Matthews (2002), J. Comput. Phys. 176, p. 430 |

The spectral anchors follow the references of `shared/canonical-equations/kh-numerical-methods.md`.

**Verification flag.** Every citation now records whether its section, equation and page have been checked against the printed book (`verified`). In this release none has been. The default audit reports this. `scripts/audit_citations.py --strict` fails until every citation is verified and nothing is pending.

## Development records and workstation names

Some code comments in v0.1.0–v0.1.2 cited development documents that are not part of this repository, by their private paths. Others named the development workstations. Both are listed under known limitations in `RELEASE-NOTES-v0.1.0.md`.

They are removed from 50 files:
- references to development documents now read "development record" with the document's subject and date;
- the numerical-methods references point to the public `shared/canonical-equations/kh-numerical-methods.md`;
- workstation names read "reference workstation", "natively" or similar.

The publish stage now fails on either, where v0.1.0–v0.1.2 only warned.

## Source changes

The comment and docstring changes are checked mechanically. After comments and docstrings are removed, the code of every changed Fortran, Pascal, C++, Rust, Scala, Python and Makefile source is identical to v0.1.2. Pascal braces nest, as in FPC's objfpc mode.

Text-only changes:
- `shared/api/openapi.yaml` (the description field);
- `KH-PHYSICS.md`;
- `shared/golden-vectors/README.md`;
- two PowerShell scripts.
