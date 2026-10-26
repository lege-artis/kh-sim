# lege-artis/kh-sim — v0.1.5

**Tag:** `v0.1.5`
**Type:** a second stage of the turbulence testbed CP-7F-B, which reaches the asymptotic enstrophy cascade. Also a fix in the Rust backend: it now reads every request the contract allows.

## CP-7F-B2: the asymptotic enstrophy cascade

Until v0.1.4, testbed B measured slopes of −4.12 … −4.16. That is the finite-Reynolds coherent-vortex regime (McWilliams, 1984), not Kraichnan's k⁻³.

**Why B1 cannot reach k⁻³.** At its parameters (forcing at k 10–12, 256², Re = 5000), the dissipation wavenumber is k_d = (η_ν/ν³)^{1/6} ≈ 83:
- the fit range k 16–64 lies mostly in the viscous range;
- no range is left between 2k_f and k_d/3;
- the enstrophy flux falls from 1.0 to 0.46 of the viscous dissipation across the fit.

B1 stays as it is: a regime and reproducibility test with a per-backend spin-up. A scan of eight 512² configurations found the cascade once the forcing moved to k 4–6.

**B2 configuration:**
- **Grid and physics:** L = 2π, 512², Re = 40 000, band forcing A = 2 on 4 ≤ |k| ≤ 6, rotating phases `forcing_omega0` 5, linear drag 0.01, dt 0.001.
- **Shared spin-up:** the oracle spins up once, to t = 440 through 128², 256² and 512², and every target starts from that field.
- **Window:** each target then runs 40 requests of 500 steps and records one spectrum per request.

**What is measured, on k 12–44 (2k_f … k_d/3):**
- the slope of the mean spectrum;
- the enstrophy flux Π_Z(k) = Σ_{q>k} T_Z(q) relative to the viscous dissipation η_ν = 2νP, where T_Z is the transfer of the dealiased nonlinear term, computed by the harness from each returned field;
- the stationarity between the two 10-unit halves;
- the flatness of the log-compensated spectrum k³E(k)[ln(k/k_f)]^{1/3}, reported.

| | oracle reference | all six implementations | gate (reference ± 3σ) |
|---|---|---|---|
| slope | −3.492 | −3.5127 | [−3.64, −3.35] |
| stationarity | 0.107 | 0.060 | ≤ 0.20 |
| Π_Z/η_ν, mean | 0.938 | 1.047 | [0.82, 1.06] |
| Π_Z/η_ν, coefficient of variation | 0.064 | 0.092 | ≤ 0.12 |
| k_d | 131 | 129 | — |

**Interpretation.** Kraichnan's log-corrected spectrum k⁻³[ln(k/k_f)]^{−1/3} (1971) reads −3.22 on this range. The remaining 0.27–0.29 is the steepening by the linear drag (Boffetta et al., 2002). Across the scan it grows with the drag parameter: ξ ≈ 0.10 + 5.3·α·η^{−1/3}.

**How the gate is set.** Each σ is the larger of two estimates:
- a block bootstrap of the 40 spectra;
- the spread of three scan runs at the same drag parameter. A 20-unit window cannot sample the condensate's time scale 1/α = 100.

**Why the six implementations agree to four decimals.** They share one spin-up. Over the 20-unit window their round-off differences grow only from about 1e-15 to about 1e-9, measured in the peak Courant number of the last request. B1, spun up per implementation over 110 time units, still separates them.

**Runtime.** B2 runs only with `KH_CP7FB2=1`: about 2 h of spin-up, then the targets; Pascal takes about 10 h.

## Rust reads every request the contract allows

The contract allows grids up to 4 194 304 cells. A 512² request with `initial_omega` is a 5.4 MB body, but the Rust backend read at most 2 MiB (axum's default) and answered 413 in plain text before the contract was consulted. The B2 run found this with its first 512² request; the B1 bodies at 256² are 1.3 MB.

Rust now:
- reads a body up to a cell-cap request with `initial_omega` (4 194 304 × 32 bytes + 1 MiB);
- answers a larger body with the contract's 400, `field: "body"`.

The other four backends read any size.

`test_conformance.py` §7 pins it: a 512² request with `initial_omega` must be computed as the oracle computes it.
