# kh-sim canonical equation — linear stability dispersion (KH-LS)

**Status:** canonical (TeX-grade)
**Version:** 0.1 — first cut, 2026-05-12
**Companion to:** `kh-linear-stability.tex`
**Cited by:** future property test (phase 7) anchoring the simulated early-time
growth rate against the analytic dispersion relation.

This document states the **linear stability dispersion relation** σ(k) for the
kh-sim canonical initial conditions (tanh shear profile). It is the analytical
anchor for one of the property tests in the phase-7 cross-backend property suite:
*"simulated early-time growth rate of a small-amplitude perturbation at mode k
must match the analytic σ(k) within tolerance."*

This is the most demanding property test in the suite — it certifies that the
implementation correctly solves the linearised Navier-Stokes equations, not just
that it conserves a few invariants.

---

## KH-LS-1: tanh-profile inviscid dispersion relation

Consider the inviscid limit (Re → ∞) of the kh-sim governing equations
(KH-GE-1 with ν = 0). For an unbounded domain with base-flow velocity profile

$$
U(y) = U_0\,\tanh\!\left(\frac{y - L_y/2}{\delta}\right),
$$

(corresponding to the kh-sim canonical IC, Eq. KH-GE-4), linearising about the
base flow and seeking normal-mode solutions of the form

$$
\omega'(x, y, t) = \tilde\omega(y)\,e^{i k_x x + \sigma t}
$$

gives the **Rayleigh stability equation** (Drazin & Reid §4):

$$
(U - c)\left(\frac{d^2 \tilde\psi}{dy^2} - k_x^2\,\tilde\psi\right) - U''\,\tilde\psi = 0
$$

where c = -iσ/k_x is the complex phase speed.

For the tanh profile specifically, the dispersion relation (Drazin & Reid §4.5)
admits an analytic solution. The fastest-growing mode has wavenumber:

$$
k_{x,\max}\,\delta \approx 0.4446
$$

and corresponding growth rate:

$$
\sigma_{\max} \approx 0.19\,\frac{U_0}{\delta}.
\quad\text{(KH-LS-1)}
$$

For modes with k_xδ < 1 (wavelength > 2πδ), σ(k_x) is real and positive
(unstable). For k_xδ > 1, σ(k_x) is purely imaginary (neutrally stable
oscillation).

## KH-LS-2: kh-sim canonical mode

The kh-sim canonical IC excites perturbation mode k = 2 in the x-direction over
domain length L_x = 1.0, giving physical wavenumber:

$$
k_x = \frac{2\pi k}{L_x} = \frac{2\pi \cdot 2}{1.0} = 4\pi \approx 12.566\;\text{rad/m}.
$$

With δ = 0.05·Ly = 0.025 m and U_0 = 1 m/s, the non-dimensional combination is

$$
k_x\,\delta = 12.566 \cdot 0.025 \approx 0.314.
$$

This is comfortably below the marginal-stability cutoff k_xδ = 1, so the mode
is in the **unstable** band and should grow exponentially at early times.

The expected growth rate σ(k_x = 4π) for k_xδ ≈ 0.314 can be read off the
tanh-profile dispersion curve (Drazin & Reid Fig. 4.4 or equivalent). For
k_xδ near 0.3, σδ/U_0 ≈ 0.18; thus:

$$
\sigma_{\text{expected}} \approx 0.18 \cdot \frac{U_0}{\delta} = 0.18 \cdot \frac{1}{0.025} = 7.2\;\text{s}^{-1}.
\quad\text{(KH-LS-2)}
$$

(Numerical value provisional — to be confirmed in phase 7 with explicit
numerical solve of the Rayleigh equation for the kh-sim default parameters.)

## KH-LS-3: property test gate (phase 7)

For a small-amplitude perturbation initial condition (A = 1e-4 instead of
the canonical A = 0.01, to keep the perturbation in the linear regime), the
simulated kinetic-energy growth rate over early times should match the
analytic σ(k_x):

$$
\frac{d\,\ln E_k(t)}{d t}\bigg|_{t \to 0^+} \approx 2\,\sigma_{\text{expected}}
\quad\text{(KH-LS-3)}
$$

(factor of 2 because energy scales as amplitude squared).

**Tolerance gate** for the property test:

- Linear regime: A·t < 0.1 (perturbation has not entered non-linear roll-up yet)
- Match: |σ_simulated - σ_analytic| / σ_analytic < 5% for canonical resolution (128 × 64)
- Match: < 1% at well-resolved settings (256 × 128 + sub-CFL Δt)

This is **the canonical property test for kh-sim**. It anchors physical
correctness, not just numerical-method correctness. Phase 7 implementation:
a new test in each backend (NumPy reference + Rust + Scala + C++ + Fortran +
Pascal) called `test_kh_linear_growth_rate.<ext>`.

---

## Outside scope of this document

- **Viscous corrections** to σ(k_x) for finite Re — covered in Drazin & Reid §4.6
  and Schmid & Henningson §3.4. Important for the kh-sim default Re = 1000
  (small but non-zero correction).
- **Non-linear billow roll-up dynamics** — the linear theory captures only the
  initial exponential growth; the visible "rolling vortex structure" of KH
  is a non-linear saturation phenomenon (vortex sheet roll-up, pairing).
- **Density-stratified KH** (Chandrasekhar, Schmid & Henningson §3.5).
  kh-sim assumes uniform density.

---

## References

See `shared/reference-bibliography/refs.bib`:
- Drazin & Reid (2004) §4.5 — tanh-profile dispersion relation
- Schmid & Henningson (2001) §3 — modern numerical treatment
- Chandrasekhar (1981) — classical stratified-fluid treatment
