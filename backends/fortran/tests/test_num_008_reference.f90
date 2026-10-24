! test_num_008_reference.f90 — TC-NUM-KH-008: canonical reference run validation
!
! TC-NUM-KH-008 specification (NM development spec §2.3):
!   Test target:  Canonical reference run output
!   Method check: Run kh_solver_run with the canonical parameters defined in
!                 kh_reference.f90; compare final diagnostics against the
!                 stored KH_REF_* constants in kh_constants.f90.
!   Pass criterion:
!     A) |KE − KH_REF_KE|           / KH_REF_KE           ≤ KH_REF_TOL  (5%)
!     B) |Enstrophy − KH_REF_ENSTROPHY| / KH_REF_ENSTROPHY ≤ KH_REF_TOL  (5%)
!     C) |max|ω| − KH_REF_MAX_VORT| / KH_REF_MAX_VORT      ≤ KH_REF_TOL  (5%)
!     D) div_rms ≤ KH_TOL_DIVERGENCE (1e-10)
!
! Canonical run: NX=64, NY=32, Re=1000, dt=0.001, T=0.1 (100 steps).
! The 5% tolerance accommodates alternative FFT libraries, compilation flags,
! and minor algorithmic variations while still certifying physical correctness.
! Byte-exact sha256 reproducibility is deferred to OQ-NUM-05.
!
! Expected outcome based on the reference-workstation gfortran validation (2026-05-04):
!   KE=1.1281e-1 vs ref 1.1281e-1 → rel err ≈ 0.003%  PASS
!   Enstrophy=4.3703e+1 vs ref 4.3705e+1 → rel err ≈ 0.005%  PASS
!   max|ω|=3.1986e+1 vs ref 3.1908e+1 → rel err ≈ 0.25%  PASS
!   div_rms=1.50e-14  PASS
!
! Exit code: 0 = PASS, 1 = FAIL
!
! Compilation (DRY-RUN — validate natively):
!   gfortran -O2 -o test_num_008 \
!       src/kh_constants.f90 src/kh_fft.f90 src/kh_grid.f90 \
!       src/kh_initial_conditions.f90 src/kh_nonlinear.f90 \
!       src/kh_diagnostics.f90 src/kh_physics.f90 \
!       src/kh_reference.f90 \
!       tests/test_num_008_reference.f90
!   ./test_num_008

program test_num_008_reference
  use kh_constants
  use kh_physics,   only: kh_simulate_rk4
  use kh_reference, only: KH_REF_NX, KH_REF_NY, KH_REF_LX, KH_REF_LY, &
                           KH_REF_NU_RUN, KH_REF_DT, KH_REF_NSTEPS,    &
                           KH_REF_AMP, KH_REF_MODE, KH_REF_RE_RUN,      &
                           kh_reference_compare
  implicit none

  ! ── Solver outputs ────────────────────────────────────────────────────────────

  real(c_double) :: ke, enstrophy, max_vort, div_rms, cfl_peak

  ! ── Comparison results ────────────────────────────────────────────────────────

  logical        :: passed_ke, passed_ens, passed_vort, passed_div
  real(c_double) :: rel_ke, rel_ens, rel_vort
  logical        :: all_pass

  ! ── Run canonical reference simulation ───────────────────────────────────────

  write(*, '(A)') "TC-NUM-KH-008: Canonical reference run validation"
  write(*, '(A,I0,A,I0,A,F8.4,A,F8.4)') &
      "  Grid: ", KH_REF_NX, " x ", KH_REF_NY, &
      "  Lx = ", KH_REF_LX, "  Ly = ", KH_REF_LY
  write(*, '(A,F10.1,A,ES10.3)') &
      "  Re = ", KH_REF_RE_RUN, "   nu = ", KH_REF_NU_RUN
  write(*, '(A,ES10.3,A,I0,A,F6.3)') &
      "  dt = ", KH_REF_DT, "  steps = ", KH_REF_NSTEPS, &
      "  T = ", KH_REF_DT * real(KH_REF_NSTEPS, c_double)

  call kh_simulate_rk4(KH_REF_NX, KH_REF_NY, KH_REF_LX, KH_REF_LY, &
                        KH_REF_NU_RUN, KH_REF_DT, KH_REF_NSTEPS,    &
                        KH_REF_AMP, KH_REF_MODE,                     &
                        ke, enstrophy, max_vort, div_rms, cfl_peak)

  ! ── Compare against stored reference values ───────────────────────────────────

  call kh_reference_compare(ke, enstrophy, max_vort, div_rms, &
                             passed_ke, passed_ens, passed_vort, passed_div, &
                             rel_ke, rel_ens, rel_vort)

  ! ── Report ───────────────────────────────────────────────────────────────────

  write(*, '(A)') "  --- Diagnostic comparison (tolerance = 5%) ---"
  write(*, '(A)') "  Quantity      Computed       Reference      Rel err     Status"
  write(*, '(A)') "  ----------  -----------    -----------    ----------    ------"

  write(*, '(A,ES13.6,A,ES13.6,A,ES10.3,A)') &
      "  KE          ", ke,         "    ", KH_REF_KE,         "    ", rel_ke, &
      merge("    PASS", "    FAIL", passed_ke)

  write(*, '(A,ES13.6,A,ES13.6,A,ES10.3,A)') &
      "  Enstrophy   ", enstrophy,  "    ", KH_REF_ENSTROPHY,  "    ", rel_ens, &
      merge("    PASS", "    FAIL", passed_ens)

  write(*, '(A,ES13.6,A,ES13.6,A,ES10.3,A)') &
      "  max|omega|  ", max_vort,   "    ", KH_REF_MAX_VORT,   "    ", rel_vort, &
      merge("    PASS", "    FAIL", passed_vort)

  write(*, '(A,ES13.6,A,ES13.6,A,ES10.3,A)') &
      "  div_rms     ", div_rms,    "    ", KH_TOL_DIVERGENCE,  "    ", div_rms, &
      merge("    PASS", "    FAIL", passed_div)

  write(*, '(A,ES12.4)') "  Peak CFL = ", cfl_peak

  all_pass = passed_ke .and. passed_ens .and. passed_vort .and. passed_div

  if (all_pass) then
    write(*, '(A)') "  RESULT: PASS"
    stop 0
  else
    write(*, '(A)') "  RESULT: FAIL"
    stop 1
  end if

end program test_num_008_reference
