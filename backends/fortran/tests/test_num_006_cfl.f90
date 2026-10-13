! test_num_006_cfl.f90 — TC-NUM-KH-006: solver stability across Reynolds numbers
!
! TC-NUM-KH-006 specification (PHYSICS-NUMERICAL-METHODS-v0.1.md §2.3):
!   Test target:  Full solver stability with CFL monitoring
!   Method check: Run kh_solver_run at Re=100, 1000, 10000; T=0.1; 64×32 grid.
!                 Verify solution remains bounded and divergence-free.
!   Pass criterion (per sub-test):
!     A) KE > 0 and KE < 10.0  (finite, positive, not exploded)
!     B) div_rms ≤ KH_TOL_DIVERGENCE (1e-10)
!     C) cfl_peak ≤ KH_CFL_FACTOR (0.4)  [i.e., within stable CFL band]
!
! Physical expectations:
!   Re=100   (ν=0.01):  strong viscous dissipation → KE decreases, enstrophy damped
!   Re=1000  (ν=0.001): moderate viscosity → reference regime (standard KH)
!   Re=10000 (ν=0.0001): nearly inviscid → more nonlinear, higher enstrophy growth
!
!   In all cases, with U₀=1.0, dx≈0.0156, dt=0.001:
!     CFL ≈ 0.001 × (1.0/0.0156 + 0.5/0.0156) ≈ 0.001 × 96 ≈ 0.096
!   This is well within the KH_CFL_FACTOR=0.4 safety limit.
!
! Grid: 64×32.  T = 0.1 (100 steps at dt=0.001).  No JSON output (out_unit=-1).
!
! Exit code: 0 = PASS (all 3 sub-tests), 1 = FAIL
!
! Compilation (DRY-RUN — validate on ThinkPad):
!   gfortran -O2 -o test_num_006 \
!       src/kh_constants.f90 src/kh_fft.f90 src/kh_grid.f90 \
!       src/kh_initial_conditions.f90 src/kh_nonlinear.f90 \
!       src/kh_diagnostics.f90 src/kh_physics.f90 \
!       tests/test_num_006_cfl.f90
!   ./test_num_006

program test_num_006_cfl
  use kh_constants
  use kh_physics, only: kh_simulate_rk4
  implicit none

  ! ── Test parameters ──────────────────────────────────────────────────────────

  integer,        parameter :: NX       = 64
  integer,        parameter :: NY       = 32
  real(c_double), parameter :: LX       = KH_LX_DEFAULT   ! 1.0
  real(c_double), parameter :: LY       = KH_LY_DEFAULT   ! 0.5
  real(c_double), parameter :: DT       = KH_DT_DEFAULT   ! 0.001
  integer,        parameter :: NSTEPS   = 100              ! T = 0.1
  real(c_double), parameter :: T_END    = DT * real(NSTEPS, c_double)
  real(c_double), parameter :: AMP      = KH_AMP_DEFAULT  ! 0.01
  integer,        parameter :: KMODE    = KH_MODE_DEFAULT  ! 2

  ! Pass criteria
  real(c_double), parameter :: KE_MIN   = 0.0_c_double
  real(c_double), parameter :: KE_MAX   = 10.0_c_double   ! explosion guard
  real(c_double), parameter :: CFL_PASS = KH_CFL_FACTOR   ! 0.4

  ! Three test Reynolds numbers
  integer,        parameter :: N_RE     = 3
  real(c_double)            :: re_vals(N_RE)
  character(len=12)         :: re_labels(N_RE)

  ! ── Per-run results ───────────────────────────────────────────────────────────

  real(c_double) :: ke(N_RE), enstrophy(N_RE), max_vort(N_RE), div_rms(N_RE)
  real(c_double) :: cfl_peak(N_RE), nu

  logical :: passed_ke(N_RE), passed_div(N_RE), passed_cfl(N_RE)
  logical :: all_pass

  integer :: r

  ! ── Initialise test matrix ───────────────────────────────────────────────────

  re_vals(1) = 100.0_c_double;   re_labels(1) = "Re=100     "
  re_vals(2) = 1000.0_c_double;  re_labels(2) = "Re=1000    "
  re_vals(3) = 10000.0_c_double; re_labels(3) = "Re=10000   "

  ! ── Header ───────────────────────────────────────────────────────────────────

  write(*, '(A)') "TC-NUM-KH-006: Solver stability across Reynolds numbers"
  write(*, '(A,I0,A,I0,A,F6.4,A,I0,A,F5.3)') &
      "  Grid: ", NX, " x ", NY, &
      "  dt = ", DT, "  steps = ", NSTEPS, "  T = ", T_END

  ! ── Run solver at each Re ─────────────────────────────────────────────────────

  do r = 1, N_RE
    nu = KH_U0_DEFAULT / re_vals(r)

    write(*, '(A,A,A,ES10.3)') "  Running ", trim(re_labels(r)), "  nu = ", nu

    call kh_simulate_rk4(NX, NY, LX, LY, nu, DT, NSTEPS, AMP, KMODE, &
                          ke(r), enstrophy(r), max_vort(r), div_rms(r), cfl_peak(r))

    ! Evaluate pass/fail for this Re
    passed_ke(r)  = (ke(r) > KE_MIN) .and. (ke(r) < KE_MAX)
    passed_div(r) = (div_rms(r) <= KH_TOL_DIVERGENCE)
    passed_cfl(r) = (cfl_peak(r) <= CFL_PASS)
  end do

  ! ── Report ───────────────────────────────────────────────────────────────────

  write(*, '(A)') ""
  write(*, '(A)') "  Results:"
  write(*, '(A)') "  Re        KE           Enstrophy    max|ω|       div_rms      CFL_peak  KE  div  CFL"
  write(*, '(A)') "  ------  -----------  -----------  -----------  -----------  ---------  --  ---  ---"

  do r = 1, N_RE
    write(*, '(A,A,5(2X,ES11.4),2X,A,2X,A,2X,A)') &
        "  ", trim(re_labels(r)), &
        ke(r), enstrophy(r), max_vort(r), div_rms(r), cfl_peak(r), &
        merge("OK", "!!", passed_ke(r)),  &
        merge("OK ", "!! ", passed_div(r)), &
        merge("OK ", "!! ", passed_cfl(r))
  end do

  ! ── Physical consistency cross-check: Re↑ → less dissipation ─────────────────

  write(*, '(A)') ""
  write(*, '(A)') "  Physical consistency (Re increases → more enstrophy expected):"
  if (enstrophy(2) >= enstrophy(1)) then
    write(*, '(A,ES10.3,A,ES10.3,A)') &
        "    Enstrophy Re=1000 (", enstrophy(2), ") >= Re=100  (", enstrophy(1), ")  OK"
  else
    write(*, '(A)') "    WARNING: Re=1000 enstrophy < Re=100 (unexpected)"
  end if
  if (enstrophy(3) >= enstrophy(2)) then
    write(*, '(A,ES10.3,A,ES10.3,A)') &
        "    Enstrophy Re=10000(", enstrophy(3), ") >= Re=1000 (", enstrophy(2), ")  OK"
  else
    write(*, '(A)') "    WARNING: Re=10000 enstrophy < Re=1000 (unexpected)"
  end if

  ! ── Final verdict ─────────────────────────────────────────────────────────────

  all_pass = all(passed_ke) .and. all(passed_div) .and. all(passed_cfl)

  write(*, '(A)') ""
  if (all_pass) then
    write(*, '(A)') "  RESULT: PASS"
    stop 0
  else
    write(*, '(A)') "  RESULT: FAIL"
    stop 1
  end if

end program test_num_006_cfl
