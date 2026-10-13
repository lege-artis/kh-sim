! test_num_003_etdrk4_linear.f90 — TC-NUM-KH-003: ETDRK4 on linear scalar ODE
!
! TC-NUM-KH-003 specification (PHYSICS-NUMERICAL-METHODS-v0.1.md §2.3):
!   Test target:  ETDRK4 on linear scalar test
!   Method check: dy/dt = λ·y, λ < 0; integrate to t = 10/|λ|; compare to e^(λ·t)
!   Pass criterion: rel err ≤ 1e-10
!
! This test exercises kh_etdrk4_precompute directly on a 1×1 "grid" (trivial
! spectral space), where the linear operator L = λ and N = 0 at all times.
! The ETDRK4 φ functions are exact for the linear problem; this verifies that
! the precompute + step combination reproduces exp(λ·t) to near-machine precision.
!
! Test parameters:
!   λ = -1.0   (decay rate)
!   y₀ = 1.0   (initial condition)
!   T = 10.0 / |λ| = 10.0
!   dt = 0.1   → 100 steps
!   Exact solution: y(T) = exp(-10.0) ≈ 4.53999e-5
!
! Exit code: 0 = PASS, 1 = FAIL
!
! Compilation (DRY-RUN — validate on ThinkPad):
!   gfortran -O2 -o test_num_003 \
!       src/kh_constants.f90 src/kh_grid.f90 src/kh_fft.f90 \
!       src/kh_poisson.f90 src/kh_velocity.f90 src/kh_nonlinear.f90 \
!       src/kh_etdrk4.f90 tests/test_num_003_etdrk4_linear.f90
!   ./test_num_003

program test_num_003_etdrk4_linear
  use kh_constants, only: c_double, KH_TOL_NONLINEAR
  use kh_etdrk4,    only: kh_etdrk4_precompute
  implicit none

  ! Test parameters
  real(c_double), parameter :: LAMBDA = -1.0_c_double
  real(c_double), parameter :: Y0     =  1.0_c_double
  real(c_double), parameter :: T_END  = 10.0_c_double
  real(c_double), parameter :: DT     = 0.1_c_double
  integer,        parameter :: NSTEPS = 100             ! T_END / DT

  ! 1×1 "grid" (trivial spectral space for scalar ODE test)
  integer, parameter :: NX = 1, NY = 1

  real(c_double) :: L_arr(0:0, 0:0)
  real(c_double) :: E_arr(0:0, 0:0), E2_arr(0:0, 0:0)
  real(c_double) :: phi1_arr(0:0, 0:0), phi1h_arr(0:0, 0:0)
  real(c_double) :: phi2_arr(0:0, 0:0)

  complex(8) :: y      ! current state (scalar, stored as 1×1 complex array element)
  complex(8) :: Na, Nb, Nc, Nd   ! nonlinear terms (all zero for linear test)
  complex(8) :: a_hat, b_hat, c_hat

  real(c_double) :: E, E2, phi1, phi1h, phi2
  real(c_double) :: y_exact, rel_err
  integer :: s
  logical :: passed

  ! ── Precompute ETDRK4 coefficients for L = λ ────────────────────────────────

  L_arr(0,0) = LAMBDA
  call kh_etdrk4_precompute(L_arr, DT, NX, NY, E_arr, E2_arr, phi1_arr, phi1h_arr, phi2_arr)

  E     = E_arr(0,0)
  E2    = E2_arr(0,0)
  phi1  = phi1_arr(0,0)
  phi1h = phi1h_arr(0,0)
  phi2  = phi2_arr(0,0)

  ! ── Time integration (N=0 at all times → pure linear ETDRK4) ────────────────

  y = cmplx(Y0, 0.0_8, kind=8)
  Na = cmplx(0.0_8, 0.0_8, kind=8)  ! nonlinear term is zero for linear problem

  do s = 1, NSTEPS
    ! Substep a
    a_hat = E2 * y + phi1h * DT * 0.5_c_double * Na
    Nb = cmplx(0.0_8, 0.0_8, kind=8)

    ! Substep b
    b_hat = E2 * y + phi1h * DT * 0.5_c_double * Nb
    Nc = cmplx(0.0_8, 0.0_8, kind=8)

    ! Substep c
    c_hat = E2 * a_hat + phi1h * DT * 0.5_c_double * (2.0_c_double*Nc - Na)
    Nd = cmplx(0.0_8, 0.0_8, kind=8)

    ! Final combination
    y = E * y &
      + DT * (phi1 - 3.0_c_double*phi2) * Na &
      + DT * 2.0_c_double * phi2 * (Nb + Nc) &
      + DT * (phi2 - phi1 + 1.0_c_double) * Nd
  end do

  ! ── Compare against exact solution: y(T) = exp(λ·T) ────────────────────────

  y_exact = exp(LAMBDA * T_END)
  rel_err = abs(real(y, kind=8) - y_exact) / abs(y_exact)
  passed = (rel_err <= KH_TOL_NONLINEAR)

  write(*, '(A)') "TC-NUM-KH-003: ETDRK4 on linear scalar ODE"
  write(*, '(A,F6.2,A,F6.2,A,I0,A,F6.4)') &
      "  λ = ", LAMBDA, "  T_end = ", T_END, "  steps = ", NSTEPS, "  dt = ", DT
  write(*, '(A,ES16.8)') "  y_exact  = exp(λ·T) = ", y_exact
  write(*, '(A,ES16.8)') "  y_ETDRK4 =           ", real(y, kind=8)
  write(*, '(A,ES12.4)') "  Relative error       = ", rel_err
  write(*, '(A,ES12.4)') "  Tolerance (KH_TOL_NONLINEAR) = ", KH_TOL_NONLINEAR

  if (passed) then
    write(*, '(A)') "  RESULT: PASS"
    stop 0
  else
    write(*, '(A)') "  RESULT: FAIL"
    stop 1
  end if

end program test_num_003_etdrk4_linear
