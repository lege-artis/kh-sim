! test_num_005_dealias.f90 — TC-NUM-KH-005: 2/3-rule de-aliasing correctness
!
! TC-NUM-KH-005 specification (PHYSICS-NUMERICAL-METHODS-v0.1.md §2.3):
!   Test target:  De-aliasing 2/3 rule
!   Method check: Fourier mode at k > 2N/3; should be zeroed by kh_dealias
!   Pass criterion: post-filter amplitude ≤ 1e-15
!
! Two sub-tests:
!   A) A mode exactly AT the cutoff (2/3·kmax) must be KEPT (amplitude unchanged).
!   B) A mode ABOVE the cutoff (> 2/3·kmax) must be ZEROED (amplitude ≤ 1e-15).
!
! Exit code: 0 = PASS, 1 = FAIL
!
! Compilation (DRY-RUN — validate on ThinkPad):
!   gfortran -O2 -o test_num_005 \
!       src/kh_constants.f90 src/kh_grid.f90 src/kh_fft.f90 \
!       src/kh_poisson.f90 src/kh_velocity.f90 src/kh_nonlinear.f90 \
!       tests/test_num_005_dealias.f90
!   ./test_num_005

program test_num_005_dealias
  use kh_constants
  use kh_grid,      only: kh_grid_make, kh_grid_wavenums, kh_grid_dealias_mask
  use kh_nonlinear, only: kh_dealias
  implicit none

  integer,        parameter :: NX = 64
  integer,        parameter :: NY = 32
  real(c_double), parameter :: LX = KH_LX_DEFAULT
  real(c_double), parameter :: LY = KH_LY_DEFAULT

  ! Mode indices to test
  ! kx_max index = NX/2 = 32.  2/3 cutoff = floor(2/3 * 32) = 21.
  ! Below cutoff: ix=20  →  kx1(20) = 2π*20/(NX*dx); 20 < 21  → kept.
  ! Above cutoff: ix=22  →  kx1(22); 22 > 21  → zeroed.
  integer, parameter :: IX_KEPT   = 20   ! inside 2/3 cutoff
  integer, parameter :: IX_ZEROED = 22   ! outside 2/3 cutoff

  real(c_double) :: x(0:NX-1), y(0:NY-1), dx, dy
  real(c_double) :: kx(0:NX-1, 0:NY-1), ky(0:NX-1, 0:NY-1), k2(0:NX-1, 0:NY-1)
  logical        :: mask(0:NX-1, 0:NY-1)
  complex(8)     :: a_hat(0:NX-1, 0:NY-1)
  real(c_double) :: amp_kept_before, amp_kept_after
  real(c_double) :: amp_zero_before, amp_zero_after
  logical        :: passed_a, passed_b

  real(c_double), parameter :: AMPLITUDE = 1.0_c_double
  real(c_double), parameter :: ZERO_TOL  = 1.0e-15_c_double

  ! Build grid and mask
  call kh_grid_make(NX, NY, LX, LY, x, y, dx, dy)
  call kh_grid_wavenums(NX, NY, dx, dy, kx, ky, k2)
  call kh_grid_dealias_mask(NX, NY, kx, ky, dx, dy, mask)

  ! ── Sub-test A: mode inside cutoff (ix=IX_KEPT, iy=0) ───────────────────────

  a_hat = cmplx(0.0_8, 0.0_8, kind=8)
  a_hat(IX_KEPT, 0) = cmplx(AMPLITUDE, 0.0_8, kind=8)
  amp_kept_before = abs(a_hat(IX_KEPT, 0))

  call kh_dealias(a_hat, mask, NX, NY)
  amp_kept_after = abs(a_hat(IX_KEPT, 0))

  ! Mode should survive: amplitude unchanged
  passed_a = (abs(amp_kept_after - amp_kept_before) <= ZERO_TOL)

  ! ── Sub-test B: mode above cutoff (ix=IX_ZEROED, iy=0) ──────────────────────

  a_hat = cmplx(0.0_8, 0.0_8, kind=8)
  a_hat(IX_ZEROED, 0) = cmplx(AMPLITUDE, 0.0_8, kind=8)
  amp_zero_before = abs(a_hat(IX_ZEROED, 0))

  call kh_dealias(a_hat, mask, NX, NY)
  amp_zero_after = abs(a_hat(IX_ZEROED, 0))

  ! Mode should be zeroed: amplitude ≤ 1e-15
  passed_b = (amp_zero_after <= ZERO_TOL)

  ! ── Report ───────────────────────────────────────────────────────────────────

  write(*, '(A)') "TC-NUM-KH-005: 2/3-rule de-aliasing"
  write(*, '(A,I0,A,I0)') "  Grid: ", NX, " x ", NY
  write(*, '(A,I0,A,I0)') "  kx_max index = ", NX/2, ";  2/3 cutoff index ≈ ", int(KH_DEALIAS_FACTOR * (NX/2))

  write(*, '(A)') "  Sub-test A: mode inside cutoff (should be KEPT)"
  write(*, '(A,I0,A,ES12.4)') "    ix = ", IX_KEPT, ";  amplitude before = ", amp_kept_before
  write(*, '(A,ES12.4,A)') "    amplitude after  = ", amp_kept_after, &
      merge("  PASS", "  FAIL", passed_a)

  write(*, '(A)') "  Sub-test B: mode above cutoff (should be ZEROED)"
  write(*, '(A,I0,A,ES12.4)') "    ix = ", IX_ZEROED, ";  amplitude before = ", amp_zero_before
  write(*, '(A,ES12.4,A,ES12.4,A)') "    amplitude after  = ", amp_zero_after, &
      "  (tol = ", ZERO_TOL, ")", merge("  PASS", "  FAIL", passed_b)

  if (passed_a .and. passed_b) then
    write(*, '(A)') "  RESULT: PASS"
    stop 0
  else
    write(*, '(A)') "  RESULT: FAIL"
    stop 1
  end if

end program test_num_005_dealias
