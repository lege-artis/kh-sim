! test_num_001_fft_roundtrip.f90 — TC-NUM-KH-001: FFT2 round-trip correctness
!
! TC-NUM-KH-001 specification (PHYSICS-NUMERICAL-METHODS-v0.1.md §2.3):
!   Test target:  FFT2 round-trip
!   Method check: IFFT2(FFT2(field)) == field
!   Pass criterion: ‖diff‖_∞ ≤ 1e-12 for a random complex field
!
! This test exercises kh_fft_forward_2d + kh_fft_inverse_2d from kh_fft.f90.
! It uses a deterministic pseudo-random field (LCG seed 42) to guarantee
! reproducible results without a platform RNG dependency.
!
! Exit code: 0 = PASS, 1 = FAIL
!
! Compilation (DRY-RUN — validate on ThinkPad with gfortran 13+):
!   gfortran -O2 -o test_num_001 \
!       src/kh_constants.f90 src/kh_grid.f90 src/kh_fft.f90 \
!       tests/test_num_001_fft_roundtrip.f90
!   ./test_num_001

program test_num_001_fft_roundtrip
  use kh_constants, only: KH_TOL_LINEAR
  use kh_fft,       only: kh_fft_forward_2d, kh_fft_inverse_2d
  implicit none

  ! ── Test parameters ──────────────────────────────────────────────────────────

  integer, parameter :: NX = 64
  integer, parameter :: NY = 32
  integer, parameter :: N_TOTAL = NX * NY

  ! ── Local variables ──────────────────────────────────────────────────────────

  complex(8) :: field_orig(0:NX-1, 0:NY-1)
  complex(8) :: field_work(0:NX-1, 0:NY-1)
  real(8)    :: diff_inf, diff_ij
  integer    :: i, j
  logical    :: passed

  ! LCG state for deterministic pseudo-random field (seed 42)
  integer(8) :: lcg_state

  ! ── Build deterministic pseudo-random complex field ──────────────────────────

  lcg_state = 42_8
  do j = 0, NY-1
    do i = 0, NX-1
      field_orig(i, j) = cmplx(lcg_next(lcg_state), lcg_next(lcg_state), kind=8)
    end do
  end do

  ! ── Copy to working array ────────────────────────────────────────────────────

  field_work = field_orig

  ! ── Forward FFT2 ─────────────────────────────────────────────────────────────

  call kh_fft_forward_2d(field_work, NX, NY)

  ! ── Inverse FFT2 ─────────────────────────────────────────────────────────────

  call kh_fft_inverse_2d(field_work, NX, NY)

  ! ── Compute ‖IFFT2(FFT2(field)) - field‖_∞ ──────────────────────────────────

  diff_inf = 0.0_8
  do j = 0, NY-1
    do i = 0, NX-1
      diff_ij = abs(field_work(i, j) - field_orig(i, j))
      if (diff_ij > diff_inf) diff_inf = diff_ij
    end do
  end do

  ! ── Evaluate pass / fail ─────────────────────────────────────────────────────

  passed = (diff_inf <= KH_TOL_LINEAR)

  write(*, '(A)') "TC-NUM-KH-001: FFT2 round-trip"
  write(*, '(A,I0,A,I0,A)') "  Grid: ", NX, " x ", NY, " (powers of 2)"
  write(*, '(A,ES12.4)') "  ‖IFFT2(FFT2(f)) - f‖_∞ = ", diff_inf
  write(*, '(A,ES12.4)') "  Tolerance (KH_TOL_LINEAR) = ", KH_TOL_LINEAR

  if (passed) then
    write(*, '(A)') "  RESULT: PASS"
    stop 0
  else
    write(*, '(A)') "  RESULT: FAIL"
    stop 1
  end if

contains

  ! ── Simple LCG pseudo-random number generator ────────────────────────────────
  !
  ! Returns values in (-1, 1) using the LCG:
  !   state = (1664525 * state + 1013904223) mod 2^32
  !
  ! Not suitable for physics — used only for reproducible test field generation.

  real(8) function lcg_next(state)
    integer(8), intent(inout) :: state
    integer(8), parameter :: A = 1664525_8
    integer(8), parameter :: C = 1013904223_8
    integer(8), parameter :: M = 4294967296_8   ! 2^32
    state = mod(A * state + C, M)
    ! Map to (-1, 1)
    lcg_next = (real(state, 8) / real(M, 8)) * 2.0_8 - 1.0_8
  end function lcg_next

end program test_num_001_fft_roundtrip
