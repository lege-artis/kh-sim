! kh_fft.f90 — KH instability: 2D FFT module
!
! REFERENCE WATERMARK — kh-sim Fortran reference, v0.1;
! canonical for: TC-NUM-KH-001..008;
! do not modify without bumping reference version + regenerating reference outputs.
!
! Wraps the hand-rolled Cooley-Tukey radix-2 DIT FFT from the monolithic
! kh_physics.f90 into a clean public API.
!
! OQ-NUM-01 status (NM development spec §9): FFT library choice (FFTW3
! vs fftpack5) is deferred — this module provides the same algorithm currently
! in kh_physics.f90 (hand-rolled, no external dependency) so NUM-KH-FOR-01
! compiles without any system library. When OQ-NUM-01 is resolved, swap only
! fft1d_inplace below; all callers of kh_fft_forward_2d / kh_fft_inverse_2d
! remain unchanged.
!
! Public interface:
!   kh_fft_forward_2d(a, nx, ny)   — in-place forward 2D FFT
!   kh_fft_inverse_2d(a, nx, ny)   — in-place inverse 2D FFT (normalised)
!
! Convention:
!   - nx and ny MUST be powers of 2.
!   - Array a(0:nx-1, 0:ny-1) is complex(8) column-major.
!   - Forward FFT: follows standard DFT sign convention (e^{-2πi k n/N}).
!   - Inverse FFT: applies 1/(nx*ny) normalisation (unitary-style).
!   - Row-then-column: forward/inverse FFT along axis-1 (ny) then axis-0 (nx).
!     This matches numpy.fft.fft2 row-major semantics when the Fortran array
!     is viewed as transposed (ix is the "slow" axis, iy the "fast" axis).
!
! Reference: kh-sim/shared/physics/KH-PHYSICS.md §4 (FFT spatial discretisation)
! Numerical methods: shared/canonical-equations/kh-numerical-methods.md

module kh_fft
  use kh_constants, only: KH_PI
  implicit none
  private

  public :: kh_fft_forward_2d
  public :: kh_fft_inverse_2d

contains

  ! ── 1D in-place Cooley-Tukey FFT (radix-2, DIT) ────────────────────────────
  !
  ! Identical algorithm to fft1d in kh_physics.f90, isolated here so it can
  ! be replaced by an FFTW3 call when OQ-NUM-01 is resolved without touching
  ! the 2D wrappers.
  !
  !   a(0:n-1)  [inout]  complex array, n must be a power of 2
  !   n         [in]     array length
  !   inverse   [in]     .false. = forward DFT; .true. = inverse (normalised)
  !
  ! Reference: Press, Teukolsky, Vetterling & Flannery (2007), Numerical Recipes,
  !            3rd ed., S12.2, eq. 12.2.2 (decimation-in-time butterfly recurrence).
  ! Sign convention: forward = exp(-2*pi*i*k*n/N), per S12.1 (standard); the
  !            inverse applies 1/n normalisation. Phase-6c sign fix: commit 8a92f8e.
  ! Complexity: O(N log N) -- log2(N) stages of N/2 butterflies; checked
  !            empirically by tests/test_complexity_c4.f90 (doctrine C4).

  subroutine fft1d_inplace(a, n, inverse)
    integer,    intent(in)    :: n
    complex(8), intent(inout) :: a(0:n-1)
    logical,    intent(in)    :: inverse

    integer    :: i, j, bit, len, k
    complex(8) :: wlen, w, u, v
    real(8)    :: ang

    ! bit-reversal permutation
    j = 0
    do i = 1, n-1
      bit = n / 2
      do while (iand(j, bit) /= 0)
        j = ieor(j, bit)
        bit = bit / 2
      end do
      j = ieor(j, bit)
      if (i < j) then
        u = a(i); a(i) = a(j); a(j) = u
      end if
    end do

    ! butterfly stages
    len = 2
    do while (len <= n)
      ang = -2.0_8 * KH_PI / real(len, 8)
      if (inverse) ang = -ang
      wlen = cmplx(cos(ang), sin(ang), kind=8)
      i = 0
      do while (i < n)
        w = cmplx(1.0_8, 0.0_8, kind=8)
        do k = 0, len/2 - 1
          u = a(i + k)
          v = a(i + k + len/2) * w
          a(i + k)         = u + v
          a(i + k + len/2) = u - v
          w = w * wlen
        end do
        i = i + len
      end do
      len = len * 2
    end do

    ! normalise for inverse (1/n per axis; 2D normalisation is 1/(nx*ny) applied
    ! once per axis so each axis divides by its own n)
    if (inverse) then
      a = a / real(n, 8)
    end if
  end subroutine fft1d_inplace


  ! ── 2D forward FFT ───────────────────────────────────────────────────────────
  !
  ! In-place forward 2D FFT of a(0:nx-1, 0:ny-1).
  ! Process: FFT along axis-1 (each row of length ny), then axis-0 (each col nx).
  !
  ! Arguments:
  !   a(0:nx-1, 0:ny-1) [inout]  complex 2D array
  !   nx, ny             [in]    dimensions (powers of 2)

  subroutine kh_fft_forward_2d(a, nx, ny)
    integer,    intent(in)    :: nx, ny
    complex(8), intent(inout) :: a(0:nx-1, 0:ny-1)

    integer    :: i, j
    complex(8) :: row(0:ny-1), col(0:nx-1)

    ! Forward FFT along axis-1 (rows of length ny)
    do i = 0, nx-1
      do j = 0, ny-1; row(j) = a(i, j); end do
      call fft1d_inplace(row, ny, .false.)
      do j = 0, ny-1; a(i, j) = row(j); end do
    end do

    ! Forward FFT along axis-0 (columns of length nx)
    do j = 0, ny-1
      do i = 0, nx-1; col(i) = a(i, j); end do
      call fft1d_inplace(col, nx, .false.)
      do i = 0, nx-1; a(i, j) = col(i); end do
    end do
  end subroutine kh_fft_forward_2d


  ! ── 2D inverse FFT ───────────────────────────────────────────────────────────
  !
  ! In-place inverse 2D FFT of a(0:nx-1, 0:ny-1).
  ! Normalisation: 1/(nx*ny) applied as 1/ny per row pass and 1/nx per col pass.
  !
  ! Arguments:
  !   a(0:nx-1, 0:ny-1) [inout]  complex 2D array (spectral input)
  !   nx, ny             [in]    dimensions (powers of 2)

  subroutine kh_fft_inverse_2d(a, nx, ny)
    integer,    intent(in)    :: nx, ny
    complex(8), intent(inout) :: a(0:nx-1, 0:ny-1)

    integer    :: i, j
    complex(8) :: row(0:ny-1), col(0:nx-1)

    ! Inverse FFT along axis-1 (rows of length ny) — divides by ny
    do i = 0, nx-1
      do j = 0, ny-1; row(j) = a(i, j); end do
      call fft1d_inplace(row, ny, .true.)
      do j = 0, ny-1; a(i, j) = row(j); end do
    end do

    ! Inverse FFT along axis-0 (columns of length nx) — divides by nx
    do j = 0, ny-1
      do i = 0, nx-1; col(i) = a(i, j); end do
      call fft1d_inplace(col, nx, .true.)
      do i = 0, nx-1; a(i, j) = col(i); end do
    end do
  end subroutine kh_fft_inverse_2d

end module kh_fft
