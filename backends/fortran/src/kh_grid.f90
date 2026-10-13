! kh_grid.f90 — KH instability: grid coordinates and wavenumber arrays
!
! REFERENCE WATERMARK — kh-sim Fortran reference, v0.1;
! canonical for: TC-NUM-KH-001..008;
! do not modify without bumping reference version + regenerating reference outputs.
!
! Provides:
!   kh_grid_make      — allocate and fill physical-space grid coords (x, y)
!   kh_grid_wavenums  — fill spectral wavenumber arrays (kx, ky, k2)
!   kh_grid_dealias_mask — 2/3-rule Boolean mask over wavenumber space
!
! Grid convention (matches kh_physics.f90 and KH-PHYSICS.md §4):
!   Fortran column-major arrays, index (ix, iy) where:
!     ix = 0 .. nx-1  (x-direction)
!     iy = 0 .. ny-1  (y-direction)
!   Physical spacing:  dx = Lx/nx,  dy = Ly/ny
!   FFT wavenumbers:   kx(ix) = 2π·fftfreq(ix, nx) / dx  (angular, rad/m)
!
! Reference: kh-sim/shared/physics/KH-PHYSICS.md §4
! Numerical methods: _config/PHYSICS-NUMERICAL-METHODS-v0.1.md §2.2

module kh_grid
  use kh_constants
  implicit none
  private

  public :: kh_grid_make
  public :: kh_grid_wavenums
  public :: kh_grid_dealias_mask

contains

  ! ── Physical-space grid coordinates ─────────────────────────────────────────
  !
  ! x(ix) = ix * dx,  y(iy) = iy * dy
  ! Outputs are 1D arrays of length nx and ny respectively.
  !
  ! Arguments:
  !   nx, ny    [in]  grid resolution
  !   lx, ly    [in]  domain dimensions
  !   x(0:nx-1) [out] x coordinates
  !   y(0:ny-1) [out] y coordinates
  !   dx        [out] grid spacing in x
  !   dy        [out] grid spacing in y

  subroutine kh_grid_make(nx, ny, lx, ly, x, y, dx, dy)
    integer,          intent(in)  :: nx, ny
    real(c_double),   intent(in)  :: lx, ly
    real(c_double),   intent(out) :: x(0:nx-1), y(0:ny-1)
    real(c_double),   intent(out) :: dx, dy
    integer :: i

    dx = lx / real(nx, c_double)
    dy = ly / real(ny, c_double)

    do i = 0, nx-1
      x(i) = real(i, c_double) * dx
    end do
    do i = 0, ny-1
      y(i) = real(i, c_double) * dy
    end do
  end subroutine kh_grid_make


  ! ── Spectral wavenumber arrays ───────────────────────────────────────────────
  !
  ! Angular wavenumber (rad/unit_length) following the FFT fftfreq convention:
  !   kx1(i) = 2π * i / (nx * dx)      for i = 0 .. nx/2
  !   kx1(i) = 2π * (i - nx) / (nx * dx) for i = nx/2+1 .. nx-1
  !
  ! 2D arrays kx(ix, iy) and ky(ix, iy) replicate the 1D vectors across each
  ! dimension, and k2(ix, iy) = kx^2 + ky^2 (k^2 for the Poisson solve).
  !
  ! Arguments:
  !   nx, ny          [in]  grid resolution
  !   dx, dy          [in]  grid spacing
  !   kx(0:nx-1,0:ny-1) [out] x-wavenumber 2D array
  !   ky(0:nx-1,0:ny-1) [out] y-wavenumber 2D array
  !   k2(0:nx-1,0:ny-1) [out] kx^2 + ky^2; k2(0,0) = 1 (Poisson convention)

  subroutine kh_grid_wavenums(nx, ny, dx, dy, kx, ky, k2)
    integer,        intent(in)  :: nx, ny
    real(c_double), intent(in)  :: dx, dy
    real(c_double), intent(out) :: kx(0:nx-1, 0:ny-1)
    real(c_double), intent(out) :: ky(0:nx-1, 0:ny-1)
    real(c_double), intent(out) :: k2(0:nx-1, 0:ny-1)

    integer        :: i, j
    real(c_double) :: kx1(0:nx-1), ky1(0:ny-1)
    real(c_double) :: base_x, base_y

    base_x = KH_TWO_PI / (real(nx, c_double) * dx)
    base_y = KH_TWO_PI / (real(ny, c_double) * dy)

    ! 1D x-wavenumbers (fftfreq convention)
    do i = 0, nx/2
      kx1(i) = base_x * real(i, c_double)
    end do
    do i = nx/2 + 1, nx-1
      kx1(i) = base_x * real(i - nx, c_double)
    end do

    ! 1D y-wavenumbers
    do j = 0, ny/2
      ky1(j) = base_y * real(j, c_double)
    end do
    do j = ny/2 + 1, ny-1
      ky1(j) = base_y * real(j - ny, c_double)
    end do

    ! Expand to 2D
    do j = 0, ny-1
      do i = 0, nx-1
        kx(i, j) = kx1(i)
        ky(i, j) = ky1(j)
        k2(i, j) = kx1(i)**2 + ky1(j)**2
      end do
    end do

    ! Zero mode: set k2(0,0) = 1 to avoid division by zero in Poisson solve.
    ! The solution at (0,0) is always forced to zero afterward (mean psi = 0).
    k2(0, 0) = 1.0_c_double

  end subroutine kh_grid_wavenums


  ! ── 2/3-rule de-aliasing mask ────────────────────────────────────────────────
  !
  ! Returns a logical mask over (nx, ny) wavenumber space.
  ! mask(i,j) = .true.  → mode is kept
  ! mask(i,j) = .false. → mode is zeroed (aliased)
  !
  ! Criterion (Orszag 1971):
  !   |kx(i,j)| > (2/3) * kx_max  OR  |ky(i,j)| > (2/3) * ky_max
  !   → zero; otherwise keep.
  !
  ! kx_max = π/dx = (nx/2) * 2π/(nx*dx)
  ! ky_max = π/dy = (ny/2) * 2π/(ny*dy)
  !
  ! Arguments:
  !   nx, ny      [in]  grid resolution
  !   kx, ky      [in]  wavenumber arrays (from kh_grid_wavenums)
  !   dx, dy      [in]  grid spacing
  !   mask        [out] logical array; .true. = keep mode

  subroutine kh_grid_dealias_mask(nx, ny, kx, ky, dx, dy, mask)
    integer,        intent(in)  :: nx, ny
    real(c_double), intent(in)  :: kx(0:nx-1, 0:ny-1)
    real(c_double), intent(in)  :: ky(0:nx-1, 0:ny-1)
    real(c_double), intent(in)  :: dx, dy
    logical,        intent(out) :: mask(0:nx-1, 0:ny-1)

    integer        :: i, j
    real(c_double) :: kx_max, ky_max, kx_cut, ky_cut

    kx_max = KH_PI / dx
    ky_max = KH_PI / dy
    kx_cut = KH_DEALIAS_FACTOR * kx_max
    ky_cut = KH_DEALIAS_FACTOR * ky_max

    do j = 0, ny-1
      do i = 0, nx-1
        mask(i, j) = (abs(kx(i,j)) <= kx_cut) .and. (abs(ky(i,j)) <= ky_cut)
      end do
    end do
  end subroutine kh_grid_dealias_mask

end module kh_grid
