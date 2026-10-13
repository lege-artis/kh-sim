! kh_nonlinear.f90 — KH instability: nonlinear term assembly + 2/3-rule de-aliasing
!
! REFERENCE WATERMARK — kh-sim Fortran reference, v0.1;
! canonical for: TC-NUM-KH-001..008;
! do not modify without bumping reference version + regenerating reference outputs.
!
! Computes the nonlinear advection term N = -(u·∂ω/∂x + v·∂ω/∂y) in spectral
! space, with Orszag 2/3-rule de-aliasing applied after the quadratic product.
!
! Algorithm (per PHYSICS-NUMERICAL-METHODS §2.2 step 4):
!   1. Forward FFT: ω̂ = FFT2(ω)
!   2. Spectral gradients: ∂̂ω/∂x = i·kx·ω̂,  ∂̂ω/∂y = i·ky·ω̂
!   3. IFFT → dox = ∂ω/∂x,  doy = ∂ω/∂y  (physical space)
!   4. Velocity already available from kh_velocity: u, v
!   5. Nonlinear product in physical space: N = -(u·dox + v·doy)
!   6. Forward FFT: N̂ = FFT2(N)
!   7. De-alias N̂: zero modes outside 2/3 cutoff (Orszag 1971)
!
! Public interface:
!   kh_nonlinear_rhs(omega, u, v, kx, ky, nx, ny, N_hat)
!   kh_dealias(a_hat, mask, nx, ny)
!
! Reference: kh-sim/shared/physics/KH-PHYSICS.md §4; Orszag 1971
! Numerical methods: _config/PHYSICS-NUMERICAL-METHODS-v0.1.md §2.1 (2/3 rule)

module kh_nonlinear
  use kh_constants, only: c_double
  use kh_fft,       only: kh_fft_forward_2d, kh_fft_inverse_2d
  implicit none
  private

  public :: kh_nonlinear_rhs
  public :: kh_dealias

  complex(8), parameter :: I_UNIT = cmplx(0.0_8, 1.0_8, kind=8)

contains

  ! ── Nonlinear RHS in spectral space ─────────────────────────────────────────
  !
  ! Computes N̂ = FFT2(-(u·∂ω/∂x + v·∂ω/∂y)) with 2/3-rule de-aliasing.
  !
  ! Arguments:
  !   omega(0:nx-1, 0:ny-1)  [in]  vorticity field (physical space)
  !   u(0:nx-1, 0:ny-1)      [in]  x-velocity (physical space, from kh_velocity)
  !   v(0:nx-1, 0:ny-1)      [in]  y-velocity (physical space)
  !   kx(0:nx-1, 0:ny-1)     [in]  x-wavenumber array
  !   ky(0:nx-1, 0:ny-1)     [in]  y-wavenumber array
  !   mask(0:nx-1, 0:ny-1)   [in]  2/3-rule de-aliasing mask (from kh_grid_dealias_mask)
  !   nx, ny                 [in]  grid dimensions
  !   N_hat(0:nx-1, 0:ny-1)  [out] de-aliased spectral nonlinear term

  subroutine kh_nonlinear_rhs(omega, u, v, kx, ky, mask, nx, ny, N_hat)
    integer,        intent(in)  :: nx, ny
    real(c_double), intent(in)  :: omega(0:nx-1, 0:ny-1)
    real(c_double), intent(in)  :: u(0:nx-1, 0:ny-1)
    real(c_double), intent(in)  :: v(0:nx-1, 0:ny-1)
    real(c_double), intent(in)  :: kx(0:nx-1, 0:ny-1)
    real(c_double), intent(in)  :: ky(0:nx-1, 0:ny-1)
    logical,        intent(in)  :: mask(0:nx-1, 0:ny-1)
    complex(8),     intent(out) :: N_hat(0:nx-1, 0:ny-1)

    complex(8) :: oh(0:nx-1, 0:ny-1)   ! spectral vorticity
    complex(8) :: dx_s(0:nx-1, 0:ny-1) ! spectral ∂ω/∂x
    complex(8) :: dy_s(0:nx-1, 0:ny-1) ! spectral ∂ω/∂y
    real(c_double) :: dox(0:nx-1, 0:ny-1) ! physical ∂ω/∂x
    real(c_double) :: doy(0:nx-1, 0:ny-1) ! physical ∂ω/∂y
    real(c_double) :: nonlin(0:nx-1, 0:ny-1) ! physical N = -(u dox + v doy)
    integer :: i, j

    ! Step 1: Forward FFT of ω
    oh = cmplx(omega, 0.0_8, kind=8)
    call kh_fft_forward_2d(oh, nx, ny)

    ! Step 2: Spectral gradient operators
    do j = 0, ny-1
      do i = 0, nx-1
        dx_s(i, j) = I_UNIT * kx(i, j) * oh(i, j)
        dy_s(i, j) = I_UNIT * ky(i, j) * oh(i, j)
      end do
    end do

    ! Step 3: IFFT → physical gradients
    call kh_fft_inverse_2d(dx_s, nx, ny);  dox = real(dx_s, kind=8)
    call kh_fft_inverse_2d(dy_s, nx, ny);  doy = real(dy_s, kind=8)

    ! Step 4: Physical-space product N = -(u·∂ω/∂x + v·∂ω/∂y)
    nonlin = -(u * dox + v * doy)

    ! Step 5: Forward FFT of nonlinear term
    N_hat = cmplx(nonlin, 0.0_8, kind=8)
    call kh_fft_forward_2d(N_hat, nx, ny)

    ! Step 6: Apply 2/3-rule de-aliasing mask
    call kh_dealias(N_hat, mask, nx, ny)

  end subroutine kh_nonlinear_rhs


  ! ── 2/3-rule de-aliasing ────────────────────────────────────────────────────
  !
  ! Zeros spectral modes outside the de-aliasing cutoff in-place.
  ! mask(i,j) = .true.  → mode kept
  ! mask(i,j) = .false. → mode set to zero
  !
  ! Arguments:
  !   a_hat(0:nx-1, 0:ny-1) [inout]  spectral array to de-alias
  !   mask(0:nx-1, 0:ny-1)  [in]     Boolean mask (from kh_grid_dealias_mask)
  !   nx, ny                [in]     grid dimensions

  subroutine kh_dealias(a_hat, mask, nx, ny)
    integer,    intent(in)    :: nx, ny
    complex(8), intent(inout) :: a_hat(0:nx-1, 0:ny-1)
    logical,    intent(in)    :: mask(0:nx-1, 0:ny-1)
    integer :: i, j

    do j = 0, ny-1
      do i = 0, nx-1
        if (.not. mask(i, j)) a_hat(i, j) = cmplx(0.0_8, 0.0_8, kind=8)
      end do
    end do
  end subroutine kh_dealias

end module kh_nonlinear
