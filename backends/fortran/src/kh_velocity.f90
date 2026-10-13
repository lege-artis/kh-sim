! kh_velocity.f90 — KH instability: spectral velocity recovery
!
! REFERENCE WATERMARK — kh-sim Fortran reference, v0.1;
! canonical for: TC-NUM-KH-001..008;
! do not modify without bumping reference version + regenerating reference outputs.
!
! Recovers the physical-space velocity (u, v) from the spectral streamfunction
! ψ̂ using the exact spectral derivative relations:
!
!   û(kx,ky) =  i·ky · ψ̂(kx,ky)     →  u(x,y) = IFFT2(û)
!   v̂(kx,ky) = -i·kx · ψ̂(kx,ky)    →  v(x,y) = IFFT2(v̂)
!
! The continuity equation ∂u/∂x + ∂v/∂y = 0 is satisfied identically by
! construction (spectral accuracy). See KH-PHYSICS §2.
!
! Public interface:
!   kh_velocity_from_psi(psi_hat, kx, ky, nx, ny, u, v)
!
! Reference: kh-sim/shared/physics/KH-PHYSICS.md §2
! Numerical methods: _config/PHYSICS-NUMERICAL-METHODS-v0.1.md §2.2

module kh_velocity
  use kh_constants, only: c_double
  use kh_fft,       only: kh_fft_inverse_2d
  implicit none
  private

  public :: kh_velocity_from_psi

  ! Imaginary unit (kind=8 complex)
  complex(8), parameter :: I_UNIT = cmplx(0.0_8, 1.0_8, kind=8)

contains

  ! ── Velocity recovery from streamfunction ────────────────────────────────────
  !
  ! Arguments:
  !   psi_hat(0:nx-1, 0:ny-1) [in]  spectral streamfunction ψ̂ (from kh_poisson)
  !   kx(0:nx-1, 0:ny-1)      [in]  x-wavenumber array (from kh_grid_wavenums)
  !   ky(0:nx-1, 0:ny-1)      [in]  y-wavenumber array
  !   nx, ny                  [in]  grid dimensions
  !   u(0:nx-1, 0:ny-1)       [out] x-velocity field (physical space)
  !   v(0:nx-1, 0:ny-1)       [out] y-velocity field (physical space)

  subroutine kh_velocity_from_psi(psi_hat, kx, ky, nx, ny, u, v)
    integer,        intent(in)  :: nx, ny
    complex(8),     intent(in)  :: psi_hat(0:nx-1, 0:ny-1)
    real(c_double), intent(in)  :: kx(0:nx-1, 0:ny-1)
    real(c_double), intent(in)  :: ky(0:nx-1, 0:ny-1)
    real(c_double), intent(out) :: u(0:nx-1, 0:ny-1)
    real(c_double), intent(out) :: v(0:nx-1, 0:ny-1)

    complex(8) :: u_hat(0:nx-1, 0:ny-1)
    complex(8) :: v_hat(0:nx-1, 0:ny-1)
    integer :: i, j

    ! Spectral velocity:  û =  i·ky·ψ̂,   v̂ = -i·kx·ψ̂
    do j = 0, ny-1
      do i = 0, nx-1
        u_hat(i, j) =  I_UNIT * ky(i, j) * psi_hat(i, j)
        v_hat(i, j) = -I_UNIT * kx(i, j) * psi_hat(i, j)
      end do
    end do

    ! Inverse FFT to physical space
    call kh_fft_inverse_2d(u_hat, nx, ny)
    call kh_fft_inverse_2d(v_hat, nx, ny)

    ! Extract real part (imaginary part is numerical noise ~1e-15)
    u = real(u_hat, kind=8)
    v = real(v_hat, kind=8)

  end subroutine kh_velocity_from_psi

end module kh_velocity
