! kh_poisson.f90 — KH instability: spectral Poisson solver
!
! REFERENCE WATERMARK — kh-sim Fortran reference, v0.1;
! canonical for: TC-NUM-KH-001..008;
! do not modify without bumping reference version + regenerating reference outputs.
!
! Solves:  ∇²ψ = -ω   →   in spectral space: ψ̂(kx,ky) = ω̂(kx,ky) / k²
!
! Convention:
!   The Poisson equation in vorticity-streamfunction form is:
!     ∇²ψ = -ω   (KH-PHYSICS §2)
!   In Fourier space (-kx²-ky²)ψ̂ = -ω̂, so ψ̂ = ω̂ / k²  (k² = kx²+ky²).
!   Zero mode: ψ̂(0,0) = 0  (mean streamfunction pinned to zero).
!   k2(0,0) is set to 1.0 by kh_grid_wavenums to avoid division by zero;
!   this routine unconditionally zeros psi_hat(0,0) afterward.
!
! Public interface:
!   kh_poisson_solve(omega_hat, k2, nx, ny, psi_hat)
!
! Reference: kh-sim/shared/physics/KH-PHYSICS.md §2, §4
! Numerical methods: _config/PHYSICS-NUMERICAL-METHODS-v0.1.md §2.2

module kh_poisson
  use kh_constants, only: c_double
  implicit none
  private

  public :: kh_poisson_solve

contains

  ! ── Spectral Poisson solve ───────────────────────────────────────────────────
  !
  ! Arguments:
  !   omega_hat(0:nx-1, 0:ny-1) [in]  Fourier transform of vorticity ω
  !   k2(0:nx-1, 0:ny-1)        [in]  kx²+ky² wavenumber-squared array
  !                                   (k2(0,0) must be 1.0 per kh_grid convention)
  !   nx, ny                    [in]  grid dimensions
  !   psi_hat(0:nx-1, 0:ny-1)   [out] Fourier transform of streamfunction ψ
  !                                   psi_hat = omega_hat / k2; psi_hat(0,0) = 0

  subroutine kh_poisson_solve(omega_hat, k2, nx, ny, psi_hat)
    integer,        intent(in)  :: nx, ny
    complex(8),     intent(in)  :: omega_hat(0:nx-1, 0:ny-1)
    real(c_double), intent(in)  :: k2(0:nx-1, 0:ny-1)
    complex(8),     intent(out) :: psi_hat(0:nx-1, 0:ny-1)

    integer :: i, j

    do j = 0, ny-1
      do i = 0, nx-1
        psi_hat(i, j) = omega_hat(i, j) / k2(i, j)
      end do
    end do

    ! Enforce mean streamfunction = 0 (zero mode)
    psi_hat(0, 0) = cmplx(0.0_8, 0.0_8, kind=8)

  end subroutine kh_poisson_solve

end module kh_poisson
