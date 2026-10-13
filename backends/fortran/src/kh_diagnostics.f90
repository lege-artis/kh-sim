! kh_diagnostics.f90 — KH instability: run-time scalar diagnostics
!
! REFERENCE WATERMARK — kh-sim Fortran reference, v0.1;
! canonical for: TC-NUM-KH-001..008;
! do not modify without bumping reference version + regenerating reference outputs.
!
! Computes scalar diagnostics from physical-space vorticity + velocity fields:
!
!   KE            = (1/2) · ⟨u² + v²⟩              (domain-averaged kinetic energy)
!   enstrophy     = (1/2) · ⟨ω²⟩                    (domain-averaged enstrophy)
!   max_vorticity = max_{i,j} |ω(i,j)|              (peak vorticity magnitude)
!   divergence_rms = rms(∂u/∂x + ∂v/∂y)            (spectral incompressibility check)
!
! Divergence is computed spectrally for maximum accuracy:
!   û  = FFT2(u),   v̂  = FFT2(v)
!   div_hat(k,l) = i·kx(k,l)·û(k,l) + i·ky(k,l)·v̂(k,l)   [should be ≈ 0]
!   div_rms = (1/N) · √(Σ|div_hat|²)                        [Parseval, N=nx·ny]
!
! For a divergence-free velocity field derived from a streamfunction (u=∂ψ/∂y,
! v=-∂ψ/∂x) the analytical div is exactly zero; div_rms measures numerical
! residual, expected ≤ KH_TOL_DIVERGENCE (1e-10).
!
! Public interface:
!   kh_diagnostics_compute(omega, u, v, kx, ky, nx, ny, ke, enstrophy, max_vort, div_rms)
!
! References:
!   KH-PHYSICS.md §6 (diagnostics definitions)
!   PHYSICS-NUMERICAL-METHODS-v0.1.md §2.5

module kh_diagnostics
  use kh_constants, only: c_double
  use kh_fft,       only: kh_fft_forward_2d
  implicit none
  private

  public :: kh_diagnostics_compute

contains

  ! ── Compute all scalar diagnostics from physical fields ──────────────────────
  !
  ! Arguments:
  !   omega(0:nx-1, 0:ny-1) [in]  vorticity field (physical space)
  !   u(0:nx-1, 0:ny-1)     [in]  x-velocity field (physical space)
  !   v(0:nx-1, 0:ny-1)     [in]  y-velocity field (physical space)
  !   kx(0:nx-1, 0:ny-1)    [in]  x angular wavenumber array
  !   ky(0:nx-1, 0:ny-1)    [in]  y angular wavenumber array
  !   nx, ny                [in]  grid dimensions
  !   ke                    [out] kinetic energy (domain-averaged)
  !   enstrophy             [out] enstrophy (domain-averaged)
  !   max_vort              [out] peak vorticity magnitude
  !   div_rms               [out] divergence rms (spectral, should be ≈ 0)

  subroutine kh_diagnostics_compute(omega, u, v, kx, ky, nx, ny, &
                                     ke, enstrophy, max_vort, div_rms)
    integer,        intent(in)  :: nx, ny
    real(c_double), intent(in)  :: omega(0:nx-1, 0:ny-1)
    real(c_double), intent(in)  :: u(0:nx-1, 0:ny-1)
    real(c_double), intent(in)  :: v(0:nx-1, 0:ny-1)
    real(c_double), intent(in)  :: kx(0:nx-1, 0:ny-1)
    real(c_double), intent(in)  :: ky(0:nx-1, 0:ny-1)
    real(c_double), intent(out) :: ke
    real(c_double), intent(out) :: enstrophy
    real(c_double), intent(out) :: max_vort
    real(c_double), intent(out) :: div_rms

    complex(8) :: u_hat(0:nx-1, 0:ny-1)
    complex(8) :: v_hat(0:nx-1, 0:ny-1)
    real(c_double) :: div_re, div_im, div_sum
    integer :: i, j, ntot

    ntot = nx * ny

    ! ── Kinetic energy: KE = (1/2) · Σ(u²+v²) / N ───────────────────────────────

    ke = 0.0_c_double
    do j = 0, ny-1
      do i = 0, nx-1
        ke = ke + u(i,j)*u(i,j) + v(i,j)*v(i,j)
      end do
    end do
    ke = 0.5_c_double * ke / real(ntot, c_double)

    ! ── Enstrophy: Ω = (1/2) · Σω² / N ──────────────────────────────────────────

    enstrophy = 0.0_c_double
    do j = 0, ny-1
      do i = 0, nx-1
        enstrophy = enstrophy + omega(i,j)*omega(i,j)
      end do
    end do
    enstrophy = 0.5_c_double * enstrophy / real(ntot, c_double)

    ! ── Peak vorticity magnitude ──────────────────────────────────────────────────

    max_vort = 0.0_c_double
    do j = 0, ny-1
      do i = 0, nx-1
        if (abs(omega(i,j)) > max_vort) max_vort = abs(omega(i,j))
      end do
    end do

    ! ── Divergence rms (spectral) ─────────────────────────────────────────────────
    !
    ! û  = FFT2(u),  v̂  = FFT2(v)  [unnormalised forward transform]
    ! div_hat = i·kx·û + i·ky·v̂
    !         → Re(div_hat) = -(kx·Im(û) + ky·Im(v̂))
    !            Im(div_hat) =  kx·Re(û) + ky·Re(v̂)
    !
    ! Parseval [our convention: Σ|f|² = (1/N)·Σ|f_hat|²]:
    !   mean(div²) = Σ|div_hat|² / N²
    !   div_rms    = √(Σ|div_hat|²) / N

    u_hat = cmplx(u, 0.0_8, kind=8)
    v_hat = cmplx(v, 0.0_8, kind=8)
    call kh_fft_forward_2d(u_hat, nx, ny)
    call kh_fft_forward_2d(v_hat, nx, ny)

    div_sum = 0.0_c_double
    do j = 0, ny-1
      do i = 0, nx-1
        div_re = -(kx(i,j)*aimag(u_hat(i,j)) + ky(i,j)*aimag(v_hat(i,j)))
        div_im =   kx(i,j)*real(u_hat(i,j), kind=8) + ky(i,j)*real(v_hat(i,j), kind=8)
        div_sum = div_sum + div_re*div_re + div_im*div_im
      end do
    end do
    div_rms = sqrt(div_sum) / real(ntot, c_double)

  end subroutine kh_diagnostics_compute

end module kh_diagnostics
