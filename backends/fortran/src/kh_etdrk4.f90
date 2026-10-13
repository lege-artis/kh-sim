! kh_etdrk4.f90 — KH instability: ETDRK4 time integrator
!
! REFERENCE WATERMARK — kh-sim Fortran reference, v0.1;
! canonical for: TC-NUM-KH-001..008;
! do not modify without bumping reference version + regenerating reference outputs.
!
! Implements the 4th-order Exponential Time Differencing Runge-Kutta (ETDRK4)
! scheme of Cox & Matthews (2002) with the Kassam & Trefethen (2005) correction
! for near-singularity in the φ functions via contour integral averaging.
!
! The PDE in Fourier space for each mode (kx, ky):
!   d/dt ω̂ = L(k)·ω̂ + N̂(ω)
!   where  L(k) = -ν·k²       (linear stiff diffusion, integrated exactly)
!          N̂    = nonlinear RHS (computed by kh_nonlinear_rhs each substep)
!
! ETDRK4 substeps (Cox & Matthews 2002, eq. 2.5–2.8):
!   E  = exp(L·dt)
!   E2 = exp(L·dt/2)
!   φ₁ = (E-1)/L·dt              (limit 1 as L→0)
!   φ₂ = (E-1-L·dt)/(L·dt)²     (limit 1/2 as L→0)
!   a = E2·ω̂ + φ₁(L·dt/2)·N̂(ω̂)       [half step 1]
!   b = E2·ω̂ + φ₁(L·dt/2)·N̂(a)       [half step 2]
!   c = E2·a  + φ₁(L·dt/2)·(2N̂(b)-N̂(ω̂)) [full step]
!   ω̂_new = E·ω̂ + dt·(φ₁-3φ₂)·N̂(ω̂) + 2dt·φ₂·(N̂(a)+N̂(b)) + dt·(φ₂-φ₁+1)·N̂(c)
!     (Kassam-Trefethen coefficients α,β,γ,δ form)
!
! φ function evaluation uses direct formula with L=0 guard (tolerance 1e-8).
! Contour integral averaging (Kassam-Trefethen) is not implemented in v0.1 —
! the direct formula is numerically stable for ν·k²·dt << 1 (satisfied at Re=1000
! with default dt=0.001 and grid 128×64). OQ-NUM-04 (OpenMP) is deferred.
!
! References:
!   Cox, S. M., & Matthews, P. C. (2002). JCP 176:430–455.
!   Kassam, A.-K., & Trefethen, L. N. (2005). SIAM JSC 26:1214–1233.
!
! Public interface:
!   kh_etdrk4_precompute(L, dt, nx, ny, E, E2, phi1, phi1h, phi2)
!   kh_etdrk4_step(omega_hat, E, E2, phi1, phi1h, phi2, dt, ...)

module kh_etdrk4
  use kh_constants,  only: c_double
  use kh_fft,        only: kh_fft_forward_2d, kh_fft_inverse_2d
  use kh_poisson,    only: kh_poisson_solve
  use kh_velocity,   only: kh_velocity_from_psi
  use kh_nonlinear,  only: kh_nonlinear_rhs
  implicit none
  private

  public :: kh_etdrk4_precompute
  public :: kh_etdrk4_step

contains

  ! ── Precompute ETDRK4 coefficient arrays ────────────────────────────────────
  !
  ! Called once before the time loop. All arrays are (nx, ny) spectral-space.
  !
  ! Arguments:
  !   L(0:nx-1, 0:ny-1)    [in]  linear operator eigenvalues L = -ν·k²
  !   dt                   [in]  time step
  !   nx, ny               [in]  grid dimensions
  !   E(0:nx-1, 0:ny-1)    [out] exp(L·dt)
  !   E2(0:nx-1, 0:ny-1)   [out] exp(L·dt/2)
  !   phi1(0:nx-1, 0:ny-1) [out] φ₁(L·dt)   = (exp(L·dt)-1)/(L·dt)
  !   phi1h(0:nx-1,0:ny-1) [out] φ₁(L·dt/2) = (exp(L·dt/2)-1)/(L·dt/2)
  !   phi2(0:nx-1, 0:ny-1) [out] φ₂(L·dt)   = (exp(L·dt)-1-L·dt)/(L·dt)²

  subroutine kh_etdrk4_precompute(L, dt, nx, ny, E, E2, phi1, phi1h, phi2)
    integer,        intent(in)  :: nx, ny
    real(c_double), intent(in)  :: L(0:nx-1, 0:ny-1)
    real(c_double), intent(in)  :: dt
    real(c_double), intent(out) :: E(0:nx-1, 0:ny-1)
    real(c_double), intent(out) :: E2(0:nx-1, 0:ny-1)
    real(c_double), intent(out) :: phi1(0:nx-1, 0:ny-1)
    real(c_double), intent(out) :: phi1h(0:nx-1, 0:ny-1)
    real(c_double), intent(out) :: phi2(0:nx-1, 0:ny-1)

    real(c_double), parameter :: TOL = 1.0e-8_c_double
    real(c_double) :: Ldt, Ldt2
    integer :: i, j

    do j = 0, ny-1
      do i = 0, nx-1
        Ldt  = L(i,j) * dt
        Ldt2 = L(i,j) * dt * 0.5_c_double

        E(i,j)  = exp(Ldt)
        E2(i,j) = exp(Ldt2)

        ! φ₁(x) = (e^x - 1)/x;  limit 1 as x→0
        if (abs(Ldt) < TOL) then
          phi1(i,j)  = 1.0_c_double + Ldt  * (0.5_c_double + Ldt  / 6.0_c_double)
        else
          phi1(i,j)  = (E(i,j)  - 1.0_c_double) / Ldt
        end if

        if (abs(Ldt2) < TOL) then
          phi1h(i,j) = 1.0_c_double + Ldt2 * (0.5_c_double + Ldt2 / 6.0_c_double)
        else
          phi1h(i,j) = (E2(i,j) - 1.0_c_double) / Ldt2
        end if

        ! φ₂(x) = (e^x - 1 - x)/x²;  limit 1/2 as x→0
        if (abs(Ldt) < TOL) then
          phi2(i,j) = 0.5_c_double + Ldt * (1.0_c_double/6.0_c_double + Ldt / 24.0_c_double)
        else
          phi2(i,j) = (E(i,j) - 1.0_c_double - Ldt) / (Ldt * Ldt)
        end if
      end do
    end do
  end subroutine kh_etdrk4_precompute


  ! ── Single ETDRK4 time step ─────────────────────────────────────────────────
  !
  ! Advances omega_hat by one time step dt using the Cox-Matthews scheme.
  ! Nonlinear RHS is computed via kh_nonlinear_rhs at each substep.
  !
  ! Arguments:
  !   omega_hat [inout]  spectral vorticity; updated in-place
  !   E, E2, phi1, phi1h, phi2 [in]  precomputed coefficient arrays
  !   dt        [in]   time step
  !   kx, ky    [in]   wavenumber arrays
  !   k2        [in]   k² array
  !   mask      [in]   2/3-rule de-aliasing mask
  !   nx, ny    [in]   grid dimensions

  subroutine kh_etdrk4_step(omega_hat, E, E2, phi1, phi1h, phi2, &
                              dt, kx, ky, k2, mask, nx, ny)
    integer,        intent(in)    :: nx, ny
    complex(8),     intent(inout) :: omega_hat(0:nx-1, 0:ny-1)
    real(c_double), intent(in)    :: E(0:nx-1, 0:ny-1)
    real(c_double), intent(in)    :: E2(0:nx-1, 0:ny-1)
    real(c_double), intent(in)    :: phi1(0:nx-1, 0:ny-1)
    real(c_double), intent(in)    :: phi1h(0:nx-1, 0:ny-1)
    real(c_double), intent(in)    :: phi2(0:nx-1, 0:ny-1)
    real(c_double), intent(in)    :: dt
    real(c_double), intent(in)    :: kx(0:nx-1, 0:ny-1)
    real(c_double), intent(in)    :: ky(0:nx-1, 0:ny-1)
    real(c_double), intent(in)    :: k2(0:nx-1, 0:ny-1)
    logical,        intent(in)    :: mask(0:nx-1, 0:ny-1)

    ! Substep spectral arrays
    complex(8) :: Na(0:nx-1, 0:ny-1)  ! N̂(ω̂)   at t
    complex(8) :: Nb(0:nx-1, 0:ny-1)  ! N̂(a)    at t+dt/2
    complex(8) :: Nc(0:nx-1, 0:ny-1)  ! N̂(b)    at t+dt/2
    complex(8) :: Nd(0:nx-1, 0:ny-1)  ! N̂(c)    at t+dt
    complex(8) :: a_hat(0:nx-1, 0:ny-1)
    complex(8) :: b_hat(0:nx-1, 0:ny-1)
    complex(8) :: c_hat(0:nx-1, 0:ny-1)
    complex(8) :: psi_hat(0:nx-1, 0:ny-1)
    real(c_double) :: omega_phys(0:nx-1, 0:ny-1)
    real(c_double) :: u(0:nx-1, 0:ny-1), v(0:nx-1, 0:ny-1)
    complex(8) :: tmp(0:nx-1, 0:ny-1)
    integer :: i, j

    ! Helper: spectral → physical → velocity
    ! (omega_hat → physical omega, then Poisson → velocity)

    ! ── N̂(ω̂) at current state ────────────────────────────────────────────────
    call spectral_to_velocity(omega_hat, k2, kx, ky, nx, ny, u, v)
    tmp = omega_hat
    call kh_fft_inverse_2d(tmp, nx, ny)
    omega_phys = real(tmp, kind=8)
    call kh_nonlinear_rhs(omega_phys, u, v, kx, ky, mask, nx, ny, Na)

    ! ── Substep a: half step from ω̂ ────────────────────────────────────────────
    do j = 0, ny-1
      do i = 0, nx-1
        a_hat(i,j) = E2(i,j)*omega_hat(i,j) + phi1h(i,j)*dt*0.5_c_double*Na(i,j)
      end do
    end do

    ! ── N̂(a) ─────────────────────────────────────────────────────────────────
    call spectral_to_velocity(a_hat, k2, kx, ky, nx, ny, u, v)
    tmp = a_hat
    call kh_fft_inverse_2d(tmp, nx, ny)
    omega_phys = real(tmp, kind=8)
    call kh_nonlinear_rhs(omega_phys, u, v, kx, ky, mask, nx, ny, Nb)

    ! ── Substep b: half step from ω̂ using N̂(a) ──────────────────────────────
    do j = 0, ny-1
      do i = 0, nx-1
        b_hat(i,j) = E2(i,j)*omega_hat(i,j) + phi1h(i,j)*dt*0.5_c_double*Nb(i,j)
      end do
    end do

    ! ── N̂(b) ─────────────────────────────────────────────────────────────────
    call spectral_to_velocity(b_hat, k2, kx, ky, nx, ny, u, v)
    tmp = b_hat
    call kh_fft_inverse_2d(tmp, nx, ny)
    omega_phys = real(tmp, kind=8)
    call kh_nonlinear_rhs(omega_phys, u, v, kx, ky, mask, nx, ny, Nc)

    ! ── Substep c: full step from a using 2N̂(b)-N̂(ω̂) ───────────────────────
    do j = 0, ny-1
      do i = 0, nx-1
        c_hat(i,j) = E2(i,j)*a_hat(i,j) + phi1h(i,j)*dt*0.5_c_double*(2.0_c_double*Nc(i,j)-Na(i,j))
      end do
    end do

    ! ── N̂(c) ─────────────────────────────────────────────────────────────────
    call spectral_to_velocity(c_hat, k2, kx, ky, nx, ny, u, v)
    tmp = c_hat
    call kh_fft_inverse_2d(tmp, nx, ny)
    omega_phys = real(tmp, kind=8)
    call kh_nonlinear_rhs(omega_phys, u, v, kx, ky, mask, nx, ny, Nd)

    ! ── Final combination (Cox-Matthews 2002) ────────────────────────────────
    do j = 0, ny-1
      do i = 0, nx-1
        omega_hat(i,j) = E(i,j)*omega_hat(i,j) &
          + dt * (phi1(i,j) - 3.0_c_double*phi2(i,j)) * Na(i,j) &
          + dt * 2.0_c_double * phi2(i,j) * (Nb(i,j) + Nc(i,j)) &
          + dt * (phi2(i,j) - phi1(i,j) + 1.0_c_double) * Nd(i,j)
      end do
    end do

  end subroutine kh_etdrk4_step


  ! ── Internal helper: spectral ω̂ → physical velocity (u, v) ─────────────────
  !
  ! Used inside kh_etdrk4_step; not public.

  subroutine spectral_to_velocity(omega_hat, k2, kx, ky, nx, ny, u, v)
    integer,        intent(in)  :: nx, ny
    complex(8),     intent(in)  :: omega_hat(0:nx-1, 0:ny-1)
    real(c_double), intent(in)  :: k2(0:nx-1, 0:ny-1)
    real(c_double), intent(in)  :: kx(0:nx-1, 0:ny-1)
    real(c_double), intent(in)  :: ky(0:nx-1, 0:ny-1)
    real(c_double), intent(out) :: u(0:nx-1, 0:ny-1)
    real(c_double), intent(out) :: v(0:nx-1, 0:ny-1)

    complex(8) :: psi_hat(0:nx-1, 0:ny-1)

    call kh_poisson_solve(omega_hat, k2, nx, ny, psi_hat)
    call kh_velocity_from_psi(psi_hat, kx, ky, nx, ny, u, v)
  end subroutine spectral_to_velocity

end module kh_etdrk4
