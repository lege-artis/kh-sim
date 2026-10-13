! kh_initial_conditions.f90 — KH instability: shared initial-condition routines
!
! OQ-MAP-05 resolution (2026-05-13): replaces duplicate IC code that existed in
! both kh_physics.f90 (HTTP-shim path) and kh_solver.f90 (ETDRK4 reference path).
! Both paths now call kh_ic_tanh_shear with explicit parameter passing.
!
! IC formula (KH-GE-4, KH-PHYSICS.md §3):
!   ω₀(x,y) = (U₀/δ)·sech²((y − Ly/2)/δ)
!            + amp·(2π·mode/Lx)·cos(2π·mode·x/Lx)
!
! Reference: kh-sim/shared/physics/KH-PHYSICS.md §3
! Canonical: kh-sim/shared/physics/kh_physics.py::initial_conditions

module kh_initial_conditions
  use kh_constants, only: c_double, KH_TWO_PI
  implicit none
  private

  public :: kh_ic_tanh_shear

contains

  ! ── Tanh shear-layer initial vorticity ─────────────────────────────────────
  !
  ! Computes the analytic KH initial vorticity field.
  !
  ! Arguments:
  !   nx, ny  [in]  grid resolution
  !   lx, ly  [in]  domain size
  !   u0      [in]  shear velocity magnitude
  !   delta   [in]  shear layer thickness (= KH_DELTA_FACTOR * ly)
  !   amp     [in]  perturbation amplitude
  !   mode    [in]  perturbation wavenumber mode
  !   omega   [out] vorticity field omega(0:nx-1, 0:ny-1)

  subroutine kh_ic_tanh_shear(nx, ny, lx, ly, u0, delta, amp, mode, omega)
    integer,        intent(in)  :: nx, ny, mode
    real(c_double), intent(in)  :: lx, ly, u0, delta, amp
    real(c_double), intent(out) :: omega(0:nx-1, 0:ny-1)

    integer        :: i, j
    real(c_double) :: dx, dy, x, y, z, pert_x

    dx = lx / real(nx, c_double)
    dy = ly / real(ny, c_double)

    do i = 0, nx-1
      x = real(i, c_double) * dx
      pert_x = amp * (KH_TWO_PI * real(mode, c_double) / lx) &
                   * cos(KH_TWO_PI * real(mode, c_double) * x / lx)
      do j = 0, ny-1
        y = real(j, c_double) * dy
        z = (y - ly * 0.5_c_double) / delta
        omega(i, j) = u0 / delta / cosh(z)**2 + pert_x
      end do
    end do
  end subroutine kh_ic_tanh_shear

end module kh_initial_conditions
