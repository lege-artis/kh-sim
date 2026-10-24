! kh_solver.f90 — KH instability: top-level ETDRK4 time-integration driver
!
! REFERENCE WATERMARK — Fortran ETDRK4 ALTERNATE-SCHEME reference, v0.2;
! certifies: TC-NUM-KH-003, TC-NUM-KH-004 (ETDRK4-specific accuracy + order).
!
! Per phase 5a Option D (2026-05-13), the canonical kh-sim time integrator
! is RK4 (KH-NM-3a). ETDRK4 (KH-NM-3b) is retained as an alternate scheme
! for stiff-problem convergence-order testing only.
!
! Canonical-scheme reference is the Python ref:
!   shared/physics/kh_physics.py
!
! Canonical-scheme Fortran HTTP-shim path:
!   backends/fortran/src/kh_physics.f90  (kh_simulate_c / kh_simulate_rk4)
!
! TC-NUM-KH-001/002/005/006/007/008 certify the canonical RK4 scheme and
! target the kh_physics.f90 path (see development record, resolutions 2026-05-13).
!
! Do not modify without bumping ETDRK4-reference version + regenerating
! TC-NUM-KH-003/004 reference outputs.
!
! Provides a single high-level entry point kh_solver_run that:
!   1. Builds grid, wavenumber arrays, and de-aliasing mask
!   2. Initialises the KH shear-layer initial condition
!   3. Precomputes ETDRK4 coefficient arrays (L = -ν·k²)
!   4. Runs the time loop, calling kh_etdrk4_step each step
!   5. Monitors CFL number after every step (peak value returned to caller)
!   6. Computes diagnostics and writes JSON snapshots at out_interval steps
!   7. Returns final scalar diagnostics to caller
!
! CFL number definition (per-step):
!   CFL = dt · (max|u| / dx + max|v| / dy)
! The caller should verify CFL_peak ≤ KH_CFL_FACTOR (0.4).
! The solver does NOT abort on CFL violation in v0.1 — it reports via
! cfl_peak so the caller can detect and handle the condition.
!
! Initial condition (matches kh_physics.f90 convention):
!   ω₀(x,y) = (U₀/δ) · sech²((y − Ly/2) / δ)
!            + amp · (2π·mode/Lx) · cos(2π·mode·x / Lx)
!   where δ = KH_DELTA_FACTOR · Ly
!
! Public interface:
!   kh_solver_run(nx, ny, lx, ly, nu, dt, nsteps, amp, mode,
!                 out_interval, out_unit,
!                 ke_out, enstrophy_out, max_vort_out, div_rms_out, cfl_peak)
!
! References:
!   KH-PHYSICS.md §3–4 (IC, flow parameters)
!   NM development spec §2.4 (solver architecture)

module kh_solver
  use kh_constants
  use kh_grid,               only: kh_grid_make, kh_grid_wavenums, kh_grid_dealias_mask
  use kh_fft,                only: kh_fft_forward_2d, kh_fft_inverse_2d
  use kh_poisson,            only: kh_poisson_solve
  use kh_velocity,           only: kh_velocity_from_psi
  use kh_etdrk4,             only: kh_etdrk4_precompute, kh_etdrk4_step
  use kh_diagnostics,        only: kh_diagnostics_compute
  use kh_io,                 only: kh_io_write_snapshot
  use kh_initial_conditions, only: kh_ic_tanh_shear
  implicit none
  private

  public :: kh_solver_run

contains

  ! ── Top-level time integration driver ────────────────────────────────────────
  !
  ! Arguments:
  !   nx, ny         [in]  grid dimensions (must be powers of 2)
  !   lx, ly         [in]  domain extents
  !   nu             [in]  kinematic viscosity (0.0 for inviscid)
  !   dt             [in]  time step
  !   nsteps         [in]  total number of steps
  !   amp            [in]  perturbation amplitude
  !   mode           [in]  perturbation wavenumber mode
  !   out_interval   [in]  steps between JSON snapshot writes (0 = no output)
  !   out_unit       [in]  file unit from kh_io_open_output (-1 = no output)
  !   ke_out         [out] final kinetic energy
  !   enstrophy_out  [out] final enstrophy
  !   max_vort_out   [out] final peak vorticity
  !   div_rms_out    [out] final divergence rms
  !   cfl_peak       [out] maximum CFL number observed across all steps

  subroutine kh_solver_run(nx, ny, lx, ly, nu, dt, nsteps, amp, mode, &
                             out_interval, out_unit, &
                             ke_out, enstrophy_out, max_vort_out, div_rms_out, cfl_peak)
    integer,        intent(in)  :: nx, ny, nsteps, mode, out_interval, out_unit
    real(c_double), intent(in)  :: lx, ly, nu, dt, amp
    real(c_double), intent(out) :: ke_out, enstrophy_out, max_vort_out, div_rms_out
    real(c_double), intent(out) :: cfl_peak

    ! ── Grid arrays (allocatable — nx, ny are runtime parameters) ───────────────

    real(c_double), allocatable :: x(:), y(:)
    real(c_double), allocatable :: kx(:,:), ky(:,:), k2(:,:)
    logical,        allocatable :: mask(:,:)

    ! ── Field arrays ─────────────────────────────────────────────────────────────

    real(c_double), allocatable :: omega_phys(:,:)
    real(c_double), allocatable :: u(:,:), v(:,:)
    complex(8),     allocatable :: omega_hat(:,:), psi_hat(:,:), tmp(:,:)

    ! ── ETDRK4 coefficient arrays ─────────────────────────────────────────────────

    real(c_double), allocatable :: L_op(:,:)
    real(c_double), allocatable :: E_arr(:,:), E2_arr(:,:)
    real(c_double), allocatable :: phi1_arr(:,:), phi1h_arr(:,:), phi2_arr(:,:)

    ! ── Scalar temporaries ────────────────────────────────────────────────────────

    real(c_double) :: dx, dy, time, cfl, delta
    real(c_double) :: ke, enstrophy, max_vort, div_rms
    integer        :: s

    ! ── Allocate all arrays ───────────────────────────────────────────────────────

    allocate(x(0:nx-1), y(0:ny-1))
    allocate(kx(0:nx-1, 0:ny-1), ky(0:nx-1, 0:ny-1), k2(0:nx-1, 0:ny-1))
    allocate(mask(0:nx-1, 0:ny-1))
    allocate(omega_phys(0:nx-1, 0:ny-1))
    allocate(u(0:nx-1, 0:ny-1), v(0:nx-1, 0:ny-1))
    allocate(omega_hat(0:nx-1, 0:ny-1))
    allocate(psi_hat(0:nx-1, 0:ny-1))
    allocate(tmp(0:nx-1, 0:ny-1))
    allocate(L_op(0:nx-1, 0:ny-1))
    allocate(E_arr(0:nx-1, 0:ny-1), E2_arr(0:nx-1, 0:ny-1))
    allocate(phi1_arr(0:nx-1, 0:ny-1), phi1h_arr(0:nx-1, 0:ny-1))
    allocate(phi2_arr(0:nx-1, 0:ny-1))

    ! ── 1. Build grid ─────────────────────────────────────────────────────────────

    call kh_grid_make(nx, ny, lx, ly, x, y, dx, dy)
    call kh_grid_wavenums(nx, ny, dx, dy, kx, ky, k2)
    call kh_grid_dealias_mask(nx, ny, kx, ky, dx, dy, mask)

    ! ── 2. Initial condition — unified via kh_initial_conditions (OQ-MAP-05) ──────
    !   Replaces inline IC that previously used KH_U0_DEFAULT constant directly.
    !   Now passes u0 = KH_U0_DEFAULT explicitly for parameter-transparent calling.

    delta = KH_DELTA_FACTOR * ly
    call kh_ic_tanh_shear(nx, ny, lx, ly, KH_U0_DEFAULT, delta, amp, mode, omega_phys)

    ! ── 3. Forward FFT → spectral vorticity ──────────────────────────────────────

    omega_hat = cmplx(omega_phys, 0.0_8, kind=8)
    call kh_fft_forward_2d(omega_hat, nx, ny)

    ! ── 4. ETDRK4 precompute: L = −ν·k² ─────────────────────────────────────────

    L_op = -nu * k2
    call kh_etdrk4_precompute(L_op, dt, nx, ny, E_arr, E2_arr, &
                               phi1_arr, phi1h_arr, phi2_arr)

    ! ── 5. Time loop ──────────────────────────────────────────────────────────────

    time     = 0.0_c_double
    cfl_peak = 0.0_c_double

    do s = 1, nsteps

      ! Advance one ETDRK4 step
      call kh_etdrk4_step(omega_hat, E_arr, E2_arr, phi1_arr, phi1h_arr, phi2_arr, &
                           dt, kx, ky, k2, mask, nx, ny)
      time = time + dt

      ! Recover velocity for CFL monitoring (exact, spectral)
      call kh_poisson_solve(omega_hat, k2, nx, ny, psi_hat)
      call kh_velocity_from_psi(psi_hat, kx, ky, nx, ny, u, v)

      ! CFL = dt · (max|u|/dx + max|v|/dy)
      cfl = dt * (maxval(abs(u)) / dx + maxval(abs(v)) / dy)
      if (cfl > cfl_peak) cfl_peak = cfl

      ! Diagnostics + JSON snapshot at out_interval
      if (out_interval > 0 .and. mod(s, out_interval) == 0) then
        tmp = omega_hat
        call kh_fft_inverse_2d(tmp, nx, ny)
        omega_phys = real(tmp, kind=8)
        call kh_diagnostics_compute(omega_phys, u, v, kx, ky, nx, ny, &
                                     ke, enstrophy, max_vort, div_rms)
        call kh_io_write_snapshot(out_unit, s, time, ke, enstrophy, max_vort, div_rms)
      end if

    end do

    ! ── 6. Final diagnostics ──────────────────────────────────────────────────────

    tmp = omega_hat
    call kh_fft_inverse_2d(tmp, nx, ny)
    omega_phys = real(tmp, kind=8)

    ! Velocity already current from last CFL step; reuse u, v
    call kh_diagnostics_compute(omega_phys, u, v, kx, ky, nx, ny, &
                                 ke_out, enstrophy_out, max_vort_out, div_rms_out)

    ! ── 7. Deallocate ─────────────────────────────────────────────────────────────

    deallocate(x, y, kx, ky, k2, mask)
    deallocate(omega_phys, u, v, omega_hat, psi_hat, tmp)
    deallocate(L_op, E_arr, E2_arr, phi1_arr, phi1h_arr, phi2_arr)

  end subroutine kh_solver_run

end module kh_solver
