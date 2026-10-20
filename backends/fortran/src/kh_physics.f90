! kh_physics.f90 — KH instability physics kernel (Fortran 2008)
!
! Physics: 2D incompressible Navier-Stokes, vorticity-streamfunction form.
! Numerics: pseudo-spectral, RK4 time integration with 2/3-rule dealiasing (KH-NM-4)
!           and CFL stability monitoring (KH-NM-5).
!
! Exported via iso_c_binding for C++ HTTP shim (kh_shim.cpp).
!
! OQ-MAP-03 (2026-05-13): private fft1d/fft2d/make_wavenumbers removed; now uses
!   kh_fft_forward_2d / kh_fft_inverse_2d from kh_fft.f90 and kh_grid_wavenums
!   from kh_grid.f90 — single shared FFT + wavenumber implementation.
! OQ-MAP-05 (2026-05-13): private initial_conditions removed; now uses
!   kh_ic_tanh_shear from kh_initial_conditions.f90.
!
! Reference: kh-sim/shared/physics/KH-PHYSICS.md
! Canonical:  kh-sim/shared/physics/kh_physics.py

module kh_physics
  use iso_c_binding,         only: c_double, c_int, c_ptr, c_f_pointer, c_null_ptr, c_associated
  use kh_constants,          only: KH_DELTA_FACTOR, KH_U0_DEFAULT
  use kh_fft,                only: kh_fft_forward_2d, kh_fft_inverse_2d
  use kh_grid,               only: kh_grid_make, kh_grid_wavenums, kh_grid_dealias_mask
  use kh_nonlinear,          only: kh_dealias
  use kh_initial_conditions, only: kh_ic_tanh_shear
  use kh_diagnostics,        only: kh_diagnostics_compute
  use kh_forcing,            only: build_mode_forcing, build_band_phi0, build_band_omegas, &
                                   MAX_BAND_MODES
  implicit none
  private

  public :: kh_simulate_c, kh_simulate_rk4, rk4_step, rk4_step_fd

  ! Forcing kinds passed by the C shim (kh_simulate_c, argument c_ftype).
  integer, parameter :: FORCING_NONE = 0, FORCING_MODE = 1, FORCING_BAND = 2

  complex(c_double), parameter :: I_UNIT = cmplx(0.0_c_double, 1.0_c_double, kind=c_double)

contains

  ! ── Poisson solve: psi_hat = omega_hat / k^2, zero mode = 0 ─────────────────
  subroutine solve_poisson(omega_hat, kx, ky, nx, ny, psi_hat)
    integer,        intent(in)  :: nx, ny
    complex(8),     intent(in)  :: omega_hat(0:nx-1, 0:ny-1)
    real(c_double), intent(in)  :: kx(0:nx-1, 0:ny-1), ky(0:nx-1, 0:ny-1)
    complex(8),     intent(out) :: psi_hat(0:nx-1, 0:ny-1)
    integer        :: i, j
    real(c_double) :: k2

    do j = 0, ny-1
      do i = 0, nx-1
        if (i == 0 .and. j == 0) then
          psi_hat(0,0) = cmplx(0.0_c_double, 0.0_c_double, kind=8)
        else
          k2 = kx(i,j)**2 + ky(i,j)**2
          psi_hat(i,j) = omega_hat(i,j) / k2
        end if
      end do
    end do
  end subroutine solve_poisson


  ! ── Velocity from streamfunction ──────────────────────────────────────────────
  ! u_hat = i*ky*psi_hat,  v_hat = -i*kx*psi_hat
  subroutine velocity_from_psi(psi_hat, kx, ky, nx, ny, u, v)
    integer,        intent(in)  :: nx, ny
    complex(8),     intent(in)  :: psi_hat(0:nx-1, 0:ny-1)
    real(c_double), intent(in)  :: kx(0:nx-1, 0:ny-1), ky(0:nx-1, 0:ny-1)
    real(c_double), intent(out) :: u(0:nx-1, 0:ny-1), v(0:nx-1, 0:ny-1)

    complex(8) :: u_hat(0:nx-1, 0:ny-1), v_hat(0:nx-1, 0:ny-1)
    integer :: i, j

    do j = 0, ny-1
      do i = 0, nx-1
        u_hat(i,j) =  I_UNIT * ky(i,j) * psi_hat(i,j)
        v_hat(i,j) = -I_UNIT * kx(i,j) * psi_hat(i,j)
      end do
    end do
    call kh_fft_inverse_2d(u_hat, nx, ny)
    call kh_fft_inverse_2d(v_hat, nx, ny)
    u = real(u_hat, kind=c_double)
    v = real(v_hat, kind=c_double)
  end subroutine velocity_from_psi


  ! ── Vorticity RHS with 2/3-rule dealiasing ────────────────────────────────────
  ! f(ω) = -u·∂ω/∂x - v·∂ω/∂y + ν·∇²ω
  ! Nonlinear term dealiased via KH-NM-4 (Orszag 1971).
  subroutine vorticity_rhs(omega, u_in, v_in, nu, kx, ky, mask, nx, ny, rhs)
    integer,        intent(in)  :: nx, ny
    real(c_double), intent(in)  :: omega(0:nx-1, 0:ny-1)
    real(c_double), intent(in)  :: u_in(0:nx-1, 0:ny-1), v_in(0:nx-1, 0:ny-1)
    real(c_double), intent(in)  :: nu
    real(c_double), intent(in)  :: kx(0:nx-1, 0:ny-1), ky(0:nx-1, 0:ny-1)
    logical,        intent(in)  :: mask(0:nx-1, 0:ny-1)
    real(c_double), intent(out) :: rhs(0:nx-1, 0:ny-1)

    complex(8) :: oh(0:nx-1,0:ny-1), dx_s(0:nx-1,0:ny-1)
    complex(8) :: dy_s(0:nx-1,0:ny-1), lap_s(0:nx-1,0:ny-1)
    complex(8) :: nl_hat(0:nx-1,0:ny-1)
    real(c_double) :: dox(0:nx-1,0:ny-1), doy(0:nx-1,0:ny-1)
    real(c_double) :: lap(0:nx-1,0:ny-1), nl_phys(0:nx-1,0:ny-1)
    integer        :: i, j
    real(c_double) :: k2

    oh = cmplx(omega, 0.0_c_double, kind=8)
    call kh_fft_forward_2d(oh, nx, ny)

    do j = 0, ny-1
      do i = 0, nx-1
        dx_s(i,j)  = I_UNIT * kx(i,j) * oh(i,j)
        dy_s(i,j)  = I_UNIT * ky(i,j) * oh(i,j)
        k2         = kx(i,j)**2 + ky(i,j)**2
        lap_s(i,j) = -k2 * oh(i,j)
      end do
    end do

    call kh_fft_inverse_2d(dx_s,  nx, ny);  dox = real(dx_s,  kind=c_double)
    call kh_fft_inverse_2d(dy_s,  nx, ny);  doy = real(dy_s,  kind=c_double)
    call kh_fft_inverse_2d(lap_s, nx, ny);  lap = real(lap_s, kind=c_double)

    ! Nonlinear term with 2/3-rule dealiasing (KH-NM-4, Orszag 1971)
    nl_hat = cmplx(-(u_in*dox + v_in*doy), 0.0_c_double, kind=8)
    call kh_fft_forward_2d(nl_hat, nx, ny)
    call kh_dealias(nl_hat, mask, nx, ny)
    call kh_fft_inverse_2d(nl_hat, nx, ny)
    nl_phys = real(nl_hat, kind=c_double)

    rhs = nl_phys + nu * lap
  end subroutine vorticity_rhs


  ! ── RK4 step ──────────────────────────────────────────────────────────────────
  ! Classical fourth-order Runge-Kutta step for d(omega)/dt = rhs(omega).
  ! Reference: Press, Teukolsky, Vetterling & Flannery (2007), Numerical Recipes,
  !            3rd ed., S17.1, eq. 17.1.3 (classical fourth-order Runge-Kutta).
  ! Convention: explicit; stages at t, t+dt/2, t+dt/2, t+dt; weights 1/6, 1/3, 1/3, 1/6.
  ! Complexity: O(N log N) per step -- four right-hand-side evaluations, each
  !            bound by the 2D FFTs of the pseudo-spectral derivatives.
  subroutine rk4_step(omega, nu, kx, ky, mask, nx, ny, dt, omega_out)
    integer,        intent(in)  :: nx, ny
    real(c_double), intent(in)  :: omega(0:nx-1, 0:ny-1)
    real(c_double), intent(in)  :: nu, dt
    real(c_double), intent(in)  :: kx(0:nx-1, 0:ny-1), ky(0:nx-1, 0:ny-1)
    logical,        intent(in)  :: mask(0:nx-1, 0:ny-1)
    real(c_double), intent(out) :: omega_out(0:nx-1, 0:ny-1)

    real(c_double) :: k1(0:nx-1,0:ny-1), k2m(0:nx-1,0:ny-1)
    real(c_double) :: k3(0:nx-1,0:ny-1), k4(0:nx-1,0:ny-1)
    real(c_double) :: tmp(0:nx-1,0:ny-1)
    complex(8)     :: wh(0:nx-1,0:ny-1), ph(0:nx-1,0:ny-1)
    real(c_double) :: u(0:nx-1,0:ny-1),  v(0:nx-1,0:ny-1)

    ! k1
    wh = cmplx(omega, 0.0_c_double, kind=8)
    call kh_fft_forward_2d(wh, nx, ny)
    call solve_poisson(wh, kx, ky, nx, ny, ph)
    call velocity_from_psi(ph, kx, ky, nx, ny, u, v)
    call vorticity_rhs(omega, u, v, nu, kx, ky, mask, nx, ny, k1)

    ! k2
    tmp = omega + 0.5_c_double * dt * k1
    wh = cmplx(tmp, 0.0_c_double, kind=8)
    call kh_fft_forward_2d(wh, nx, ny)
    call solve_poisson(wh, kx, ky, nx, ny, ph)
    call velocity_from_psi(ph, kx, ky, nx, ny, u, v)
    call vorticity_rhs(tmp, u, v, nu, kx, ky, mask, nx, ny, k2m)

    ! k3
    tmp = omega + 0.5_c_double * dt * k2m
    wh = cmplx(tmp, 0.0_c_double, kind=8)
    call kh_fft_forward_2d(wh, nx, ny)
    call solve_poisson(wh, kx, ky, nx, ny, ph)
    call velocity_from_psi(ph, kx, ky, nx, ny, u, v)
    call vorticity_rhs(tmp, u, v, nu, kx, ky, mask, nx, ny, k3)

    ! k4
    tmp = omega + dt * k3
    wh = cmplx(tmp, 0.0_c_double, kind=8)
    call kh_fft_forward_2d(wh, nx, ny)
    call solve_poisson(wh, kx, ky, nx, ny, ph)
    call velocity_from_psi(ph, kx, ky, nx, ny, u, v)
    call vorticity_rhs(tmp, u, v, nu, kx, ky, mask, nx, ny, k4)

    omega_out = omega + (dt / 6.0_c_double) * (k1 + 2.0_c_double*k2m + 2.0_c_double*k3 + k4)
  end subroutine rk4_step


  ! ── RK4 step with additive forcing and linear drag (2026-10-10) ───────────────
  ! As rk4_step, every stage right-hand side becomes (rhs - drag*w) + F, the order of
  ! the oracle's kh_physics.rk4_step: F at t (f1), t + dt/2 (f2, stages 2 and 3) and
  ! t + dt (f4). Used only when forcing or drag is requested, so an unforced run
  ! stays bit-identical to rk4_step (golden vectors).
  subroutine rk4_step_fd(omega, nu, kx, ky, mask, nx, ny, dt, drag, f1, f2, f4, omega_out)
    integer,        intent(in)  :: nx, ny
    real(c_double), intent(in)  :: omega(0:nx-1, 0:ny-1)
    real(c_double), intent(in)  :: nu, dt, drag
    real(c_double), intent(in)  :: kx(0:nx-1, 0:ny-1), ky(0:nx-1, 0:ny-1)
    logical,        intent(in)  :: mask(0:nx-1, 0:ny-1)
    real(c_double), intent(in)  :: f1(0:nx-1, 0:ny-1), f2(0:nx-1, 0:ny-1), f4(0:nx-1, 0:ny-1)
    real(c_double), intent(out) :: omega_out(0:nx-1, 0:ny-1)

    real(c_double) :: k1(0:nx-1,0:ny-1), k2m(0:nx-1,0:ny-1)
    real(c_double) :: k3(0:nx-1,0:ny-1), k4(0:nx-1,0:ny-1)
    real(c_double) :: tmp(0:nx-1,0:ny-1)

    call stage(omega, f1, k1)
    tmp = omega + 0.5_c_double * dt * k1
    call stage(tmp, f2, k2m)
    tmp = omega + 0.5_c_double * dt * k2m
    call stage(tmp, f2, k3)
    tmp = omega + dt * k3
    call stage(tmp, f4, k4)
    omega_out = omega + (dt / 6.0_c_double) * (k1 + 2.0_c_double*k2m + 2.0_c_double*k3 + k4)

  contains

    subroutine stage(w, f, r)
      real(c_double), intent(in)  :: w(0:nx-1, 0:ny-1), f(0:nx-1, 0:ny-1)
      real(c_double), intent(out) :: r(0:nx-1, 0:ny-1)
      complex(8)     :: wh(0:nx-1,0:ny-1), ph(0:nx-1,0:ny-1)
      real(c_double) :: u(0:nx-1,0:ny-1),  v(0:nx-1,0:ny-1)
      wh = cmplx(w, 0.0_c_double, kind=8)
      call kh_fft_forward_2d(wh, nx, ny)
      call solve_poisson(wh, kx, ky, nx, ny, ph)
      call velocity_from_psi(ph, kx, ky, nx, ny, u, v)
      call vorticity_rhs(w, u, v, nu, kx, ky, mask, nx, ny, r)
      r = (r - drag * w) + f
    end subroutine stage

  end subroutine rk4_step_fd


  ! ── Band forcing at time t from separable tables (2026-10-10) ─────────────────
  ! F(x,y,t) = (A/sqrt(M)) * sum_m cos(a_m x + b_m y + phi0_m + Om_m t),
  ! a_m = 2 pi ix_m / lx, b_m = 2 pi iy_m / ly (kh_forcing / oracle convention), split as
  ! cos(a x) cos(b y + th) - sin(a x) sin(b y + th): M*(nx + ny) cosines per evaluation
  ! instead of M*nx*ny (kh_forcing.eval_band_forcing evaluates every point).
  subroutine band_at(nx, ny, nmode, fscale, om, cx, sx, cy, sy, t, f)
    integer,        intent(in)  :: nx, ny, nmode
    real(c_double), intent(in)  :: fscale, om(nmode), t
    real(c_double), intent(in)  :: cx(0:nx-1, nmode), sx(0:nx-1, nmode)
    real(c_double), intent(in)  :: cy(0:ny-1, nmode), sy(0:ny-1, nmode)
    real(c_double), intent(out) :: f(0:nx-1, 0:ny-1)
    real(c_double) :: ct, st, cb, sb
    integer        :: m, i, j
    f = 0.0_c_double
    do m = 1, nmode
      ct = cos(om(m) * t)
      st = sin(om(m) * t)
      do j = 0, ny-1
        cb = cy(j,m) * ct - sy(j,m) * st
        sb = sy(j,m) * ct + cy(j,m) * st
        do i = 0, nx-1
          f(i,j) = f(i,j) + (cx(i,m) * cb - sx(i,m) * sb)
        end do
      end do
    end do
    f = fscale * f
  end subroutine band_at


  ! ── C-interop entry point ─────────────────────────────────────────────────────
  ! Called from kh_shim.cpp via iso_c_binding.
  ! All arrays passed as flat C double* (row-major from C side).

  subroutine kh_simulate_c(        &
      c_nx, c_ny,                  &
      c_lx, c_ly,                  &
      c_dt, c_steps,               &
      c_re, c_u0, c_amp, c_mode,   &
      c_init_omega,                &
      c_ftype, c_famp,             &
      c_fkx, c_fky,                &
      c_fkmin, c_fkmax,            &
      c_fomega0, c_drag, c_t_start, &
      c_u, c_v, c_omega, c_psi,    &
      c_ke, c_enstrophy,           &
      c_max_vort, c_div_rms,       &
      c_cfl_peak                   &
    ) bind(c, name="kh_simulate_c")

    integer(c_int),  value, intent(in)  :: c_nx, c_ny, c_steps, c_mode
    real(c_double),  value, intent(in)  :: c_lx, c_ly, c_dt, c_re, c_u0, c_amp
    type(c_ptr),     value, intent(in)  :: c_init_omega
    ! Forcing (2026-10-10, ruling 8): c_ftype 0 none, 1 mode, 2 band; the shim has
    ! validated the request (contract KH-VAL-2) and passes c_ftype 0 when the
    ! amplitude is 0. c_t_start: the time of the initial field.
    integer(c_int),  value, intent(in)  :: c_ftype, c_fkx, c_fky
    real(c_double),  value, intent(in)  :: c_famp, c_fkmin, c_fkmax, c_fomega0, c_drag, c_t_start
    real(c_double),         intent(out) :: c_u(0:c_nx*c_ny-1)
    real(c_double),         intent(out) :: c_v(0:c_nx*c_ny-1)
    real(c_double),         intent(out) :: c_omega(0:c_nx*c_ny-1)
    real(c_double),         intent(out) :: c_psi(0:c_nx*c_ny-1)
    real(c_double),         intent(out) :: c_ke, c_enstrophy, c_max_vort, c_div_rms
    real(c_double),         intent(out) :: c_cfl_peak

    integer        :: nx, ny, steps, mode, i, j, s
    real(c_double) :: lx, ly, dt, nu, delta, dx, dy
    real(c_double), allocatable :: omega(:,:), omega_next(:,:)
    real(c_double), allocatable :: u(:,:), v(:,:), psi_r(:,:)
    real(c_double), allocatable :: kx(:,:), ky(:,:), k2(:,:)
    complex(8),     allocatable :: oh(:,:), ph(:,:)
    logical,        allocatable :: mask(:,:)
    real(c_double), pointer     :: init_ptr(:) => null()
    ! CFL monitoring
    real(c_double) :: dt_cfl_adv, dt_cfl_diff, dt_cfl, u_max, v_max, cfl_step
    ! Diagnostics
    real(c_double) :: ke_sum, ens_sum, max_w, div_sum, n_inv
    complex(8),     allocatable :: uh(:,:), vh(:,:), div_s(:,:)
    real(c_double) :: div_r(0:c_nx-1, 0:c_ny-1)
    ! Forcing
    logical        :: forced, rotating
    real(c_double) :: t, two_pi, fscale, a, b
    real(c_double), allocatable :: f1(:,:), f2(:,:), f4(:,:)
    real(c_double), allocatable :: cx(:,:), sx(:,:), cy(:,:), sy(:,:), om(:)
    integer        :: nmode, m
    integer        :: ixa(MAX_BAND_MODES), iya(MAX_BAND_MODES)
    real(c_double) :: phi0a(MAX_BAND_MODES), oma(MAX_BAND_MODES)

    nx = c_nx; ny = c_ny; steps = c_steps; mode = c_mode
    lx = c_lx; ly = c_ly; dt = c_dt
    ! nu = |U| / Re: a reversed shear is the mirrored flow, not negative viscosity
    ! (modelling rules v0.1, K3; ruling 2026-10-07). Identical for U >= 0.
    nu    = abs(c_u0) / c_re
    delta = 0.05_c_double * ly
    dx    = lx / real(nx, c_double)
    dy    = ly / real(ny, c_double)

    allocate(omega(0:nx-1, 0:ny-1), omega_next(0:nx-1, 0:ny-1))
    allocate(kx(0:nx-1, 0:ny-1), ky(0:nx-1, 0:ny-1), k2(0:nx-1, 0:ny-1))
    allocate(mask(0:nx-1, 0:ny-1))
    allocate(u(0:nx-1, 0:ny-1), v(0:nx-1, 0:ny-1), psi_r(0:nx-1, 0:ny-1))
    allocate(oh(0:nx-1, 0:ny-1), ph(0:nx-1, 0:ny-1))
    allocate(uh(0:nx-1, 0:ny-1), vh(0:nx-1, 0:ny-1), div_s(0:nx-1, 0:ny-1))

    ! Build wavenumbers and 2/3-rule dealiasing mask
    call kh_grid_wavenums(nx, ny, dx, dy, kx, ky, k2)
    call kh_grid_dealias_mask(nx, ny, kx, ky, dx, dy, mask)

    ! Initial vorticity: from C pointer (row-major i*ny+j) or analytic IC.
    ! c_f_pointer gives init_ptr the default lower bound 1, so C element i*ny+j is
    ! init_ptr(i*ny + j + 1). Until 2026-10-03 the +1 was missing: every
    ! initial_omega was read shifted by one cell and the last read fell outside
    ! the array, which -fcheck=bounds turned into a process abort (finding v0.5, F-2).
    if (c_associated(c_init_omega)) then
      call c_f_pointer(c_init_omega, init_ptr, [nx*ny])
      do j = 0, ny-1
        do i = 0, nx-1
          omega(i,j) = init_ptr(i*ny + j + 1)
        end do
      end do
    else
      call kh_ic_tanh_shear(nx, ny, lx, ly, c_u0, delta, c_amp, mode, omega)
    end if

    ! CFL guard (KH-NM-5): warn once at start if dt exceeds stability bound
    oh = cmplx(omega, 0.0_c_double, kind=8)
    call kh_fft_forward_2d(oh, nx, ny)
    call solve_poisson(oh, kx, ky, nx, ny, ph)
    call velocity_from_psi(ph, kx, ky, nx, ny, u, v)
    u_max = maxval(abs(u))
    v_max = maxval(abs(v))
    dt_cfl_adv  = min(dx / (u_max + 1.0e-30_c_double), dy / (v_max + 1.0e-30_c_double))
    ! nu = 0 (inviscid) is a valid limit: no diffusive bound (modelling rules v0.1, D).
    if (nu == 0.0_c_double) then
      dt_cfl_diff = huge(1.0_c_double)
    else
      dt_cfl_diff = min(dx**2 / (2.0_c_double * nu), dy**2 / (2.0_c_double * nu))
    end if
    dt_cfl      = min(dt_cfl_adv, dt_cfl_diff)
    if (dt > dt_cfl) then
      write(*,'(a,e12.4,a,e12.4)') &
        'KH-NM-5: dt=', dt, ' > CFL bound ', dt_cfl
      write(*,*) '  simulation may be unstable'
    end if

    ! Forcing tables (none unless requested): static field for "mode" and for a band
    ! without rotation, separable cos/sin tables for a rotating band.
    forced   = (c_ftype /= FORCING_NONE) .or. (c_drag /= 0.0_c_double)
    rotating = .false.
    nmode    = 0
    if (forced) then
      allocate(f1(0:nx-1, 0:ny-1), f2(0:nx-1, 0:ny-1), f4(0:nx-1, 0:ny-1))
      f1 = 0.0_c_double
      if (c_ftype == FORCING_MODE) then
        call build_mode_forcing(nx, ny, lx, ly, c_famp, c_fkx, c_fky, f1)
      else if (c_ftype == FORCING_BAND) then
        call build_band_phi0(c_fkmin, c_fkmax, nmode, ixa, iya, phi0a)
        call build_band_omegas(nmode, ixa, iya, c_fomega0, oma)
        if (nmode > 0) then
          two_pi = 8.0_c_double * atan(1.0_c_double)
          fscale  = c_famp / sqrt(real(nmode, c_double))
          allocate(cx(0:nx-1, nmode), sx(0:nx-1, nmode), cy(0:ny-1, nmode), sy(0:ny-1, nmode))
          allocate(om(nmode))
          do m = 1, nmode
            do i = 0, nx-1
              a = two_pi * ixa(m) * (i * (lx / nx)) / lx
              cx(i,m) = cos(a);  sx(i,m) = sin(a)
            end do
            do j = 0, ny-1
              b = two_pi * iya(m) * (j * (ly / ny)) / ly + phi0a(m)
              cy(j,m) = cos(b);  sy(j,m) = sin(b)
            end do
            om(m) = oma(m)
          end do
          rotating = (c_fomega0 /= 0.0_c_double)
          if (.not. rotating) call band_at(nx, ny, nmode, fscale, om, cx, sx, cy, sy, 0.0_c_double, f1)
        end if
      end if
      f2 = f1
      f4 = f1
    end if

    ! Time integration — pre-step CFL tracking (matches Python oracle timing).
    ! t starts at t_start: a continued run keeps the forcing phase.
    t = c_t_start
    c_cfl_peak = 0.0_c_double
    do s = 1, steps
      oh = cmplx(omega, 0.0_c_double, kind=8)
      call kh_fft_forward_2d(oh, nx, ny)
      call solve_poisson(oh, kx, ky, nx, ny, ph)
      call velocity_from_psi(ph, kx, ky, nx, ny, u, v)
      u_max    = maxval(abs(u))
      v_max    = maxval(abs(v))
      cfl_step = dt * (u_max / dx + v_max / dy)
      if (cfl_step > c_cfl_peak) c_cfl_peak = cfl_step
      if (forced) then
        if (rotating) then
          call band_at(nx, ny, nmode, fscale, om, cx, sx, cy, sy, t, f1)
          call band_at(nx, ny, nmode, fscale, om, cx, sx, cy, sy, t + 0.5_c_double * dt, f2)
          call band_at(nx, ny, nmode, fscale, om, cx, sx, cy, sy, t + dt, f4)
        end if
        call rk4_step_fd(omega, nu, kx, ky, mask, nx, ny, dt, c_drag, f1, f2, f4, omega_next)
      else
        call rk4_step(omega, nu, kx, ky, mask, nx, ny, dt, omega_next)
      end if
      omega = omega_next
      t = t + dt
    end do

    ! Final field recovery
    oh = cmplx(omega, 0.0_c_double, kind=8)
    call kh_fft_forward_2d(oh, nx, ny)
    call solve_poisson(oh, kx, ky, nx, ny, ph)
    call kh_fft_inverse_2d(ph, nx, ny)
    psi_r = real(ph, kind=c_double)
    ! recompute ph for velocity (it was overwritten by ifft)
    oh = cmplx(omega, 0.0_c_double, kind=8)
    call kh_fft_forward_2d(oh, nx, ny)
    call solve_poisson(oh, kx, ky, nx, ny, ph)
    call velocity_from_psi(ph, kx, ky, nx, ny, u, v)

    ! Diagnostics
    n_inv = 1.0_c_double / real(nx*ny, c_double)
    ke_sum  = 0.0_c_double; ens_sum = 0.0_c_double; max_w = 0.0_c_double
    do j = 0, ny-1
      do i = 0, nx-1
        ke_sum  = ke_sum  + 0.5_c_double*(u(i,j)**2 + v(i,j)**2)
        ens_sum = ens_sum + 0.5_c_double*omega(i,j)**2
        if (abs(omega(i,j)) > max_w) max_w = abs(omega(i,j))
      end do
    end do
    c_ke        = ke_sum  * n_inv
    c_enstrophy = ens_sum * n_inv
    c_max_vort  = max_w

    ! Divergence rms
    uh = cmplx(u, 0.0_c_double, kind=8);  call kh_fft_forward_2d(uh, nx, ny)
    vh = cmplx(v, 0.0_c_double, kind=8);  call kh_fft_forward_2d(vh, nx, ny)
    do j = 0, ny-1
      do i = 0, nx-1
        div_s(i,j) = I_UNIT*kx(i,j)*uh(i,j) + I_UNIT*ky(i,j)*vh(i,j)
      end do
    end do
    call kh_fft_inverse_2d(div_s, nx, ny)
    div_r = real(div_s, kind=c_double)
    div_sum = 0.0_c_double
    do j = 0, ny-1
      do i = 0, nx-1
        div_sum = div_sum + div_r(i,j)**2
      end do
    end do
    c_div_rms = sqrt(div_sum * n_inv)

    ! Copy fields to flat C arrays (row-major: index = i*ny + j)
    do j = 0, ny-1
      do i = 0, nx-1
        c_u(i*ny+j)     = u(i,j)
        c_v(i*ny+j)     = v(i,j)
        c_omega(i*ny+j) = omega(i,j)
        c_psi(i*ny+j)   = psi_r(i,j)
      end do
    end do

    deallocate(omega, omega_next, kx, ky, k2, mask, u, v, psi_r, oh, ph, uh, vh, div_s)
    if (allocated(f1)) deallocate(f1, f2, f4)
    if (allocated(cx)) deallocate(cx, sx, cy, sy, om)
  end subroutine kh_simulate_c


  ! ── Fortran-friendly RK4 driver ───────────────────────────────────────────────
  ! Pure-Fortran parameters — no iso_c_binding required from callers.
  ! Mirrors kh_solver_run's diagnostic outputs so tests 006/008 can swap in.
  ! The C HTTP-shim still calls kh_simulate_c directly.

  subroutine kh_simulate_rk4(nx, ny, lx, ly, nu, dt, nsteps, amp, mode, &
                               ke_out, enstrophy_out, max_vort_out, div_rms_out, cfl_peak)
    integer,        intent(in)  :: nx, ny, nsteps, mode
    real(c_double), intent(in)  :: lx, ly, nu, dt, amp
    real(c_double), intent(out) :: ke_out, enstrophy_out, max_vort_out, div_rms_out
    real(c_double), intent(out) :: cfl_peak

    real(c_double), allocatable :: x(:), y(:)
    real(c_double), allocatable :: kx(:,:), ky(:,:), k2(:,:)
    logical,        allocatable :: mask(:,:)
    real(c_double), allocatable :: omega(:,:), omega_next(:,:)
    real(c_double), allocatable :: u(:,:), v(:,:)
    complex(8),     allocatable :: oh(:,:), ph(:,:)

    real(c_double) :: dx, dy, delta
    real(c_double) :: dt_cfl_adv, dt_cfl_diff, dt_cfl, u_max, v_max, cfl
    integer        :: s

    allocate(x(0:nx-1), y(0:ny-1))
    allocate(kx(0:nx-1, 0:ny-1), ky(0:nx-1, 0:ny-1), k2(0:nx-1, 0:ny-1))
    allocate(mask(0:nx-1, 0:ny-1))
    allocate(omega(0:nx-1, 0:ny-1), omega_next(0:nx-1, 0:ny-1))
    allocate(u(0:nx-1, 0:ny-1), v(0:nx-1, 0:ny-1))
    allocate(oh(0:nx-1, 0:ny-1), ph(0:nx-1, 0:ny-1))

    call kh_grid_make(nx, ny, lx, ly, x, y, dx, dy)
    call kh_grid_wavenums(nx, ny, dx, dy, kx, ky, k2)
    call kh_grid_dealias_mask(nx, ny, kx, ky, dx, dy, mask)

    delta = KH_DELTA_FACTOR * ly
    call kh_ic_tanh_shear(nx, ny, lx, ly, KH_U0_DEFAULT, delta, amp, mode, omega)

    ! CFL guard (KH-NM-5): warn once at start if dt exceeds stability bound
    oh = cmplx(omega, 0.0_c_double, kind=8)
    call kh_fft_forward_2d(oh, nx, ny)
    call solve_poisson(oh, kx, ky, nx, ny, ph)
    call velocity_from_psi(ph, kx, ky, nx, ny, u, v)
    u_max = maxval(abs(u))
    v_max = maxval(abs(v))
    dt_cfl_adv  = min(dx / (u_max + 1.0e-30_c_double), dy / (v_max + 1.0e-30_c_double))
    ! nu = 0 (inviscid) is a valid limit: no diffusive bound (modelling rules v0.1, D).
    if (nu == 0.0_c_double) then
      dt_cfl_diff = huge(1.0_c_double)
    else
      dt_cfl_diff = min(dx**2 / (2.0_c_double * nu), dy**2 / (2.0_c_double * nu))
    end if
    dt_cfl = min(dt_cfl_adv, dt_cfl_diff)
    if (dt > dt_cfl) then
      write(*,'(a,e12.4,a,e12.4)') &
        'KH-NM-5: dt=', dt, ' > CFL bound ', dt_cfl
      write(*,*) '  simulation may be unstable'
    end if

    ! Time integration — RK4 in physical space; track peak CFL per step
    cfl_peak = 0.0_c_double
    do s = 1, nsteps
      call rk4_step(omega, nu, kx, ky, mask, nx, ny, dt, omega_next)
      omega = omega_next
      oh = cmplx(omega, 0.0_c_double, kind=8)
      call kh_fft_forward_2d(oh, nx, ny)
      call solve_poisson(oh, kx, ky, nx, ny, ph)
      call velocity_from_psi(ph, kx, ky, nx, ny, u, v)
      cfl = dt * (maxval(abs(u)) / dx + maxval(abs(v)) / dy)
      if (cfl > cfl_peak) cfl_peak = cfl
    end do

    ! Final diagnostics — u, v already current from last CFL step
    call kh_diagnostics_compute(omega, u, v, kx, ky, nx, ny, &
                                 ke_out, enstrophy_out, max_vort_out, div_rms_out)

    deallocate(x, y, kx, ky, k2, mask, omega, omega_next, u, v, oh, ph)
  end subroutine kh_simulate_rk4

end module kh_physics
