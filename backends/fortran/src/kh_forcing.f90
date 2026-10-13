! kh_forcing.f90 -- phase-7f forcing field construction (Fortran 2008)
! Ref: PHASE-7F-FORCING-API-DESIGN-v0.1.md §3+§7
! Lesieur(2008) Turbulence in Fluids 4th ed ch.9 §9.3  -- CP-7F-B anchor
! Brachet et al.(1983) J. Fluid Mech. 130:411-452       -- CP-7F-D anchor
! 2026-10-10 (ruling 8, CP-7F-B/D harness): kh_physics.f90 uses the builders below,
! so the HTTP server computes forcing too (it refused forcing until then).
module kh_forcing
  use iso_c_binding, only: c_double
  implicit none
  private
  public :: build_mode_forcing, build_band_phi0
  public :: build_band_omegas, eval_band_forcing, rk4_step_forced
  public :: MAX_BAND_MODES
  integer, parameter :: MAX_BAND_MODES = 4096
contains

  ! build_mode_forcing -- CP-7F-D steady Taylor-Green-mode forcing
  ! F(x,y) = A*cos(2*pi*kx_mode*x/lx)*cos(2*pi*ky_mode*y/ly)
  ! Ref: PHASE-7F-FORCING-API-DESIGN-v0.1.md §3 "mode"
  ! Brachet et al.(1983) J. Fluid Mech. 130:411 -- forced Taylor-Green vortex
  subroutine build_mode_forcing(nlenx, nleny, lx, ly, amp, kx_mode, ky_mode, forcing)
    integer,        intent(in)  :: nlenx, nleny, kx_mode, ky_mode
    real(c_double), intent(in)  :: lx, ly, amp
    real(c_double), intent(out) :: forcing(0:nlenx-1, 0:nleny-1)
    integer        :: ii, jj
    real(c_double) :: two_pi, xi, yj
    two_pi = 8.0_c_double * atan(1.0_c_double)
    do jj = 0, nleny-1
      yj = jj * ly / nleny
      do ii = 0, nlenx-1
        xi = ii * lx / nlenx
        forcing(ii,jj) = amp * cos(two_pi * kx_mode * xi / lx) &
                             * cos(two_pi * ky_mode * yj / ly)
      end do
    end do
  end subroutine build_mode_forcing

  ! build_band_phi0 -- deterministic phase seeds for band forcing
  ! phi0(ix,iy) = 2*pi*frac(ix*sqrt(2) + iy*sqrt(3))
  ! Ref: PHASE-7F-FORCING-API-DESIGN-v0.1.md §3 "band"
  ! Half-plane: ix>=0; skip (ix==0 and iy<=0).
  subroutine build_band_phi0(kmin, kmax, nmode_out, ix_arr, iy_arr, phi0_arr)
    real(c_double), intent(in)  :: kmin, kmax
    integer,        intent(out) :: nmode_out
    integer,        intent(out) :: ix_arr(MAX_BAND_MODES)
    integer,        intent(out) :: iy_arr(MAX_BAND_MODES)
    real(c_double), intent(out) :: phi0_arr(MAX_BAND_MODES)
    real(c_double) :: two_pi, sqrt2, sqrt3, kr, fv
    integer        :: ix_max, ix, iy, mm
    two_pi = 8.0_c_double * atan(1.0_c_double)
    sqrt2  = sqrt(2.0_c_double)
    sqrt3  = sqrt(3.0_c_double)
    ix_max = int(kmax) + 1
    mm = 0
    do ix = 0, ix_max
      do iy = -(ix_max), ix_max
        if (ix == 0 .and. iy <= 0) cycle
        kr = sqrt(real(ix,c_double)**2 + real(iy,c_double)**2)
        if (kr < kmin .or. kr > kmax) cycle
        mm = mm + 1
        if (mm > MAX_BAND_MODES) then
          write(*,'(A)') 'ERROR: build_band_phi0 MAX_BAND_MODES exceeded'
          stop 1
        end if
        ix_arr(mm)   = ix
        iy_arr(mm)   = iy
        fv           = ix * sqrt2 + iy * sqrt3
        fv           = fv - floor(fv)
        phi0_arr(mm) = two_pi * fv
      end do
    end do
    nmode_out = mm
  end subroutine build_band_phi0

  ! build_band_omegas -- per-mode rotation rates for v0.2 rotating-phase forcing
  ! Omega(ix,iy) = 2*pi*frac(ix*sqrt(11) + iy*sqrt(13))*forcing_omega0
  ! Ref: PHASE-7F-FORCING-API-DESIGN-v0.1.md §7
  ! Set forcing_omega0=0 for v0.1 time-independent behaviour.
  subroutine build_band_omegas(nmode, ix_arr, iy_arr, forcing_omega0, omega_arr)
    integer,        intent(in)  :: nmode
    integer,        intent(in)  :: ix_arr(MAX_BAND_MODES), iy_arr(MAX_BAND_MODES)
    real(c_double), intent(in)  :: forcing_omega0
    real(c_double), intent(out) :: omega_arr(MAX_BAND_MODES)
    real(c_double) :: two_pi, s11, s13, fv
    integer        :: mm
    two_pi = 8.0_c_double * atan(1.0_c_double)
    s11    = sqrt(11.0_c_double)
    s13    = sqrt(13.0_c_double)
    do mm = 1, nmode
      fv = ix_arr(mm) * s11 + iy_arr(mm) * s13
      fv = fv - floor(fv)
      omega_arr(mm) = two_pi * fv * forcing_omega0
    end do
  end subroutine build_band_omegas

  ! eval_band_forcing -- evaluate F(x,y,t) for one time instant
  ! F = (A/sqrt(M))*sum_m cos(2*pi*(ix*x/lx + iy*y/ly) + phi0_m + Omega_m*t)
  ! Ref: PHASE-7F-FORCING-API-DESIGN-v0.1.md §3+§7
  ! Set omg_arr(:)=0 for v0.1 time-independent behaviour.
  subroutine eval_band_forcing(nlenx, nleny, lx, ly, amp, nmode, &
                               ix_arr, iy_arr, phi0_arr, omg_arr, t_now, forcing)
    integer,        intent(in)  :: nlenx, nleny, nmode
    real(c_double), intent(in)  :: lx, ly, amp, t_now
    integer,        intent(in)  :: ix_arr(MAX_BAND_MODES), iy_arr(MAX_BAND_MODES)
    real(c_double), intent(in)  :: phi0_arr(MAX_BAND_MODES), omg_arr(MAX_BAND_MODES)
    real(c_double), intent(out) :: forcing(0:nlenx-1, 0:nleny-1)
    real(c_double) :: two_pi, scale, xi, yj, phase
    integer        :: ii, jj, mm
    two_pi = 8.0_c_double * atan(1.0_c_double)
    scale  = amp / sqrt(real(nmode, c_double))
    forcing = 0.0_c_double
    do mm = 1, nmode
      do jj = 0, nleny-1
        yj = jj * ly / nleny
        do ii = 0, nlenx-1
          xi = ii * lx / nlenx
          phase = two_pi * (ix_arr(mm)*xi/lx + iy_arr(mm)*yj/ly) &
                 + phi0_arr(mm) + omg_arr(mm)*t_now
          forcing(ii,jj) = forcing(ii,jj) + scale * cos(phase)
        end do
      end do
    end do
  end subroutine eval_band_forcing

  ! rk4_step_forced -- RK4 step with additive vorticity forcing + linear drag
  ! Mirrors rk4_step from kh_physics.f90 but adds:
  !   - band forcing evaluated at t, t+dt/2, t+dt
  !   - linear drag -alpha*omega at each stage
  ! use_static=.true.  => add forcing_static each stage (mode / CP-7F-D)
  ! use_static=.false. => evaluate band forcing dynamically at each stage time
  ! Ref: Press et al.(2007) Numerical Recipes 3rd ed §17.1 -- RK4 formula
  !      PHASE-7F-FORCING-API-DESIGN-v0.1.md §2, §7        -- forcing convention
  subroutine rk4_step_forced(omega, nu, kx, ky, k2arr, mask, nlenx, nleny, dt, &
                             alpha, nmode, ix_arr, iy_arr, phi0_arr, omg_arr, &
                             forcing_static, use_static, amp, lx, ly, t_now, &
                             omega_out)
    use kh_fft,      only: kh_fft_forward_2d, kh_fft_inverse_2d
    use kh_poisson,  only: kh_poisson_solve
    use kh_velocity, only: kh_velocity_from_psi
    use kh_nonlinear, only: kh_dealias
    integer,        intent(in)  :: nlenx, nleny, nmode
    real(c_double), intent(in)  :: omega(0:nlenx-1, 0:nleny-1)
    real(c_double), intent(in)  :: nu, dt, alpha, amp, lx, ly, t_now
    real(c_double), intent(in)  :: kx(0:nlenx-1, 0:nleny-1)
    real(c_double), intent(in)  :: ky(0:nlenx-1, 0:nleny-1)
    real(c_double), intent(in)  :: k2arr(0:nlenx-1, 0:nleny-1)
    logical,        intent(in)  :: mask(0:nlenx-1, 0:nleny-1)
    integer,        intent(in)  :: ix_arr(MAX_BAND_MODES), iy_arr(MAX_BAND_MODES)
    real(c_double), intent(in)  :: phi0_arr(MAX_BAND_MODES), omg_arr(MAX_BAND_MODES)
    real(c_double), intent(in)  :: forcing_static(0:nlenx-1, 0:nleny-1)
    logical,        intent(in)  :: use_static
    real(c_double), intent(out) :: omega_out(0:nlenx-1, 0:nleny-1)
    real(c_double) :: s1(0:nlenx-1,0:nleny-1), s2(0:nlenx-1,0:nleny-1)
    real(c_double) :: s3(0:nlenx-1,0:nleny-1), s4(0:nlenx-1,0:nleny-1)
    real(c_double) :: tmp(0:nlenx-1,0:nleny-1), fb(0:nlenx-1,0:nleny-1)

    call get_fb(t_now,                   fb)
    call rhs_loc(omega,                  fb, s1)
    tmp = omega + 0.5_c_double*dt*s1
    call get_fb(t_now + 0.5_c_double*dt, fb)
    call rhs_loc(tmp,                    fb, s2)
    tmp = omega + 0.5_c_double*dt*s2
    call get_fb(t_now + 0.5_c_double*dt, fb)
    call rhs_loc(tmp,                    fb, s3)
    tmp = omega + dt*s3
    call get_fb(t_now + dt,              fb)
    call rhs_loc(tmp,                    fb, s4)
    omega_out = omega + (dt/6.0_c_double) * &
                (s1 + 2.0_c_double*s2 + 2.0_c_double*s3 + s4)

  contains

    subroutine get_fb(ts, fo)
      real(c_double), intent(in)  :: ts
      real(c_double), intent(out) :: fo(0:nlenx-1, 0:nleny-1)
      if (use_static) then
        fo = forcing_static
      else if (nmode > 0) then
        call eval_band_forcing(nlenx, nleny, lx, ly, amp, nmode, &
                               ix_arr, iy_arr, phi0_arr, omg_arr, ts, fo)
      else
        fo = 0.0_c_double
      end if
    end subroutine get_fb

    subroutine rhs_loc(om, frc, rout)
      real(c_double), intent(in)  :: om(0:nlenx-1, 0:nleny-1)
      real(c_double), intent(in)  :: frc(0:nlenx-1, 0:nleny-1)
      real(c_double), intent(out) :: rout(0:nlenx-1, 0:nleny-1)
      complex(8) :: oh(0:nlenx-1,0:nleny-1), ph(0:nlenx-1,0:nleny-1)
      complex(8) :: dxh(0:nlenx-1,0:nleny-1), dyh(0:nlenx-1,0:nleny-1)
      complex(8) :: laph(0:nlenx-1,0:nleny-1), nlh(0:nlenx-1,0:nleny-1)
      real(c_double) :: ul(0:nlenx-1,0:nleny-1), vl(0:nlenx-1,0:nleny-1)
      real(c_double) :: dox(0:nlenx-1,0:nleny-1), doy(0:nlenx-1,0:nleny-1)
      real(c_double) :: lapr(0:nlenx-1,0:nleny-1), nlp(0:nlenx-1,0:nleny-1)
      complex(8), parameter :: IU = cmplx(0.0_c_double, 1.0_c_double, kind=8)
      integer        :: ii2, jj2
      oh = cmplx(om, 0.0_c_double, kind=8)
      call kh_fft_forward_2d(oh, nlenx, nleny)
      do jj2 = 0, nleny-1
        do ii2 = 0, nlenx-1
          dxh(ii2,jj2)  = IU * kx(ii2,jj2) * oh(ii2,jj2)
          dyh(ii2,jj2)  = IU * ky(ii2,jj2) * oh(ii2,jj2)
          laph(ii2,jj2) = -k2arr(ii2,jj2) * oh(ii2,jj2)
        end do
      end do
      call kh_fft_inverse_2d(dxh,  nlenx, nleny); dox  = real(dxh,  kind=c_double)
      call kh_fft_inverse_2d(dyh,  nlenx, nleny); doy  = real(dyh,  kind=c_double)
      call kh_fft_inverse_2d(laph, nlenx, nleny); lapr = real(laph, kind=c_double)
      call kh_poisson_solve(oh, k2arr, nlenx, nleny, ph)
      call kh_velocity_from_psi(ph, kx, ky, nlenx, nleny, ul, vl)
      nlh = cmplx(-(ul*dox + vl*doy), 0.0_c_double, kind=8)
      call kh_fft_forward_2d(nlh, nlenx, nleny)
      call kh_dealias(nlh, mask, nlenx, nleny)
      call kh_fft_inverse_2d(nlh, nlenx, nleny)
      nlp  = real(nlh, kind=c_double)
      rout = nlp + nu*lapr - alpha*om + frc
    end subroutine rhs_loc

  end subroutine rk4_step_forced

end module kh_forcing
