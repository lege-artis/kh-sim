! kh_io.f90 — KH instability: namelist parameter reader + JSON snapshot writer
!
! REFERENCE WATERMARK — kh-sim Fortran reference, v0.1;
! canonical for: TC-NUM-KH-001..008;
! do not modify without bumping reference version + regenerating reference outputs.
!
! Provides two facilities:
!
!  (A) Parameter I/O — reads simulation configuration from a Fortran namelist
!      file (default: kh_params.nml). Falls back to kh_constants.f90 defaults
!      if file is absent or malformed — never aborts.
!
!      Namelist format (file: kh_params.nml):
!        &KH_PARAMS
!          nx           = 128,     ! grid points x  (power of 2)
!          ny           = 64,      ! grid points y  (power of 2)
!          lx           = 1.0,     ! domain length x [non-dim]
!          ly           = 0.5,     ! domain length y [non-dim]
!          re           = 1000.0,  ! Reynolds number (0 = inviscid)
!          dt           = 0.001,   ! time step
!          nsteps       = 1000,    ! total time steps
!          amp          = 0.01,    ! perturbation amplitude
!          mode         = 2,       ! perturbation wavenumber
!          out_interval = 10,      ! steps between snapshot writes
!        /
!
!  (B) JSON snapshot writer — writes one JSON object per line to an output
!      file.  Each line contains step, time, and the four scalar diagnostics
!      computed by kh_diagnostics_compute.
!
!      Output line format:
!        {"step":N,"time":T,"ke":K,"enstrophy":E,"max_vort":M,"div_rms":D}
!      where numeric values are in ES20.12 format.
!
! Public interface:
!   kh_io_read_params(filename, nx, ny, lx, ly, re, dt, nsteps, amp, mode, out_interval)
!   kh_io_open_output(filename, unit)
!   kh_io_write_snapshot(unit, step, time, ke, enstrophy, max_vort, div_rms)
!   kh_io_close_output(unit)
!
! References:
!   PHYSICS-NUMERICAL-METHODS-v0.1.md §2.5 (I/O spec)
!   kh-sim/shared/physics/KH-PHYSICS.md §3 (default parameters)

module kh_io
  use kh_constants, only: c_double, &
      KH_NX_DEFAULT, KH_NY_DEFAULT, &
      KH_LX_DEFAULT, KH_LY_DEFAULT, &
      KH_RE_DEFAULT, KH_DT_DEFAULT, KH_STEPS_DEFAULT, &
      KH_AMP_DEFAULT, KH_MODE_DEFAULT
  implicit none
  private

  public :: kh_io_read_params
  public :: kh_io_open_output
  public :: kh_io_write_snapshot
  public :: kh_io_close_output

contains

  ! ── Read simulation parameters from namelist file ────────────────────────────
  !
  ! If filename cannot be opened, or the namelist block is malformed, a warning
  ! is printed to stdout and all parameters take their kh_constants defaults.
  ! Re=0 is permitted (inviscid run; L=0 throughout, ETDRK4 → standard RK4).
  !
  ! Arguments:
  !   filename       [in]  path to namelist file (e.g. "kh_params.nml")
  !   nx, ny         [out] grid dimensions
  !   lx, ly         [out] domain extents
  !   re             [out] Reynolds number (0 = inviscid)
  !   dt             [out] time step
  !   nsteps         [out] total steps
  !   amp            [out] perturbation amplitude
  !   mode           [out] perturbation wavenumber mode
  !   out_interval   [out] steps between snapshot writes

  subroutine kh_io_read_params(filename, nx, ny, lx, ly, re, dt, &
                                nsteps, amp, mode, out_interval)
    character(len=*), intent(in)  :: filename
    integer,          intent(out) :: nx, ny, nsteps, mode, out_interval
    real(c_double),   intent(out) :: lx, ly, re, dt, amp

    ! Namelist variables (need separate names to avoid aliasing)
    integer        :: nx_nml, ny_nml, nsteps_nml, mode_nml, out_interval_nml
    real(c_double) :: lx_nml, ly_nml, re_nml, dt_nml, amp_nml
    integer        :: nml_unit, ios

    namelist /KH_PARAMS/ nx_nml, ny_nml, lx_nml, ly_nml, re_nml, &
                         dt_nml, nsteps_nml, amp_nml, mode_nml, out_interval_nml

    ! Initialise to defaults (applied if file absent or malformed)
    nx_nml           = KH_NX_DEFAULT
    ny_nml           = KH_NY_DEFAULT
    lx_nml           = KH_LX_DEFAULT
    ly_nml           = KH_LY_DEFAULT
    re_nml           = KH_RE_DEFAULT
    dt_nml           = KH_DT_DEFAULT
    nsteps_nml       = KH_STEPS_DEFAULT
    amp_nml          = KH_AMP_DEFAULT
    mode_nml         = KH_MODE_DEFAULT
    out_interval_nml = 10

    open(newunit=nml_unit, file=trim(filename), status='old', &
         action='read', iostat=ios)
    if (ios /= 0) then
      write(*, '(A,A,A)') &
          "kh_io: cannot open '", trim(filename), "' — using defaults"
    else
      read(nml_unit, nml=KH_PARAMS, iostat=ios)
      if (ios /= 0) then
        write(*, '(A,A,A)') &
            "kh_io: namelist read error in '", trim(filename), "' — using defaults"
      end if
      close(nml_unit)
    end if

    nx           = nx_nml
    ny           = ny_nml
    lx           = lx_nml
    ly           = ly_nml
    re           = re_nml
    dt           = dt_nml
    nsteps       = nsteps_nml
    amp          = amp_nml
    mode         = mode_nml
    out_interval = out_interval_nml

  end subroutine kh_io_read_params


  ! ── Open JSON output file ────────────────────────────────────────────────────
  !
  ! Creates or truncates the file.  Returns unit > 0 on success, -1 on failure.
  ! Caller passes the returned unit to kh_io_write_snapshot and kh_io_close_output.

  subroutine kh_io_open_output(filename, unit)
    character(len=*), intent(in)  :: filename
    integer,          intent(out) :: unit
    integer :: ios

    open(newunit=unit, file=trim(filename), status='replace', &
         action='write', iostat=ios)
    if (ios /= 0) then
      write(*, '(A,A,A)') &
          "kh_io: cannot open output '", trim(filename), "'"
      unit = -1
    end if
  end subroutine kh_io_open_output


  ! ── Write one JSON snapshot line ─────────────────────────────────────────────
  !
  ! Appends a single JSON object (one line) to the open file unit.
  ! No-ops if unit <= 0 (allows callers to skip opening gracefully).
  !
  ! Output format (ES20.12 for all reals):
  !   {"step":N,"time":T,"ke":K,"enstrophy":E,"max_vort":M,"div_rms":D}

  subroutine kh_io_write_snapshot(unit, step, time, ke, enstrophy, max_vort, div_rms)
    integer,        intent(in) :: unit, step
    real(c_double), intent(in) :: time, ke, enstrophy, max_vort, div_rms

    if (unit <= 0) return

    write(unit, '(A,I0,A,ES20.12,A,ES20.12,A,ES20.12,A,ES20.12,A,ES20.12,A)') &
        '{"step":', step,       &
        ',"time":', time,       &
        ',"ke":', ke,           &
        ',"enstrophy":', enstrophy, &
        ',"max_vort":', max_vort,   &
        ',"div_rms":', div_rms,     &
        '}'

  end subroutine kh_io_write_snapshot


  ! ── Close JSON output file ───────────────────────────────────────────────────

  subroutine kh_io_close_output(unit)
    integer, intent(in) :: unit
    if (unit > 0) close(unit)
  end subroutine kh_io_close_output

end module kh_io
