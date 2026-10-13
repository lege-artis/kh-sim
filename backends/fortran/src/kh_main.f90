! kh_main.f90 — KH instability: command-line entry point
!
! REFERENCE WATERMARK — kh-sim Fortran reference, v0.1;
! canonical for: TC-NUM-KH-001..008;
! do not modify without bumping reference version + regenerating reference outputs.
!
! Reads simulation parameters from kh_params.nml (default) or from the first
! command-line argument if provided.  Runs the solver and writes a JSON
! snapshot file (kh_out.json by default, or second CLI arg).
!
! Usage:
!   ./kh_main                                # reads kh_params.nml, writes kh_out.json
!   ./kh_main my_run.nml output.json         # explicit paths
!
! Namelist file (kh_params.nml):
!   &KH_PARAMS
!     nx=128, ny=64, lx=1.0, ly=0.5, re=1000.0, dt=0.001,
!     nsteps=1000, amp=0.01, mode=2, out_interval=10
!   /
!
! Output (kh_out.json): one JSON object per line for each snapshot step:
!   {"step":10,"time":0.010000...,"ke":...,"enstrophy":...,"max_vort":...,"div_rms":...}
!
! Exit code: 0 = normal completion
!
! Compilation:
!   gfortran -O2 -o kh_main \
!       src/kh_constants.f90 src/kh_grid.f90 src/kh_fft.f90 \
!       src/kh_poisson.f90 src/kh_velocity.f90 src/kh_nonlinear.f90 \
!       src/kh_etdrk4.f90 src/kh_diagnostics.f90 src/kh_io.f90 \
!       src/kh_solver.f90 src/kh_main.f90

program kh_main
  use kh_constants
  use kh_io,     only: kh_io_read_params, kh_io_open_output, kh_io_close_output
  use kh_solver, only: kh_solver_run
  implicit none

  ! ── Simulation parameters (filled by kh_io_read_params) ──────────────────────

  integer        :: nx, ny, nsteps, mode, out_interval
  real(c_double) :: lx, ly, re, nu, dt, amp

  ! ── I/O ──────────────────────────────────────────────────────────────────────

  integer            :: out_unit
  character(len=256) :: param_file, out_file

  ! ── Final diagnostics ─────────────────────────────────────────────────────────

  real(c_double) :: ke_final, enstrophy_final, max_vort_final, div_rms_final, cfl_peak

  ! ── Command-line argument handling ───────────────────────────────────────────

  integer :: nargs
  nargs = command_argument_count()

  if (nargs >= 1) then
    call get_command_argument(1, param_file)
  else
    param_file = 'kh_params.nml'
  end if

  if (nargs >= 2) then
    call get_command_argument(2, out_file)
  else
    out_file = 'kh_out.json'
  end if

  ! ── Read parameters ───────────────────────────────────────────────────────────

  call kh_io_read_params(trim(param_file), nx, ny, lx, ly, re, dt, &
                          nsteps, amp, mode, out_interval)

  ! Kinematic viscosity: ν = U₀ / Re  (re=0 → inviscid)
  if (re > 0.0_c_double) then
    nu = KH_U0_DEFAULT / re
  else
    nu = 0.0_c_double
  end if

  ! ── Open JSON output ──────────────────────────────────────────────────────────

  call kh_io_open_output(trim(out_file), out_unit)

  ! ── Print run header ──────────────────────────────────────────────────────────

  write(*, '(A)') "kh_main: KH instability pseudo-spectral reference solver"
  write(*, '(A,A)') "  Param file: ", trim(param_file)
  write(*, '(A,A)') "  Output:     ", trim(out_file)
  write(*, '(A,I0,A,I0)') "  Grid:  ", nx, " x ", ny
  write(*, '(A,F8.4,A,F8.4)') "  Domain: Lx = ", lx, "  Ly = ", ly
  write(*, '(A,ES12.4,A,ES12.4)') "  Re = ", re, "   nu = ", nu
  write(*, '(A,ES12.4,A,I0,A,ES10.4)') &
      "  dt = ", dt, "  steps = ", nsteps, &
      "  T_end = ", dt * real(nsteps, c_double)
  write(*, '(A,I0,A,I0)') "  mode = ", mode, "  out_interval = ", out_interval
  write(*, '(A)') "  Running ..."

  ! ── Run solver ────────────────────────────────────────────────────────────────

  call kh_solver_run(nx, ny, lx, ly, nu, dt, nsteps, amp, mode, &
                      out_interval, out_unit, &
                      ke_final, enstrophy_final, max_vort_final, div_rms_final, cfl_peak)

  ! ── Close output file ─────────────────────────────────────────────────────────

  call kh_io_close_output(out_unit)

  ! ── Print final diagnostics ───────────────────────────────────────────────────

  write(*, '(A)') "  Done."
  write(*, '(A)') "  --- Final diagnostics ---"
  write(*, '(A,ES14.6)') "  KE              = ", ke_final
  write(*, '(A,ES14.6)') "  Enstrophy       = ", enstrophy_final
  write(*, '(A,ES14.6)') "  max|ω|          = ", max_vort_final
  write(*, '(A,ES14.6)') "  div_rms         = ", div_rms_final
  write(*, '(A,ES14.6,A,ES12.4,A)') &
      "  Peak CFL        = ", cfl_peak, &
      "  (limit = ", KH_CFL_FACTOR, ")", &
      merge("  OK", "  !!", cfl_peak <= KH_CFL_FACTOR)

end program kh_main
