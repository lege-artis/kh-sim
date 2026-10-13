! test_num_002_poisson.f90 — TC-NUM-KH-002: spectral Poisson solver correctness
!
! TC-NUM-KH-002 specification (PHYSICS-NUMERICAL-METHODS-v0.1.md §2.3):
!   Test target:  Spectral Poisson solver
!   Method check: Solve ∇²ψ = -ω where ω is a single Fourier mode cos(2πk·x/Lx);
!                 compare numerical ψ against analytical solution.
!   Pass criterion: rel err ≤ 1e-12
!
! Analytical derivation:
!   ω(x, y) = cos(2π·k·x / Lx)   (x-only mode, no y-dependence)
!
!   In Fourier space, mode (±k, 0) carries the content:
!     ω̂(±k, 0) = Nx·Ny / 2    (DFT of real cosine, with our normalisation)
!
!   Poisson solve: ψ̂ = ω̂ / k²,  where k² = kx(k,0)² + ky(k,0)² = kx_k²
!     kx_k = 2π·k / Lx
!
!   Inverse FFT → ψ(x, y) = cos(2π·k·x / Lx) / kx_k²
!
!   Relative error metric:
!     rel_err = max |ψ_numerical(i,j) - ψ_analytical(i,j)| / max|ψ_analytical|
!
! Exit code: 0 = PASS, 1 = FAIL
!
! Compilation (DRY-RUN — validate on ThinkPad):
!   gfortran -O2 -o test_num_002 \
!       src/kh_constants.f90 src/kh_grid.f90 src/kh_fft.f90 \
!       src/kh_poisson.f90 src/kh_velocity.f90 \
!       tests/test_num_002_poisson.f90
!   ./test_num_002

program test_num_002_poisson
  use kh_constants
  use kh_grid,    only: kh_grid_make, kh_grid_wavenums
  use kh_fft,     only: kh_fft_forward_2d, kh_fft_inverse_2d
  use kh_poisson, only: kh_poisson_solve
  implicit none

  ! ── Test parameters ──────────────────────────────────────────────────────────

  integer,        parameter :: NX   = 64
  integer,        parameter :: NY   = 32
  real(c_double), parameter :: LX   = KH_LX_DEFAULT   ! 1.0
  real(c_double), parameter :: LY   = KH_LY_DEFAULT   ! 0.5
  integer,        parameter :: KMODE = 3               ! perturbation wavenumber

  ! ── Local arrays ─────────────────────────────────────────────────────────────

  real(c_double) :: x(0:NX-1), y(0:NY-1), dx, dy
  real(c_double) :: kx(0:NX-1, 0:NY-1), ky(0:NX-1, 0:NY-1), k2(0:NX-1, 0:NY-1)
  real(c_double) :: omega(0:NX-1, 0:NY-1)
  real(c_double) :: psi_numerical(0:NX-1, 0:NY-1)
  real(c_double) :: psi_analytical(0:NX-1, 0:NY-1)
  complex(8)     :: omega_hat(0:NX-1, 0:NY-1)
  complex(8)     :: psi_hat(0:NX-1, 0:NY-1)
  complex(8)     :: psi_work(0:NX-1, 0:NY-1)

  real(c_double) :: kx_mode, psi_scale
  real(c_double) :: abs_err, rel_err, psi_max
  integer        :: i, j
  logical        :: passed

  ! ── Build grid ───────────────────────────────────────────────────────────────

  call kh_grid_make(NX, NY, LX, LY, x, y, dx, dy)
  call kh_grid_wavenums(NX, NY, dx, dy, kx, ky, k2)

  ! ── Construct test vorticity field ───────────────────────────────────────────
  ! ω(x, y) = cos(2π·KMODE·x / Lx)

  kx_mode = KH_TWO_PI * real(KMODE, c_double) / LX

  do j = 0, NY-1
    do i = 0, NX-1
      omega(i, j) = cos(kx_mode * x(i))
    end do
  end do

  ! ── Analytical streamfunction ─────────────────────────────────────────────────
  ! ψ_exact(x, y) = cos(2π·KMODE·x / Lx) / kx_mode²

  psi_scale = 1.0_c_double / (kx_mode * kx_mode)
  do j = 0, NY-1
    do i = 0, NX-1
      psi_analytical(i, j) = cos(kx_mode * x(i)) * psi_scale
    end do
  end do

  ! ── Numerical Poisson solve ───────────────────────────────────────────────────

  ! Forward FFT of ω
  omega_hat = cmplx(omega, 0.0_8, kind=8)
  call kh_fft_forward_2d(omega_hat, NX, NY)

  ! Solve ψ̂ = ω̂ / k²  (psi_hat(0,0) zeroed internally)
  call kh_poisson_solve(omega_hat, k2, NX, NY, psi_hat)

  ! Inverse FFT → ψ in physical space
  psi_work = psi_hat
  call kh_fft_inverse_2d(psi_work, NX, NY)
  psi_numerical = real(psi_work, kind=8)

  ! ── Compute relative error ────────────────────────────────────────────────────

  psi_max = maxval(abs(psi_analytical))
  abs_err = 0.0_c_double
  do j = 0, NY-1
    do i = 0, NX-1
      abs_err = max(abs_err, abs(psi_numerical(i,j) - psi_analytical(i,j)))
    end do
  end do
  rel_err = abs_err / psi_max

  ! ── Evaluate pass / fail ─────────────────────────────────────────────────────

  passed = (rel_err <= KH_TOL_LINEAR)

  write(*, '(A)') "TC-NUM-KH-002: Spectral Poisson solver"
  write(*, '(A,I0,A,I0,A)') "  Grid: ", NX, " x ", NY
  write(*, '(A,I0,A,F8.4)') "  Mode k = ", KMODE, ",  kx_mode = ", kx_mode
  write(*, '(A,ES12.4)') "  max |ψ_analytical|           = ", psi_max
  write(*, '(A,ES12.4)') "  max |ψ_numerical - ψ_exact|  = ", abs_err
  write(*, '(A,ES12.4)') "  Relative error               = ", rel_err
  write(*, '(A,ES12.4)') "  Tolerance (KH_TOL_LINEAR)    = ", KH_TOL_LINEAR

  if (passed) then
    write(*, '(A)') "  RESULT: PASS"
    stop 0
  else
    write(*, '(A)') "  RESULT: FAIL"
    stop 1
  end if

end program test_num_002_poisson
