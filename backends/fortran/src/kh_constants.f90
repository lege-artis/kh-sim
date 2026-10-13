! kh_constants.f90 — KH instability: physical and numerical constants
!
! REFERENCE WATERMARK — kh-sim Fortran reference, v0.1;
! canonical for: TC-NUM-KH-001..008;
! do not modify without bumping reference version + regenerating reference outputs.
!
! All constants derive from KH-PHYSICS.md §3 defaults and
! PHYSICS-NUMERICAL-METHODS-v0.1.md §2 (ETDRK4 stack spec).
!
! Usage:
!   use kh_constants
!   print *, KH_RE_DEFAULT, KH_LX_DEFAULT
!
! Reference: kh-sim/shared/physics/KH-PHYSICS.md
! Numerical methods: _config/PHYSICS-NUMERICAL-METHODS-v0.1.md §2

module kh_constants
  use iso_c_binding, only: c_double
  implicit none
  public

  ! ── Mathematical constants ───────────────────────────────────────────────────

  real(c_double), parameter :: KH_PI = &
      3.14159265358979323846_c_double

  real(c_double), parameter :: KH_TWO_PI = &
      6.28318530717958647692_c_double

  ! ── Default domain geometry (KH-PHYSICS §3) ─────────────────────────────────

  ! Domain length in x [non-dimensional]
  real(c_double), parameter :: KH_LX_DEFAULT = 1.0_c_double

  ! Domain length in y [non-dimensional]
  real(c_double), parameter :: KH_LY_DEFAULT = 0.5_c_double

  ! ── Default grid resolution ──────────────────────────────────────────────────

  ! Grid points in x (power of 2; default per KH-PHYSICS §4 and config)
  integer, parameter :: KH_NX_DEFAULT = 128

  ! Grid points in y (power of 2)
  integer, parameter :: KH_NY_DEFAULT = 64

  ! ── Default flow parameters (KH-PHYSICS §3) ─────────────────────────────────

  ! Shear velocity magnitude U₀ [non-dimensional]
  real(c_double), parameter :: KH_U0_DEFAULT = 1.0_c_double

  ! Reynolds number Re = U₀ · Lx / ν
  real(c_double), parameter :: KH_RE_DEFAULT = 1000.0_c_double

  ! Kinematic viscosity ν = U₀ / Re (derived; stored for reference)
  real(c_double), parameter :: KH_NU_DEFAULT = &
      KH_U0_DEFAULT / KH_RE_DEFAULT

  ! Density ratio ρ₂/ρ₁ (1.0 = homogeneous)
  real(c_double), parameter :: KH_DENSITY_RATIO_DEFAULT = 1.0_c_double

  ! ── Shear layer geometry ─────────────────────────────────────────────────────

  ! Shear layer thickness δ = 0.05 · Ly
  real(c_double), parameter :: KH_DELTA_FACTOR = 0.05_c_double

  ! ── Default perturbation (KH-PHYSICS §3) ────────────────────────────────────

  ! Initial perturbation amplitude A
  real(c_double), parameter :: KH_AMP_DEFAULT = 0.01_c_double

  ! Perturbation wavenumber mode k (integer)
  integer, parameter :: KH_MODE_DEFAULT = 2

  ! ── Default time integration ─────────────────────────────────────────────────

  ! Time step dt (conservative for Re=1000, 128×64 grid)
  real(c_double), parameter :: KH_DT_DEFAULT = 0.001_c_double

  ! Default number of time steps
  integer, parameter :: KH_STEPS_DEFAULT = 100

  ! ── ETDRK4 numerical method parameters (PHYSICS-NUMERICAL-METHODS §2.2) ─────

  ! CFL safety factor (used in adaptive dt: dt = KH_CFL_FACTOR * min(dx,dy) / |u|_max)
  real(c_double), parameter :: KH_CFL_FACTOR = 0.4_c_double

  ! 2/3-rule de-aliasing threshold factor (Orszag 1971)
  ! Modes |kx| > KH_DEALIAS_FACTOR * kx_max are zeroed
  real(c_double), parameter :: KH_DEALIAS_FACTOR = 2.0_c_double / 3.0_c_double

  ! ── Validation tolerances (PHYSICS-NUMERICAL-METHODS §1.5) ──────────────────

  ! Linear-stage tolerance: FFT round-trip, Poisson solve
  real(c_double), parameter :: KH_TOL_LINEAR = 1.0e-12_c_double

  ! Time-integrated nonlinear tolerance: field diagnostics at t_final
  real(c_double), parameter :: KH_TOL_NONLINEAR = 1.0e-6_c_double

  ! Incompressibility acceptance: divergence_rms < this
  real(c_double), parameter :: KH_TOL_DIVERGENCE = 1.0e-10_c_double

  ! ── Reference output diagnostics (KH-PHYSICS §7; 64×32 grid, 100 steps) ─────

  real(c_double), parameter :: KH_REF_KE       = 0.112810_c_double
  real(c_double), parameter :: KH_REF_ENSTROPHY = 43.705142_c_double
  real(c_double), parameter :: KH_REF_MAX_VORT  = 31.907572_c_double
  real(c_double), parameter :: KH_REF_DIV_RMS   = 1.20e-14_c_double

  ! Reference tolerance for integration-level validation (5% per KH-PHYSICS §7)
  real(c_double), parameter :: KH_REF_TOL = 0.05_c_double

end module kh_constants
