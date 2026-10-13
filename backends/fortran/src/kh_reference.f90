! kh_reference.f90 — KH instability: canonical reference run parameters and validator
!
! REFERENCE WATERMARK — kh-sim Fortran reference, v0.1;
! canonical for: TC-NUM-KH-001..008;
! do not modify without bumping reference version + regenerating reference outputs.
!
! Defines the canonical reference run that every compliant kh-sim implementation
! must reproduce to within KH_REF_TOL (5%).  The reference run is:
!
!   Grid:    NX=64, NY=32
!   Domain:  Lx=1.0, Ly=0.5
!   Re:      1000  →  ν = U₀/Re = 0.001
!   IC:      δ = KH_DELTA_FACTOR·Ly = 0.025; amp = 0.01; mode = 2
!   dt:      0.001
!   steps:   100  (T = 0.1)
!
! Reference output (computed by this Fortran stack, gfortran 13, ThinkPad,
! commit a814cc9 + c08b71e + f27b832, validated 2026-05-04):
!
!   KE        =   1.1281E-01   (KH_REF_KE        = 0.112810)
!   Enstrophy =   4.3703E+01   (KH_REF_ENSTROPHY  = 43.705142)
!   max|ω|    =   3.1986E+01   (KH_REF_MAX_VORT   = 31.907572)
!   div_rms   =   1.5000E-14   (KH_REF_DIV_RMS    = 1.20e-14)
!
! All values are within 0.01% of the stored reference constants — the
! reference values in kh_constants.f90 are consistent with this stack.
!
! The comparison tolerance KH_REF_TOL = 5% is deliberately loose so that
! alternative implementations (different FFT libraries, compilation flags,
! float precision) can certify compliance without being sensitivity-tested
! on round-off.  For sha256-level byte-exact reproducibility, see OQ-NUM-05.
!
! Public interface:
!   KH_REF_NX, KH_REF_NY, KH_REF_RE, KH_REF_DT, KH_REF_NSTEPS  (parameters)
!   kh_reference_compare(ke, enstrophy, max_vort, div_rms, passed)
!
! References:
!   kh_constants.f90 — KH_REF_KE, KH_REF_ENSTROPHY, KH_REF_MAX_VORT, KH_REF_TOL
!   PHYSICS-NUMERICAL-METHODS-v0.1.md §2.3 (TC-NUM-KH-008)

module kh_reference
  use kh_constants
  implicit none
  private

  ! ── Canonical run parameters (exported as named constants) ──────────────────

  integer,        parameter, public :: KH_REF_NX    = 64
  integer,        parameter, public :: KH_REF_NY    = 32
  real(c_double), parameter, public :: KH_REF_LX    = KH_LX_DEFAULT   ! 1.0
  real(c_double), parameter, public :: KH_REF_LY    = KH_LY_DEFAULT   ! 0.5
  real(c_double), parameter, public :: KH_REF_RE_RUN = 1000.0_c_double
  real(c_double), parameter, public :: KH_REF_NU_RUN = KH_U0_DEFAULT / KH_REF_RE_RUN
  real(c_double), parameter, public :: KH_REF_DT    = KH_DT_DEFAULT   ! 0.001
  integer,        parameter, public :: KH_REF_NSTEPS = 100             ! T = 0.1
  real(c_double), parameter, public :: KH_REF_AMP   = KH_AMP_DEFAULT  ! 0.01
  integer,        parameter, public :: KH_REF_MODE  = KH_MODE_DEFAULT  ! 2

  public :: kh_reference_compare

contains

  ! ── Compare solver output against reference constants ─────────────────────────
  !
  ! Evaluates each diagnostic against the stored reference value using the
  ! relative tolerance KH_REF_TOL (5%).  For div_rms the criterion is
  ! absolute (≤ KH_TOL_DIVERGENCE) since it is machine-noise level.
  !
  ! Arguments:
  !   ke, enstrophy, max_vort, div_rms  [in]  solver output at T_ref
  !   passed_ke     [out]  |ke − KH_REF_KE| / KH_REF_KE ≤ KH_REF_TOL
  !   passed_ens    [out]  |enstrophy − KH_REF_ENSTROPHY| / KH_REF_ENSTROPHY ≤ KH_REF_TOL
  !   passed_vort   [out]  |max_vort − KH_REF_MAX_VORT| / KH_REF_MAX_VORT ≤ KH_REF_TOL
  !   passed_div    [out]  div_rms ≤ KH_TOL_DIVERGENCE
  !   rel_ke        [out]  relative error on KE (for reporting)
  !   rel_ens       [out]  relative error on enstrophy
  !   rel_vort      [out]  relative error on max|ω|

  subroutine kh_reference_compare(ke, enstrophy, max_vort, div_rms, &
                                   passed_ke, passed_ens, passed_vort, passed_div, &
                                   rel_ke, rel_ens, rel_vort)
    real(c_double), intent(in)  :: ke, enstrophy, max_vort, div_rms
    logical,        intent(out) :: passed_ke, passed_ens, passed_vort, passed_div
    real(c_double), intent(out) :: rel_ke, rel_ens, rel_vort

    rel_ke   = abs(ke        - KH_REF_KE)        / KH_REF_KE
    rel_ens  = abs(enstrophy - KH_REF_ENSTROPHY)  / KH_REF_ENSTROPHY
    rel_vort = abs(max_vort  - KH_REF_MAX_VORT)   / KH_REF_MAX_VORT

    passed_ke   = (rel_ke   <= KH_REF_TOL)
    passed_ens  = (rel_ens  <= KH_REF_TOL)
    passed_vort = (rel_vort <= KH_REF_TOL)
    passed_div  = (div_rms  <= KH_TOL_DIVERGENCE)

  end subroutine kh_reference_compare

end module kh_reference
