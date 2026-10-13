package khsim

import io.circe.{Decoder, Encoder}
import io.circe.generic.semiauto.*

// ── Request ───────────────────────────────────────────────────────────────────

final case class SimulationRequest(
  gridNx:               Option[Int]    = None,
  gridNy:               Option[Int]    = None,
  domainLx:             Option[Double] = None,
  domainLy:             Option[Double] = None,
  dt:                   Option[Double] = None,
  steps:                Option[Int]    = None,
  reynoldsNumber:       Option[Double] = None,
  velocityShear:        Option[Double] = None,
  perturbationAmplitude:Option[Double] = None,
  perturbationMode:     Option[Int]    = None,
  initialOmega:         Option[Vector[Double]] = None,
  // ── v0.2 forcing API (phase-7f; PHASE-7F-FORCING-API-DESIGN-v0.1 §1+§7) ──
  // All default to None/off => exact pre-7f behaviour (no forcing, no drag),
  // so phase-6c golden vectors + phase-7d testbeds are bit-identically unaffected.
  forcingType:      Option[String] = None,   // "none" | "mode" | "band"
  forcingAmplitude: Option[Double] = None,    // A
  forcingKx:        Option[Int]    = None,    // "mode" only
  forcingKy:        Option[Int]    = None,    // "mode" only
  forcingKmin:      Option[Double] = None,    // "band" only
  forcingKmax:      Option[Double] = None,    // "band" only
  forcingOmega0:    Option[Double] = None,    // "band" only: phase-rotation rate scale (v0.2)
  linearDrag:       Option[Double] = None,    // Ekman linear drag alpha; -alpha*omega on RHS (all modes)
  // Time of the initial field (2026-10-10, CP-7F-B/D harness): forcing runs from
  // t_start and t_final = t_start + steps*dt, so a continued run keeps its phase.
  tStart:           Option[Double] = None,
):
  // Resolved values (openapi defaults)
  def nx:    Int    = gridNx.getOrElse(128)
  def ny:    Int    = gridNy.getOrElse(64)
  def lx:    Double = domainLx.getOrElse(1.0)
  def ly:    Double = domainLy.getOrElse(0.5)
  def dtVal: Double = dt.getOrElse(0.001)
  def nSteps:Int    = steps.getOrElse(100)
  def re:    Double = reynoldsNumber.getOrElse(1000.0)
  def u0:    Double = velocityShear.getOrElse(1.0)
  def amp:   Double = perturbationAmplitude.getOrElse(0.01)
  def mode:  Int    = perturbationMode.getOrElse(2)
  // Forcing-resolved (defaults => off; "none" path skips all forcing code)
  def forceType:   String = forcingType.getOrElse("none")
  def forceAmp:    Double  = forcingAmplitude.getOrElse(0.0)
  def forceKx:     Int     = forcingKx.getOrElse(1)
  def forceKy:     Int     = forcingKy.getOrElse(1)
  def forceKmin:   Double  = forcingKmin.getOrElse(0.0)
  def forceKmax:   Double  = forcingKmax.getOrElse(0.0)
  def forceOmega0: Double  = forcingOmega0.getOrElse(0.0)
  def drag:        Double  = linearDrag.getOrElse(0.0)
  def tStartVal:   Double  = tStart.getOrElse(0.0)

// Circe snake_case codec (matches openapi.yaml field names).
// initial_omega: optional flat row-major array (length nx*ny); when Some, overrides analytic IC.
// forcing_* / linear_drag: optional v0.2 forcing fields (absent => off; backward compatible).
object SimulationRequest:
  given Decoder[SimulationRequest] = Decoder.instance { c =>
    for
      nx    <- c.get[Option[Int]]("grid_nx")
      ny    <- c.get[Option[Int]]("grid_ny")
      lx    <- c.get[Option[Double]]("domain_lx")
      ly    <- c.get[Option[Double]]("domain_ly")
      dt    <- c.get[Option[Double]]("dt")
      steps <- c.get[Option[Int]]("steps")
      re    <- c.get[Option[Double]]("reynolds_number")
      u0    <- c.get[Option[Double]]("velocity_shear")
      amp   <- c.get[Option[Double]]("perturbation_amplitude")
      mode  <- c.get[Option[Int]]("perturbation_mode")
    yield SimulationRequest(
      nx, ny, lx, ly, dt, steps, re, u0, amp, mode,
      c.downField("initial_omega").as[Vector[Double]].toOption,
      c.downField("forcing_type").as[String].toOption,
      c.downField("forcing_amplitude").as[Double].toOption,
      c.downField("forcing_kx").as[Int].toOption,
      c.downField("forcing_ky").as[Int].toOption,
      c.downField("forcing_kmin").as[Double].toOption,
      c.downField("forcing_kmax").as[Double].toOption,
      c.downField("forcing_omega0").as[Double].toOption,
      c.downField("linear_drag").as[Double].toOption,
      c.downField("t_start").as[Double].toOption)
  }

// ── Response ──────────────────────────────────────────────────────────────────

final case class Diagnostics(
  kineticEnergy:  Double,
  enstrophy:      Double,
  maxVorticity:   Double,
  divergenceRms:  Double,
  cflPeak:        Double,   // max over all steps of dt*(max|u|/dx + max|v|/dy) — KH-NM-5
)
object Diagnostics:
  given Encoder[Diagnostics] = Encoder.forProduct5(
    "kinetic_energy", "enstrophy", "max_vorticity", "divergence_rms", "cfl_peak")(d =>
    (d.kineticEnergy, d.enstrophy, d.maxVorticity, d.divergenceRms, d.cflPeak))

final case class SimulationResult(
  backend:        String,
  language:       String,
  stepsCompleted: Int,
  tFinal:         Double,
  gridNx:         Int,
  gridNy:         Int,
  uVelocity:      Vector[Double],
  vVelocity:      Vector[Double],
  vorticity:      Vector[Double],
  pressure:       Vector[Double],
  diagnostics:    Diagnostics,
  computeTimeMs:  Double,
  warnings:       Vector[String] = Vector.empty,   // KH-VAL-2 K5 (filled by Routes)
)
object SimulationResult:
  given Encoder[SimulationResult] = Encoder.forProduct13(
    "backend", "language", "steps_completed", "t_final",
    "grid_nx", "grid_ny",
    "u_velocity", "v_velocity", "vorticity", "pressure",
    "diagnostics", "compute_time_ms", "warnings")(r =>
    (r.backend, r.language, r.stepsCompleted, r.tFinal,
     r.gridNx, r.gridNy,
     r.uVelocity, r.vVelocity, r.vorticity, r.pressure,
     r.diagnostics, r.computeTimeMs, r.warnings))

// ── Health / Info ─────────────────────────────────────────────────────────────

final case class HealthResponse(status: String, backend: String, port: Int)
object HealthResponse:
  given Encoder[HealthResponse] = deriveEncoder

final case class InfoResponse(
  backend:        String,
  language:       String,
  framework:      String,
  fftLibrary:     String,
  gridConstraint: String,   // "any length >= 2": JTransforms, axis of 1 refused
  buildDate:      String,   // YYYY-MM-DD, the packaged classes' timestamp
  port:           Int,
  openapiSpec:    String,
)
object InfoResponse:
  given Encoder[InfoResponse] = Encoder.forProduct8(
    "backend", "language", "framework", "fft_library", "grid_constraint", "build_date",
    "port", "openapi_spec")(r =>
    (r.backend, r.language, r.framework, r.fftLibrary, r.gridConstraint, r.buildDate,
     r.port, r.openapiSpec))
