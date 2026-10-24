package khsim

import io.circe.{Json, JsonObject}

/** A request the backend refuses: HTTP 400 {"error": message, "field": field}. */
final case class KhInvalid(field: String, message: String)

/** POST /simulate request contract KH-VAL-1 (Scala implementation).
  *
  * Same rules as backends/cpp/src/kh_http_contract.hpp and the Pascal server
  * (finding v0.5, development record, finding 2026-09-16):
  *   - body must be a JSON object                          -> 400 field "body"
  *   - a key present with the wrong JSON type               -> 400 field <key>
  *   - grid_nx, grid_ny in [2, 65536]: JTransforms takes any length, so there is
  *     no power-of-two rule, but DoubleFFT_2D refuses an axis of 1, which the
  *     other backends accept; grid_nx*grid_ny <= 4194304   -> 400 field grid_nx|grid_ny
  *   - domain_lx, domain_ly, reynolds_number non-zero      -> 400 field <key>
  *   - initial_omega (legacy alias init_omega, not both)
  *     has exactly grid_nx*grid_ny entries                 -> 400 field initial_omega
  *   - forcing_type in {none, mode, band}; band with forcing_amplitude != 0
  *     needs |domain_lx - domain_ly| < 1e-12 (the test Solver.buildBand applies;
  *     C++ uses a relative 1e-10)                          -> 400 field forcing_type
  * Until 2026-10-03 an unknown forcing_type was silently computed as "none", and
  * a malformed body or an impossible grid produced a 500 with an empty body.
  *
  * Contract KH-VAL-2 (development record, modelling rules v0.1, ruling 2026-10-07),
  * after the KH-VAL-1 checks, same rules and order as backends/cpp/src/kh_http_contract.hpp:
  *   K3 outside the model: domain_lx, domain_ly, reynolds_number < 0; dt < 0; steps < 0;
  *      linear_drag < 0; perturbation_mode < 0; band forcing with kmin > kmax -> 400 field
  *   K4 not representable (analytic field, amplitude != 0): perturbation_mode 0 or
  *      beyond the 2/3-rule band (3*mode > 2*(grid_nx/2))                     -> 400 field
  *   K4b (ruling 2026-10-09) a non-zero perturbation_amplitude sent with
  *      initial_omega (it would be silently ignored)                         -> 400 field
  *   K4 forcing (amplitude != 0, 2026-10-10): mode forcing_kx/ky beyond the 2/3-rule
  *      band, or a band with 3*floor(kmax) > 2*(min(nx, ny)/2)                -> 400 field
  *   K3 a band with no integer mode kmin <= |k| <= kmax                  -> 400 forcing_kmax
  *   K5 warnings in the response (warnings).
  *   K7 (2026-10-24) a key the contract does not define, checked first; with several,
  *      the first in byte order; foreign names get a hint                  -> 400 field
  * t_start (2026-10-10): the time of the initial field (forcing phase, t_final).
  */
object Contract:
  val MaxAxis: Long  = 65536L
  val MaxCells: Long = 4194304L

  private def typeName(j: Json): String =
    j.fold("null", _ => "boolean", _ => "number", _ => "string", _ => "array", _ => "object")

  private def present(o: JsonObject, key: String): Option[Json] =
    o(key).filterNot(_.isNull)

  private def check(ok: Boolean, field: String, message: => String): Either[KhInvalid, Unit] =
    if ok then Right(()) else Left(KhInvalid(field, message))

  private def readInt(o: JsonObject, key: String): Either[KhInvalid, Option[Int]] =
    present(o, key) match
      case None => Right(None)
      case Some(j) =>
        j.asNumber match
          case None => Left(KhInvalid(key, s"$key must be an integer, got ${typeName(j)}"))
          case Some(n) =>
            n.toLong match
              case Some(l) if l >= -2147483647L && l <= 2147483647L => Right(Some(l.toInt))
              case Some(_) => Left(KhInvalid(key, s"$key is outside the 32-bit integer range"))
              case None    => Left(KhInvalid(key, s"$key must be an integer"))

  private def readDouble(o: JsonObject, key: String): Either[KhInvalid, Option[Double]] =
    present(o, key) match
      case None => Right(None)
      case Some(j) =>
        j.asNumber match
          case Some(n) => Right(Some(n.toDouble))
          case None    => Left(KhInvalid(key, s"$key must be a number, got ${typeName(j)}"))

  private def readString(o: JsonObject, key: String): Either[KhInvalid, Option[String]] =
    present(o, key) match
      case None => Right(None)
      case Some(j) =>
        j.asString match
          case Some(s) => Right(Some(s))
          case None    => Left(KhInvalid(key, s"$key must be a string, got ${typeName(j)}"))

  private def readField(o: JsonObject, key: String): Either[KhInvalid, Option[Vector[Double]]] =
    present(o, key) match
      case None => Right(None)
      case Some(j) =>
        j.asArray match
          case None => Left(KhInvalid("initial_omega", s"$key must be an array of numbers"))
          case Some(arr) if arr.length.toLong > MaxCells =>
            Left(KhInvalid("initial_omega", s"$key has more entries than the $MaxCells-cell limit"))
          case Some(arr) =>
            val nums = arr.flatMap(_.asNumber).map(_.toDouble)
            if nums.length == arr.length then Right(Some(nums))
            else Left(KhInvalid("initial_omega", s"$key must contain numbers only"))

  // DoubleFFT_2D(rows, columns) throws for an axis below 2, so 1 is refused here.
  private def checkAxis(key: String, n: Int): Either[KhInvalid, Unit] =
    check(n >= 2 && n <= MaxAxis, key,
      s"$key = $n is out of range [2, $MaxAxis]; JTransforms needs at least 2 points per axis")

  private def validate(r: SimulationRequest): Either[KhInvalid, Unit] =
    val cells = r.nx.toLong * r.ny.toLong
    val omegaOk = r.initialOmega.forall(_.length.toLong == cells)
    val omegaLen = r.initialOmega.map(_.length).getOrElse(0)
    val forcingCheck: Either[KhInvalid, Unit] = r.forceType match
      case "none" | "mode" => Right(())
      case "band" =>
        // Same test as Solver.Forcing.buildBand's require, so a request the
        // solver would abort on is refused here instead.
        check(r.forceAmp == 0.0 || math.abs(r.lx - r.ly) < 1e-12,
          "forcing_type", "band forcing requires domain_lx == domain_ly")
      case other =>
        Left(KhInvalid("forcing_type", s"unknown forcing_type '$other'; expected none, mode or band"))
    for
      _ <- checkAxis("grid_nx", r.nx)
      _ <- checkAxis("grid_ny", r.ny)
      _ <- check(cells <= MaxCells, "grid_nx",
             s"grid ${r.nx}x${r.ny} has $cells cells, more than the $MaxCells-cell limit")
      _ <- check(r.lx != 0.0, "domain_lx", "domain_lx must be non-zero")
      _ <- check(r.ly != 0.0, "domain_ly", "domain_ly must be non-zero")
      _ <- check(r.re != 0.0, "reynolds_number", "reynolds_number must be non-zero")
      _ <- check(omegaOk, "initial_omega",
             s"initial_omega has $omegaLen entries; grid ${r.nx}x${r.ny} needs $cells")
      _ <- forcingCheck
      // KH-VAL-2, K3: outside the model.
      _ <- check(r.lx >= 0.0, "domain_lx", "domain_lx must be positive")
      _ <- check(r.ly >= 0.0, "domain_ly", "domain_ly must be positive")
      _ <- check(r.re >= 0.0, "reynolds_number",
             "reynolds_number must be positive (a negative value is negative viscosity)")
      _ <- check(r.dtVal >= 0.0, "dt", "dt must be >= 0 (backward integration of a viscous flow is ill-posed)")
      _ <- check(r.nSteps >= 0, "steps", "steps must be >= 0")
      _ <- check(r.drag >= 0.0, "linear_drag", "linear_drag must be >= 0 (negative drag injects energy)")
      _ <- check(r.mode >= 0, "perturbation_mode", "perturbation_mode must be >= 0")
      _ <- check(!(r.forceType == "band" && r.forceAmp != 0.0 && r.forceKmin > r.forceKmax),
             "forcing_kmin", "forcing_kmin > forcing_kmax: the forcing band is empty")
      // KH-VAL-2, K4 (forcing): the forced modes must lie inside the 2/3-rule band.
      _ <- check(!(r.forceType == "mode" && r.forceAmp != 0.0 && 3L * math.abs(r.forceKx.toLong) > 2L * (r.nx / 2)),
             "forcing_kx", s"forcing_kx ${r.forceKx} lies beyond the 2/3-rule band of grid_nx ${r.nx}")
      _ <- check(!(r.forceType == "mode" && r.forceAmp != 0.0 && 3L * math.abs(r.forceKy.toLong) > 2L * (r.ny / 2)),
             "forcing_ky", s"forcing_ky ${r.forceKy} lies beyond the 2/3-rule band of grid_ny ${r.ny}")
      _ <- check(!(r.forceType == "band" && r.forceAmp != 0.0
                   && 3.0 * math.floor(r.forceKmax) > 2.0 * (math.min(r.nx, r.ny) / 2)),
             "forcing_kmax", s"forcing_kmax reaches beyond the 2/3-rule band of the ${r.nx}x${r.ny} grid")
      // K3: the band must contain an integer mode (half-plane, kmin <= |k| <= kmax).
      _ <- check(!(r.forceType == "band" && r.forceAmp != 0.0 && !bandHasMode(r.forceKmin, r.forceKmax)),
             "forcing_kmax", "the forcing band contains no integer mode: nothing would be forced")
      // KH-VAL-2, K4: the analytic perturbation must exist on the grid, inside the 2/3-rule band.
      _ <- check(r.initialOmega.isDefined || r.amp == 0.0 || r.mode != 0, "perturbation_mode",
             "perturbation_mode 0 with a non-zero amplitude perturbs nothing")
      _ <- check(r.initialOmega.isDefined || r.amp == 0.0 || 3L * r.mode <= 2L * (r.nx / 2),
             "perturbation_mode",
             s"perturbation_mode ${r.mode} lies beyond the 2/3-rule band of grid_nx ${r.nx} " +
             s"(largest resolved mode ${2 * (r.nx / 2) / 3})")
      // KH-VAL-2, K4b: initial_omega replaces the analytic field; an amplitude sent with it would be ignored.
      _ <- check(!(r.initialOmega.isDefined && r.perturbationAmplitude.exists(_ != 0.0)),
             "perturbation_amplitude", "perturbation_amplitude has no effect with initial_omega (send 0 or omit it)")
    yield ()

  /** A half-plane integer mode (ix, iy) with kmin <= |k| <= kmax exists. Called only
    * after the 2/3-rule check, so kmax is bounded by the grid. */
  private def bandHasMode(kmin: Double, kmax: Double): Boolean =
    val km = math.floor(kmax).toInt
    (0 to km).exists(ix => (-km to km).exists { iy =>
      !(ix == 0 && iy <= 0) && {
        val kr = math.sqrt((ix * ix + iy * iy).toDouble)
        kr >= kmin && kr <= kmax
      }
    })

  private val Known: Set[String] = Set(
    "grid_nx", "grid_ny", "domain_lx", "domain_ly", "dt", "steps", "reynolds_number",
    "velocity_shear", "perturbation_amplitude", "perturbation_mode", "initial_omega",
    "init_omega", "forcing_type", "forcing_amplitude", "forcing_kx", "forcing_ky",
    "forcing_kmin", "forcing_kmax", "forcing_omega0", "linear_drag", "t_start")

  private val Foreign: Map[String, String] = Map(
    "nsteps" -> "steps", "n_steps" -> "steps", "num_steps" -> "steps",
    "nx" -> "grid_nx", "ny" -> "grid_ny", "lx" -> "domain_lx",
    "ly" -> "domain_ly", "re" -> "reynolds_number", "Re" -> "reynolds_number",
    "reynolds" -> "reynolds_number", "amplitude" -> "perturbation_amplitude", "mode" -> "perturbation_mode")

  /** KH-VAL-2, K7 (step 2, 2026-10-24): a key the contract does not define is refused
    * before any value is read; with several, the first in byte order is named. */
  private def checkUnknown(o: JsonObject): Either[KhInvalid, Unit] =
    o.keys.filterNot(Known.contains).minOption match
      case None => Right(())
      case Some(u) =>
        System.err.println(s"KH-VAL-2 K7: unknown key '$u' refused")
        val msg = Foreign.get(u) match
          case Some(to) => s"unknown key '$u' (did you mean '$to'?)"
          case None     => s"unknown key '$u' (not part of the request contract)"
        Left(KhInvalid(u, msg))

  /** KH-VAL-2, K5: computed, but possibly not resolved. Same conditions, order and
    * prefixes in every backend and the oracle (kh_physics.request_warnings). */
  def warnings(r: SimulationRequest, cflPeak: Double): Vector[String] =
    val nu   = math.abs(r.u0) / r.re
    val dmin = math.min(r.lx / r.nx, r.ly / r.ny)
    val diff = 2.0 * nu * r.dtVal / (dmin * dmin)
    Vector(
      Option.when(cflPeak > 1.0)(f"cfl_advective: peak Courant number $cflPeak%.3g exceeds 1"),
      Option.when(diff > 1.0)(f"cfl_diffusive: 2*nu*dt/min(dx,dy)^2 = $diff%.3g exceeds 1"),
      Option.when(nu == 0.0)("inviscid: nu = 0, no physical dissipation; enstrophy accumulates at the 2/3 cutoff"),
    ).flatten

  /** Parse and validate a POST /simulate body. */
  def parse(body: String): Either[KhInvalid, SimulationRequest] =
    for
      json  <- io.circe.parser.parse(body).left.map(_ => KhInvalid("body", "request body is not valid JSON"))
      o     <- json.asObject.toRight(KhInvalid("body", "request body must be a JSON object"))
      _     <- checkUnknown(o)
      nx    <- readInt(o, "grid_nx")
      ny    <- readInt(o, "grid_ny")
      lx    <- readDouble(o, "domain_lx")
      ly    <- readDouble(o, "domain_ly")
      dt    <- readDouble(o, "dt")
      steps <- readInt(o, "steps")
      re    <- readDouble(o, "reynolds_number")
      u0    <- readDouble(o, "velocity_shear")
      amp   <- readDouble(o, "perturbation_amplitude")
      mode  <- readInt(o, "perturbation_mode")
      ftype <- readString(o, "forcing_type")
      famp  <- readDouble(o, "forcing_amplitude")
      fkx   <- readInt(o, "forcing_kx")
      fky   <- readInt(o, "forcing_ky")
      fkmin <- readDouble(o, "forcing_kmin")
      fkmax <- readDouble(o, "forcing_kmax")
      fom0  <- readDouble(o, "forcing_omega0")
      drag  <- readDouble(o, "linear_drag")
      tst   <- readDouble(o, "t_start")
      canonical <- readField(o, "initial_omega")
      legacy    <- readField(o, "init_omega")
      _     <- check(canonical.isEmpty || legacy.isEmpty, "initial_omega",
                 "send initial_omega or the legacy init_omega, not both")
      req = SimulationRequest(
              nx, ny, lx, ly, dt, steps, re, u0, amp, mode,
              canonical.orElse(legacy),
              ftype, famp, fkx, fky, fkmin, fkmax, fom0, drag, tst)
      _     <- validate(req)
    yield req

  /** True when every number in the result is finite; a blown-up run is a 500. */
  def isFinite(r: SimulationResult): Boolean =
    val d = r.diagnostics
    val scalars = Seq(d.kineticEnergy, d.enstrophy, d.maxVorticity, d.divergenceRms, d.cflPeak)
    (scalars.iterator ++ r.uVelocity.iterator ++ r.vVelocity.iterator ++
      r.vorticity.iterator ++ r.pressure.iterator).forall(x => java.lang.Double.isFinite(x))
