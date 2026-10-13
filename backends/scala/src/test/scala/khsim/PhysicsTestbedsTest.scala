package khsim

/** PhysicsTestbedsTest.scala — phase 7d physics testbeds, Scala backend.
  *
  * Mode 1: HTTP POST /simulate to running Scala server on port 8002.
  * BACKEND_NAME = "kh-scala"; backend_factor = 2.0 (JTransforms FFT
  * accumulated rounding, calibrated phase 6c).
  *
  * Tier 1: TB-1a (Taylor-Green exact), TB-1b (pure diffusion mode),
  *         TB-1c (sigma(k) linear-stability dispersion)
  * Tier 2: TB-2a (total vorticity conservation), TB-2b (inviscid KE),
  *         TB-2c (inviscid enstrophy), TB-2d (viscous energy decay DC-corrected),
  *         TB-2e (incompressibility)
  * Tier 3: TB-3a (spatial convergence TG), TB-3b (temporal convergence RK4)
  * PENDING: TB-3c (MMS, PENDING-OPUS-7c), TB-4a (Lyapunov, PENDING-OPUS-7e),
  *          TB-4b (invariant measure, STRETCH-DEFER)
  *
  * Run:
  *   cd backends/scala
  *   sbt run &
  *   sleep 10  # JVM cold start
  *   sbt "testOnly khsim.PhysicsTestbedsTest"
  *   kill %1
  */

import io.circe.parser.*
import io.circe.*
import org.scalatest.funsuite.AnyFunSuite
import org.scalatest.matchers.should.Matchers

import java.net.URI
import java.net.http.{HttpClient, HttpRequest, HttpResponse}
import scala.util.Try

class PhysicsTestbedsTest extends AnyFunSuite with Matchers:

  val BackendUrl: String = "http://localhost:8002/simulate"

  // ── Server availability ────────────────────────────────────────────────────

  private val serverAvailable: Boolean =
    Try {
      val client = HttpClient.newHttpClient()
      val req = HttpRequest.newBuilder()
        .uri(URI.create("http://localhost:8002/health"))
        .GET().build()
      val resp = client.send(req, HttpResponse.BodyHandlers.ofString())
      resp.statusCode() == 200
    }.getOrElse(false)

  private def checkServer(): Unit =
    if !serverAvailable then
      cancel("Scala server not available on port 8002 — start with `sbt run` first")

  // ── HTTP helper ────────────────────────────────────────────────────────────

  private def postSimulate(body: String): Json =
    val client = HttpClient.newHttpClient()
    val req = HttpRequest.newBuilder()
      .uri(URI.create(BackendUrl))
      .header("Content-Type", "application/json")
      .POST(HttpRequest.BodyPublishers.ofString(body))
      .build()
    Try {
      val resp = client.send(req, HttpResponse.BodyHandlers.ofString())
      if resp.statusCode() / 100 != 2 then
        fail(s"HTTP ${resp.statusCode()}: ${resp.body().take(300)}")
      parse(resp.body()).getOrElse(fail(s"JSON parse error: ${resp.body().take(200)}"))
    }.recover {
      case e: java.net.ConnectException =>
        cancel(s"Scala server not reachable at $BackendUrl: ${e.getMessage}")
    }.get

  // ── Request builder ────────────────────────────────────────────────────────

  private def mkReq(
    nx: Int, ny: Int,
    lx: Double, ly: Double,
    dt: Double, steps: Int,
    re: Double, shear: Double,
    amp: Double, mode: Int,
    initOmega: Option[Array[Double]] = None,
  ): String =
    val base =
      s"""{"grid_nx":$nx,"grid_ny":$ny,"domain_lx":$lx,"domain_ly":$ly,""" +
      s""""dt":$dt,"steps":$steps,"reynolds_number":$re,""" +
      s""""velocity_shear":$shear,"perturbation_amplitude":$amp,"perturbation_mode":$mode"""
    initOmega match
      case Some(arr) => base + s""","initial_omega":[${arr.mkString(",")}]}"""
      case None      => base + "}"

  // ── JSON accessors ─────────────────────────────────────────────────────────

  private def diagDbl(json: Json, field: String): Double =
    json.hcursor.downField("diagnostics").downField(field).as[Double]
      .getOrElse(fail(s"missing diagnostics.$field in response"))

  private def vortArr(json: Json): Array[Double] =
    json.hcursor.downField("vorticity").as[Vector[Double]]
      .getOrElse(fail("missing vorticity in response"))
      .toArray

  // ── IC builders ────────────────────────────────────────────────────────────

  /** Taylor-Green IC: omega_0(x,y) = -2*cos(x)*cos(y) on [0,2pi]^2. */
  private def tgIC(nx: Int, ny: Int, lx: Double, ly: Double): Array[Double] =
    val arr = new Array[Double](nx * ny)
    val dx  = lx / nx
    val dy  = ly / ny
    for i <- 0 until nx; j <- 0 until ny do
      arr(i * ny + j) = -2.0 * math.cos(i * dx) * math.cos(j * dy)
    arr

  /** Single-mode diffusion IC: omega_0(x,y) = sin(kx*x)*sin(ky*y) on [0,2pi]^2. */
  private def sineIC(nx: Int, ny: Int, lx: Double, ly: Double, kx: Int, ky: Int): Array[Double] =
    val arr = new Array[Double](nx * ny)
    val dx  = lx / nx
    val dy  = ly / ny
    for i <- 0 until nx; j <- 0 until ny do
      arr(i * ny + j) = math.sin(kx * i * dx) * math.sin(ky * j * dy)
    arr

  /** Canonical KH IC: tanh shear layer + cosine perturbation.
    * Mirrors Solver.scala initialConditions exactly (row-major, x-outer). */
  private def canonicalIC(
    nx: Int, ny: Int, lx: Double, ly: Double,
    u0: Double, delta: Double, amp: Double, mode: Int,
  ): Array[Double] =
    val arr = new Array[Double](nx * ny)
    val dx  = lx / nx
    val dy  = ly / ny
    for i <- 0 until nx do
      val x     = i * dx
      val pertX = amp * (2.0 * math.Pi * mode / lx) * math.cos(2.0 * math.Pi * mode * x / lx)
      for j <- 0 until ny do
        val y   = j * dy
        val z   = (y - ly / 2.0) / delta
        val shr = u0 / delta / (math.cosh(z) * math.cosh(z))
        arr(i * ny + j) = shr + pertX
    arr

  // ── Array helpers ──────────────────────────────────────────────────────────

  private def rms(a: Array[Double]): Double =
    math.sqrt(a.map(x => x * x).sum / a.length)

  private def mean(a: Array[Double]): Double =
    a.sum / a.length

  // ── Tier 1: Exact-solution tests ───────────────────────────────────────────

  test("TB-1a: Taylor-Green vortex") {
    checkServer()
    // Exact solution: omega(x,y,t) = -2*cos(x)*cos(y)*exp(-2*nu*t)
    // Nonlinear term vanishes identically for TG on [0,2pi]^2 (KH-PHYSICS §2).
    // Pass: rms_rel < 2e-6 (2× factor; Python ref: 6.7e-16).
    val nx = 64; val ny = 64
    val lx = 2.0 * math.Pi; val ly = 2.0 * math.Pi
    val nu = 0.01
    val T  = 1.0
    val dt = 0.001
    val nsteps = (T / dt).toInt  // 1000

    val omega0     = tgIC(nx, ny, lx, ly)
    val omegaExact = omega0.map(_ * math.exp(-2.0 * nu * T))
    val rmsExact   = rms(omegaExact)

    val body = mkReq(nx, ny, lx, ly, dt, nsteps,
      re = 1.0 / nu, shear = 1.0, amp = 0.0, mode = 0,
      initOmega = Some(omega0))
    val resp     = postSimulate(body)
    val omegaSim = vortArr(resp)

    val diff   = omegaSim.zip(omegaExact).map { case (s, e) => s - e }
    val rmsRel = rms(diff) / rmsExact

    assert(rmsRel < 2e-6,
      f"TB-1a Taylor-Green rms_rel=$rmsRel%.3e >= 2e-6 " +
      f"(rms_diff=${rms(diff)}%.3e, rms_exact=$rmsExact%.3e)")
  }

  test("TB-1b: pure diffusion single mode") {
    checkServer()
    // Single-mode IC omega_0 = sin(x)*sin(y) on [0,2pi]^2.
    // Mode amplitude should decay as exp(-nu*(kx^2+ky^2)*t).
    // Nonlinear advection is non-zero for this IC (unlike TG) but small.
    // Tolerance 2e-4 accounts for nonlinear correction (Python ref: 0.0).
    val nx = 32; val ny = 32
    val lx = 2.0 * math.Pi; val ly = 2.0 * math.Pi
    val kxTest = 1; val kyTest = 1
    val nu     = 0.1
    val dt     = 0.001
    val nsteps = 500
    val T      = dt * nsteps

    val omega0        = sineIC(nx, ny, lx, ly, kxTest, kyTest)
    val k2            = (kxTest * kxTest + kyTest * kyTest).toDouble
    val expectedDecay = math.exp(-nu * k2 * T)

    val body = mkReq(nx, ny, lx, ly, dt, nsteps,
      re = 1.0 / nu, shear = 1.0, amp = 0.0, mode = 0,
      initOmega = Some(omega0))
    val resp     = postSimulate(body)
    val omegaSim = vortArr(resp)

    // Mode projection: sum(omega * sin(x)*sin(y)); avoids FFT in test harness.
    // For single-mode IC, proj(omega_0) = Nx*Ny/4 (sum of sin^2 terms).
    // ratio = proj(omega_sim) / proj(omega_0) approximates amplitude decay.
    val dx = lx / nx; val dy = ly / ny
    var proj0 = 0.0; var projT = 0.0
    for i <- 0 until nx; j <- 0 until ny do
      val basis = math.sin(kxTest * i * dx) * math.sin(kyTest * j * dy)
      proj0 += omega0(i * ny + j) * basis
      projT += omegaSim(i * ny + j) * basis

    assert(math.abs(proj0) > 1e-10, f"TB-1b: initial mode projection too small: $proj0%.3e")
    val ratio  = projT / proj0
    val relErr = math.abs(ratio - expectedDecay) / expectedDecay

    // 2× factor: 1e-4 * 2 = 2e-4
    assert(relErr < 2e-4,
      f"TB-1b pure diffusion: ratio=$ratio%.8f expected=$expectedDecay%.8f " +
      f"rel_err=$relErr%.3e >= 2e-4 (T=$T, nu=$nu, k2=$k2)")
  }

  test("TB-1c: sigma(k) linear-stability dispersion") {
    checkServer()
    // Two-run subtraction method (single-run amplitude fails: base-flow omega
    // at kx=2 dominates perturbation).  See STATUS-python.md §6.
    // sigma_expected = 0.18 * U0/delta = 7.2 /s (KH-LS-2, inviscid theory).
    // Pass: rel_err < 20% (2× factor; Python ref: 8.0%).
    val nx = 128; val ny = 64
    val lx = 1.0; val ly = 0.5
    val dt        = 0.001
    val nstepsMeas = 50      // T_meas = 0.05 s
    val ampPert   = 1e-4
    val mode      = 2
    val re        = 1000.0; val shear = 1.0
    val delta     = 0.05 * ly
    val sigmaExpected = 0.18 * shear / delta  // = 7.2 /s

    // t=0 ICs for energy_pert_0
    val omegaBase0 = canonicalIC(nx, ny, lx, ly, shear, delta, 0.0, 0)
    val omegaPert0 = canonicalIC(nx, ny, lx, ly, shear, delta, ampPert, mode)
    val diff0      = omegaBase0.zip(omegaPert0).map { case (b, p) => p - b }
    val ePert0     = 0.5 * diff0.map(x => x * x).sum / (nx * ny)

    if ePert0 < 1e-20 then
      cancel(f"TB-1c: perturbation energy at t=0 too small: $ePert0%.2e")

    // Base run (amp=0): no perturbation
    val bodyBase  = mkReq(nx, ny, lx, ly, dt, nstepsMeas, re, shear, 0.0, 0)
    val omegaBaseT = vortArr(postSimulate(bodyBase))

    // Perturbed run
    val bodyPert  = mkReq(nx, ny, lx, ly, dt, nstepsMeas, re, shear, ampPert, mode)
    val omegaPertT = vortArr(postSimulate(bodyPert))

    val diffT  = omegaBaseT.zip(omegaPertT).map { case (b, p) => p - b }
    val ePertT = 0.5 * diffT.map(x => x * x).sum / (nx * ny)

    if ePertT < 1e-20 then
      cancel(f"TB-1c: perturbation energy at T too small: $ePertT%.2e")

    val tMeas    = dt * nstepsMeas
    val sigmaSim = math.log(ePertT / ePert0) / (2.0 * tMeas)
    val relErr   = math.abs(sigmaSim - sigmaExpected) / sigmaExpected

    // 2× factor: 10% * 2 = 20%
    assert(relErr < 0.20,
      f"TB-1c sigma(k): sigma_sim=$sigmaSim%.4f /s expected=$sigmaExpected%.4f /s " +
      f"rel_err=$relErr%.3e >= 20%% " +
      f"(E_pert_0=$ePert0%.3e E_pert_T=$ePertT%.3e T_meas=$tMeas%.3f)")
  }

  // ── Tier 2: Conservation tests ─────────────────────────────────────────────

  test("TB-2a: total vorticity conservation") {
    checkServer()
    // Spectral method conserves DC mode (sum(omega)) exactly.
    // Canonical KH IC has sum(omega) ~ 8192 (positive shear layer; sech^2 > 0).
    // Test: |sum_T - sum_0| < 2e-8 (2× factor; Python ref: 0.0 exactly).
    val nx = 64; val ny = 32
    val lx = 1.0; val ly = 0.5
    val delta = 0.05 * ly

    val omega0 = canonicalIC(nx, ny, lx, ly, 1.0, delta, 0.01, 2)
    val sum0   = omega0.sum

    val body = mkReq(nx, ny, lx, ly, 0.001, 100, 1000.0, 1.0, 0.01, 2)
    val sumT = vortArr(postSimulate(body)).sum

    val drift = math.abs(sumT - sum0)
    assert(drift < 2e-8,
      f"TB-2a vorticity conservation: drift=$drift%.3e >= 2e-8 " +
      f"(sum_0=$sum0%.6e, sum_T=$sumT%.6e)")
  }

  test("TB-2b: inviscid kinetic energy conservation") {
    checkServer()
    // Near-inviscid (nu=1e-8): |E(T)-E(0)| / E(0) < 2e-4 (2× factor).
    // Python ref: 7.2e-7. Tolerance accounts for 2/3-rule dealiasing energy leak.
    // E(0) approximated via nsteps=1 (E(dt) ~ E(0) to O(dt^4) for RK4).
    val nx = 64; val ny = 32
    val lx = 1.0; val ly = 0.5
    val nu    = 1e-8; val shear = 1.0
    val re    = shear / nu

    val body1 = mkReq(nx, ny, lx, ly, 0.0005, 1, re, shear, 0.01, 2)
    val E0    = diagDbl(postSimulate(body1), "kinetic_energy")

    val body  = mkReq(nx, ny, lx, ly, 0.0005, 200, re, shear, 0.01, 2)
    val ET    = diagDbl(postSimulate(body), "kinetic_energy")

    if E0 < 1e-15 then cancel(f"TB-2b: E0 too small ($E0%.2e)")
    val relErr = math.abs(ET - E0) / E0

    // 2× factor: 1e-4 * 2 = 2e-4
    assert(relErr < 2e-4,
      f"TB-2b inviscid KE: E0=$E0%.8e ET=$ET%.8e rel_err=$relErr%.3e >= 2e-4")
  }

  test("TB-2c: inviscid enstrophy conservation") {
    checkServer()
    // Near-inviscid (nu=1e-8): |Z(T)-Z(0)| / Z(0) < 2e-4 (2× factor).
    // Python ref: 1.2e-6. Same rationale as TB-2b.
    val nx = 64; val ny = 32
    val lx = 1.0; val ly = 0.5
    val nu    = 1e-8; val shear = 1.0
    val re    = shear / nu

    val body1 = mkReq(nx, ny, lx, ly, 0.0005, 1, re, shear, 0.01, 2)
    val Z0    = diagDbl(postSimulate(body1), "enstrophy")

    val body  = mkReq(nx, ny, lx, ly, 0.0005, 200, re, shear, 0.01, 2)
    val ZT    = diagDbl(postSimulate(body), "enstrophy")

    if Z0 < 1e-15 then cancel(f"TB-2c: Z0 too small ($Z0%.2e)")
    val relErr = math.abs(ZT - Z0) / Z0

    // 2× factor: 1e-4 * 2 = 2e-4
    assert(relErr < 2e-4,
      f"TB-2c inviscid enstrophy: Z0=$Z0%.8e ZT=$ZT%.8e rel_err=$relErr%.3e >= 2e-4")
  }

  test("TB-2d: viscous energy decay rate (DC-corrected)") {
    checkServer()
    // Energy balance: dE/dt = -2*nu*Z_ac, where Z_ac = Z - Z_dc.
    // Z_dc = 0.5 * mean(omega)^2 — DC mode excluded (psi_hat[0,0]=0 in solver;
    // DC vorticity has no kinetic energy and is not dissipated by viscosity).
    // Pass: rel_err < 10% (2× factor; Python ref: 0.31% with DC correction).
    val nx = 64; val ny = 32
    val lx = 1.0; val ly = 0.5
    val dt = 0.001; val re = 1000.0; val shear = 1.0
    val nu = shear / re

    val bodyT1 = mkReq(nx, ny, lx, ly, dt, 100, re, shear, 0.01, 2)
    val bodyT2 = mkReq(nx, ny, lx, ly, dt, 200, re, shear, 0.01, 2)
    val bodyTm = mkReq(nx, ny, lx, ly, dt, 150, re, shear, 0.01, 2)

    val respT1 = postSimulate(bodyT1)
    val respT2 = postSimulate(bodyT2)
    val respTm = postSimulate(bodyTm)

    val E1     = diagDbl(respT1, "kinetic_energy")
    val E2     = diagDbl(respT2, "kinetic_energy")
    val ZTotal = diagDbl(respTm, "enstrophy")

    val omegaMid = vortArr(respTm)
    val meanOm   = mean(omegaMid)
    val Z_dc     = 0.5 * meanOm * meanOm
    val Z_ac     = ZTotal - Z_dc

    val T1 = dt * 100; val T2 = dt * 200
    val dEdtFd     = (E2 - E1) / (T2 - T1)
    val dEdtTheory = -2.0 * nu * Z_ac
    val denom      = math.abs(dEdtTheory)

    if denom < 1e-15 then cancel(f"TB-2d: |2*nu*Z_ac| too small ($denom%.2e)")
    val relErr = math.abs(dEdtFd - dEdtTheory) / denom

    // 2× factor: 5% * 2 = 10%
    assert(relErr < 0.10,
      f"TB-2d energy decay: dE/dt_fd=$dEdtFd%.6e theory=$dEdtTheory%.6e " +
      f"Z_total=$ZTotal%.6e Z_dc=$Z_dc%.6e Z_ac=$Z_ac%.6e rel_err=$relErr%.3e >= 10%%")
  }

  test("TB-2e: incompressibility") {
    checkServer()
    // Spectral vorticity-streamfunction formulation enforces div(u)=0 by construction.
    // Pass: divergence_rms < 2e-12 (2× factor; Python ref: 1.1e-14).
    val body   = mkReq(64, 32, 1.0, 0.5, 0.001, 100, 1000.0, 1.0, 0.01, 2)
    val divRms = diagDbl(postSimulate(body), "divergence_rms")

    // 2× factor: 1e-12 * 2 = 2e-12
    assert(divRms < 2e-12,
      f"TB-2e incompressibility: divergence_rms=$divRms%.3e >= 2e-12")
  }

  // ── Tier 3: Convergence tests ──────────────────────────────────────────────

  test("TB-3a: spatial convergence Taylor-Green") {
    checkServer()
    // TG vortex on [0,2pi]^2 at 4 grids: 32^2, 64^2, 128^2, 256^2.
    // T=0.1, nu=0.01, dt=0.0001 (1000 steps each).
    // Expected: all errors at machine-eps floor (spectral convergence
    // achieved at 32^2; TG uses only (1,1) Fourier mode).
    // Pass: (a) max_rms_err < 2e-12; (b) errors monotone non-increasing.
    // Same floor-test as Python (OQ-7D-PY-3 resolution); 2× factor on bound.
    val nu     = 0.01
    val T      = 0.1
    val dt     = 0.0001
    val nsteps = (T / dt).toInt  // 1000
    val grids  = List(32, 64, 128, 256)
    val errors = collection.mutable.ArrayBuffer[Double]()

    for N <- grids do
      val lx = 2.0 * math.Pi; val ly = 2.0 * math.Pi
      val omega0     = tgIC(N, N, lx, ly)
      val omegaExact = omega0.map(_ * math.exp(-2.0 * nu * T))
      val body       = mkReq(N, N, lx, ly, dt, nsteps,
        re = 1.0 / nu, shear = 1.0, amp = 0.0, mode = 0,
        initOmega = Some(omega0))
      val omegaSim = vortArr(postSimulate(body))
      val diff     = omegaSim.zip(omegaExact).map { case (s, e) => s - e }
      errors += rms(diff)

    val maxErr = errors.max
    // 2× factor: 1e-12 * 2 = 2e-12
    assert(maxErr < 2e-12,
      f"TB-3a spatial convergence: max_rms_err=$maxErr%.3e >= 2e-12 " +
      s"(grids: ${grids.zip(errors.toList).map { case (n, e) => f"$n^2: $e%.3e" }.mkString(", ")})")

    // Monotone non-increasing (10% slack for machine-eps noise)
    for i <- 0 until errors.length - 1 do
      assert(errors(i + 1) <= errors(i) * 1.1,
        f"TB-3a: error increased from ${grids(i)}^2 (${errors(i)}%.3e) " +
        f"to ${grids(i+1)}^2 (${errors(i+1)}%.3e) by more than 10%%")
  }

  test("TB-3b: temporal convergence RK4") {
    checkServer()
    // 4 dt values: 1e-2, 5e-3, 2e-3, 1e-3; reference at dt_ref=5e-4.
    // T_final = 0.1 for all runs; nu=0.001 (Re=1000) for stable dt range.
    // Fit log(error) vs log(dt); assert slope >= 3.5 (RK4 theoretical: 4).
    // Python ref: slope = 4.11. Not tolerance-scaled (order test).
    val nx = 64; val ny = 32
    val lx = 1.0; val ly = 0.5
    val nu = 0.001; val shear = 1.0
    val re = shear / nu
    val T  = 0.1

    val bodyRef   = mkReq(nx, ny, lx, ly, 5e-4, 200, re, shear, 0.01, 2)
    val omegaRef  = vortArr(postSimulate(bodyRef))

    val dtValues = List(1e-2, 5e-3, 2e-3, 1e-3)
    val errors   = collection.mutable.ArrayBuffer[Double]()
    for dtTest <- dtValues do
      val ns      = (T / dtTest).round.toInt
      val bodyT   = mkReq(nx, ny, lx, ly, dtTest, ns, re, shear, 0.01, 2)
      val omegaT  = vortArr(postSimulate(bodyT))
      val diff    = omegaT.zip(omegaRef).map { case (a, b) => a - b }
      errors += rms(diff)

    // Log-log linear regression for convergence order
    val logDt  = dtValues.map(math.log).toArray
    val logErr = errors.map(math.log).toArray
    val n      = logDt.length.toDouble
    val sumX   = logDt.sum; val sumY = logErr.sum
    val sumXX  = logDt.map(x => x * x).sum
    val sumXY  = logDt.zip(logErr).map { case (x, y) => x * y }.sum
    val slope  = (n * sumXY - sumX * sumY) / (n * sumXX - sumX * sumX)

    assert(slope >= 3.5,
      f"TB-3b temporal convergence RK4: slope=$slope%.3f < 3.5 " +
      s"(dt_values=$dtValues, errors=${errors.map(e => f"$e%.3e").mkString(",")})")
  }

  // ── Tier 3 PENDING ─────────────────────────────────────────────────────────

  test("TB-3c: method of manufactured solutions") {
    cancel("PENDING-OPUS-7c: awaits Opus MMS design (manufactured solution + source term)")
  }

  // ── Tier 4 PENDING / DEFER ─────────────────────────────────────────────────

  test("TB-4a: Lyapunov time estimate") {
    cancel("PENDING-OPUS-7e: awaits Opus Lyapunov estimator design (T_L auto-detection algorithm)")
  }

  test("TB-4b: invariant measure agreement") {
    cancel("STRETCH-DEFER per MASTER-BRIEF §2.5 (long-time ensemble; high complexity)")
  }
