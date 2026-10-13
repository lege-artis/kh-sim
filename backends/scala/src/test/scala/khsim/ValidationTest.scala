package khsim

import io.circe.parser.*
import io.circe.{ACursor, Json}
import org.scalatest.funsuite.AnyFunSuite
import org.scalatest.matchers.should.Matchers

import java.nio.file.{Files, Paths}

/** KH-004 acceptance test — validates Scala diagnostics against the oracle
  * golden vector V1 (kh-sim/shared/golden-vectors/v1-canonical-reference.json,
  * class A2), which replaced the retired
  * kh-sim/shared/physics/kh_reference_output.json on 2026-10-08.
  *
  * The request is read from V1 "parameters"; the tolerance is V1's A2 band:
  * |got - ref| <= max(scalar_abs * f, scalar_rel * f * |ref|), f = the
  * "kh-scala" scalar_factor override (1 when absent).
  * Incompressibility: divergence_rms < 1e-10.
  *
  * Run: sbt test
  */
class ValidationTest extends AnyFunSuite with Matchers:

  private val Backend = "kh-scala"

  // ── Load golden vector V1 ───────────────────────────────────────────────────

  private def loadV1(): Json =
    // Path relative to sbt project root (kh-sim/backends/scala/)
    val path = Paths.get("../../shared/golden-vectors/v1-canonical-reference.json")
    assert(Files.exists(path), s"Golden vector not found: $path")
    parse(new String(Files.readAllBytes(path))) match
      case Right(json) => json
      case Left(e)     => throw new RuntimeException(s"Bad JSON: ${e.getMessage}")

  private def num(c: ACursor, key: String): Double =
    c.downField(key).as[Double] match
      case Right(v) => v
      case Left(e)  => throw new RuntimeException(s"golden vector: $key: ${e.getMessage}")

  private def int(c: ACursor, key: String): Int =
    c.downField(key).as[Int] match
      case Right(v) => v
      case Left(e)  => throw new RuntimeException(s"golden vector: $key: ${e.getMessage}")

  // ── Tests ────────────────────────────────────────────────────────────────────

  test("diagnostics within the A2 band of golden vector V1") {
    val v1   = loadV1().hcursor
    val p    = v1.downField("parameters")
    val rd   = v1.downField("diagnostics")
    val tols = v1.downField("tolerances")
    val factor = tols.downField("backend_overrides").downField(Backend)
      .downField("scalar_factor").as[Double].getOrElse(1.0)
    val sAbs = num(tols, "scalar_abs") * factor
    val sRel = num(tols, "scalar_rel") * factor

    val req = SimulationRequest(
      gridNx               = Some(int(p, "grid_nx")),
      gridNy               = Some(int(p, "grid_ny")),
      domainLx             = Some(num(p, "domain_lx")),
      domainLy             = Some(num(p, "domain_ly")),
      dt                   = Some(num(p, "dt")),
      steps                = Some(int(p, "nsteps")),
      reynoldsNumber       = Some(num(p, "reynolds_number")),
      velocityShear        = Some(num(p, "velocity_shear")),
      perturbationAmplitude= Some(num(p, "perturbation_amplitude")),
      perturbationMode     = Some(int(p, "perturbation_mode")),
    )

    val result = physics.Solver.simulate(req)
    val got    = result.diagnostics

    println(s"--- Scala vs golden vector V1 (A2 band, factor $factor) ---")
    val failures = Seq(
      "kinetic_energy" -> got.kineticEnergy,
      "enstrophy"      -> got.enstrophy,
      "max_vorticity"  -> got.maxVorticity,
    ).flatMap { (name, value) =>
      val ref  = num(rd, name)
      val diff = math.abs(value - ref)
      val band = math.max(sAbs, sRel * math.abs(ref))
      println(f"$name%-15s: got=$value%.17g  ref=$ref%.17g  |diff|=$diff%.3e  band=$band%.3e")
      if diff <= band then None else Some(f"$name: |diff| $diff%.3e > band $band%.3e")
    }
    println(f"divergence_rms : got=${got.divergenceRms}%.2e")

    failures shouldBe empty
    got.divergenceRms should be < 1e-10
  }

  test("output field shapes are correct") {
    val req = SimulationRequest(
      gridNx = Some(32), gridNy = Some(16),
      dt     = Some(0.001), steps = Some(5),
    )
    val result = physics.Solver.simulate(req)
    val n = 32 * 16
    result.uVelocity.length shouldBe n
    result.vVelocity.length shouldBe n
    result.vorticity.length shouldBe n
    result.pressure.length  shouldBe n
    result.tFinal            should be (0.005 +- 1e-9)
  }
