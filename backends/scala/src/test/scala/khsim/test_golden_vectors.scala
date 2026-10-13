package khsim

/** test_golden_vectors.scala — per-backend golden-vector regression tests, Scala backend.
  *
  * Mode 1: HTTP POST /simulate to running server on port 8002.
  * Compares response against oracle diagnostics + field_summaries using
  * Higham §4.2 hybrid tolerance with backend_factor = 2.0
  * ("kh-scala" override in each vector JSON; JTransforms butterfly order
  * differs from NumPy Cooley-Tukey → ~2x accumulated rounding).
  *
  * Run (requires server on :8002):
  *   cd backends/scala
  *   sbt run &
  *   sleep 5
  *   sbt test
  *   kill %1
  *
  * Each test covers all three check types for that vector:
  *   1. Scalar diagnostics (KE, enstrophy, max_vort, div_rms; cfl_peak if present)
  *   2. Field summary statistics (min/max/mean per field u, v, ω, ψ)
  *   3. Field element comparison (first_10 / last_10 per field)
  *
  * If the server is unreachable the test is cancelled (not failed).
  * cfl_peak is validated only if present in the response (backward-compatible
  * with pre-Stage-0 server build).
  *
  * No new sbt dependencies required:
  *   - ScalaTest 3.2.19 (already in build.sbt % Test)
  *   - circe-parser (already in build.sbt main + % Test scope)
  *   - java.net.http.HttpClient (JDK 11+, no dep needed)
  */

import io.circe.parser.*
import io.circe.*
import org.scalatest.funsuite.AnyFunSuite
import org.scalatest.matchers.should.Matchers

import java.net.URI
import java.net.http.{HttpClient, HttpRequest, HttpResponse}
import java.nio.file.{Files, Paths}
import scala.jdk.CollectionConverters.*
import scala.util.Try

class GoldenVectorsTest extends AnyFunSuite with Matchers:

  // ── Configuration ──────────────────────────────────────────────────────────

  /** Backend name as it appears in vector JSON tolerances.backend_overrides */
  val BackendName: String = "kh-scala"
  val BackendUrl:  String = "http://localhost:8002/simulate"

  /** Path relative to sbt project root (backends/scala/) */
  val GoldenDir = Paths.get("../../shared/golden-vectors")

  val MachineEps: Double = 2.22e-16   // IEEE 754 binary64

  // ── Vector loading ─────────────────────────────────────────────────────────

  private def loadVectors(): List[(String, Json)] =
    val stream = Files.newDirectoryStream(GoldenDir, "v*.json")
    try
      stream.iterator().asScala.toList
        .sortBy(_.getFileName.toString)
        .map { p =>
          val content = new String(Files.readAllBytes(p))
          (p.getFileName.toString, parse(content).getOrElse(Json.Null))
        }
    finally
      stream.close()

  // ── Tolerance helpers (Higham §4.2 hybrid) ────────────────────────────────

  private def backendFactor(tols: Json, kind: String): Double =
    tols.hcursor
      .downField("backend_overrides")
      .downField(BackendName)
      .downField(s"${kind}_factor")
      .as[Double]
      .getOrElse(1.0)

  private def scalarTol(oracle: Double, absT: Double, relT: Double, factor: Double): Double =
    math.max(absT * factor, relT * factor * math.abs(oracle))

  private def fieldTol(oracle: Double, absT: Double, relT: Double,
                       factor: Double, sqrtN: Double, useSqrtN: Boolean): Double =
    val epsTerm = if useSqrtN then sqrtN * MachineEps * math.abs(oracle) else 0.0
    math.max(absT * factor, relT * factor * math.abs(oracle) + epsTerm)

  // ── HTTP helper ────────────────────────────────────────────────────────────

  /** POST request body JSON to BackendUrl. Returns Left("SKIP:...") on connection error. */
  private def postSimulate(reqBody: String): Either[String, Json] =
    val httpClient = HttpClient.newHttpClient()
    val req = HttpRequest.newBuilder()
      .uri(URI.create(BackendUrl))
      .header("Content-Type", "application/json")
      .POST(HttpRequest.BodyPublishers.ofString(reqBody))
      .build()
    Try {
      val resp = httpClient.send(req, HttpResponse.BodyHandlers.ofString())
      if resp.statusCode() / 100 != 2 then
        Left(s"HTTP ${resp.statusCode()}: ${resp.body().take(200)}")
      else
        parse(resp.body()).left.map(_.message)
    }.recover {
      case e: java.net.ConnectException =>
        Left(s"SKIP: ${e.getMessage}")
      case e: java.io.IOException if e.getMessage != null && e.getMessage.contains("Connection refused") =>
        Left(s"SKIP: ${e.getMessage}")
    }.toEither.left.map(_.getMessage).flatten

  // ── Request builder ────────────────────────────────────────────────────────

  /** Build REST request body from vector parameters; maps "nsteps" → "steps". */
  private def paramsToRequest(params: Json): String =
    val c = params.hcursor
    def raw(field: String): String =
      c.downField(field).focus.map(_.noSpaces).getOrElse("null")
    // REST endpoint uses "steps"; golden-vector JSON uses "nsteps"
    s"""{
       |  "grid_nx": ${raw("grid_nx")},
       |  "grid_ny": ${raw("grid_ny")},
       |  "domain_lx": ${raw("domain_lx")},
       |  "domain_ly": ${raw("domain_ly")},
       |  "dt": ${raw("dt")},
       |  "steps": ${raw("nsteps")},
       |  "reynolds_number": ${raw("reynolds_number")},
       |  "velocity_shear": ${raw("velocity_shear")},
       |  "perturbation_amplitude": ${raw("perturbation_amplitude")},
       |  "perturbation_mode": ${raw("perturbation_mode")}
       |}""".stripMargin

  // ── Core comparison ────────────────────────────────────────────────────────

  private def compareVector(filename: String, body: Json, vector: Json): Seq[String] =
    val failures = collection.mutable.ArrayBuffer[String]()

    val tols     = vector.hcursor.downField("tolerances").focus.getOrElse(Json.Null)
    val oracle   = vector.hcursor.downField("diagnostics").focus.getOrElse(Json.Null)
    val oracleFw = vector.hcursor.downField("field_summaries").focus.getOrElse(Json.Null)
    val params   = vector.hcursor.downField("parameters").focus.getOrElse(Json.Null)

    val nx       = params.hcursor.downField("grid_nx").as[Int].getOrElse(64)
    val ny       = params.hcursor.downField("grid_ny").as[Int].getOrElse(32)
    val sqrtN    = math.sqrt((nx * ny).toDouble)
    val useSqrtN = tols.hcursor.downField("sqrt_n_scaling").as[Boolean].getOrElse(false)
    val sf       = backendFactor(tols, "scalar")
    val ff       = backendFactor(tols, "field")

    def dbl(j: Json, field: String): Option[Double] =
      j.hcursor.downField(field).as[Double].toOption

    val sAbs  = dbl(tols, "scalar_abs").getOrElse(1e-12)
    val sRel  = dbl(tols, "scalar_rel").getOrElse(1e-9)
    val sdAbs = dbl(tols, "scalar_div_rms_abs").getOrElse(1e-11)
    val sdRel = dbl(tols, "scalar_div_rms_rel").getOrElse(1e-8)
    val fsAbs = dbl(tols, "field_summary_abs").getOrElse(1e-10)
    val fsRel = dbl(tols, "field_summary_rel").getOrElse(1e-7)
    val feAbs = dbl(tols, "field_element_abs").getOrElse(1e-9)
    val feRel = dbl(tols, "field_element_rel").getOrElse(1e-7)

    val bDiag = body.hcursor.downField("diagnostics").focus.getOrElse(body)

    // ── 1. Scalar diagnostics ─────────────────────────────────────────────────
    val scalarChecks = List(
      ("kinetic_energy", sAbs,  sRel),
      ("enstrophy",      sAbs,  sRel),
      ("max_vorticity",  sAbs,  sRel),
      ("divergence_rms", sdAbs, sdRel),
    )
    for (key, absT, relT) <- scalarChecks do
      for
        bv <- dbl(bDiag, key)
        ov <- dbl(oracle, key)
      do
        val tol  = scalarTol(ov, absT, relT, sf)
        val diff = math.abs(bv - ov)
        if diff >= tol then
          failures += f"$key: |diff|=$diff%.3e >= tol=$tol%.3e (actual=$bv%.6e, oracle=$ov%.6e)"

    // cfl_peak: optional — skip gracefully if backend does not emit it yet
    for
      bv <- dbl(bDiag, "cfl_peak")
      ov <- dbl(oracle, "cfl_peak")
    do
      val tol  = scalarTol(ov, sAbs, sRel, sf)
      val diff = math.abs(bv - ov)
      if diff >= tol then
        failures += f"cfl_peak: |diff|=$diff%.3e >= tol=$tol%.3e (actual=$bv%.6e, oracle=$ov%.6e)"

    // ── 2. Field summaries + 3. Field elements ───────────────────────────────
    val fieldMapping = List(
      ("u_velocity", "u_velocity"),
      ("v_velocity", "v_velocity"),
      ("vorticity",  "vorticity"),
      ("pressure",   "pressure"),
    )
    for (oracleKey, bodyKey) <- fieldMapping do
      val oFs = oracleFw.hcursor.downField(oracleKey).focus.getOrElse(Json.Null)
      body.hcursor.downField(bodyKey).as[Vector[Double]].toOption match
        case None =>
          failures += s"$bodyKey: missing or not a float array in response"
        case Some(flat) if flat.isEmpty =>
          failures += s"$bodyKey: empty array"
        case Some(flat) =>
          // Summary stats
          val bMin  = flat.min
          val bMax  = flat.max
          val bMean = flat.sum / flat.size
          for (stat, bv, oKey) <- List(("min", bMin, "min"), ("max", bMax, "max"), ("mean", bMean, "mean")) do
            dbl(oFs, oKey).foreach { ov =>
              val tol  = fieldTol(ov, fsAbs, fsRel, ff, sqrtN, useSqrtN)
              val diff = math.abs(bv - ov)
              if diff >= tol then
                failures += f"$oracleKey.$stat: |diff|=$diff%.3e >= tol=$tol%.3e (actual=$bv%.6e, oracle=$ov%.6e)"
            }
          // Element-wise: first_10 / last_10
          val oFirst = oFs.hcursor.downField("first_10").as[Vector[Double]].getOrElse(Vector.empty)
          val oLast  = oFs.hcursor.downField("last_10").as[Vector[Double]].getOrElse(Vector.empty)
          val bFirst = flat.take(oFirst.size)
          val bLast  = flat.takeRight(oLast.size)
          for (slot, bElems, oElems) <- List(("first_10", bFirst, oFirst), ("last_10", bLast, oLast)) do
            for ((bv, ov), i) <- bElems.zip(oElems).zipWithIndex do
              val tol  = fieldTol(ov, feAbs, feRel, ff, sqrtN, useSqrtN)
              val diff = math.abs(bv - ov)
              if diff >= tol then
                failures += f"$oracleKey.$slot[$i]: |diff|=$diff%.3e >= tol=$tol%.3e (actual=$bv%.6e, oracle=$ov%.6e)"

    failures.toSeq

  // ── Per-vector test runner ─────────────────────────────────────────────────

  private def runVectorTest(filename: String): Unit =
    val vectors = loadVectors()
    val (_, vector) = vectors
      .find(_._1 == filename)
      .getOrElse(fail(s"golden vector file not found in $GoldenDir: $filename"))
    val params  = vector.hcursor.downField("parameters").focus.getOrElse(Json.Null)
    val reqBody = paramsToRequest(params)
    postSimulate(reqBody) match
      case Left(msg) if msg.startsWith("SKIP:") =>
        cancel(s"Scala server not reachable at $BackendUrl — $msg")
      case Left(msg) =>
        fail(s"$filename: HTTP error: $msg")
      case Right(body) =>
        val failures = compareVector(filename, body, vector)
        if failures.nonEmpty then
          fail(
            s"$filename: ${failures.size} comparison failure(s):\n" +
            failures.map(f => s"  $f").mkString("\n")
          )

  // ── 9 parametrised test cases (one per golden vector) ─────────────────────

  test("V1-canonical-reference")  { runVectorTest("v1-canonical-reference.json")  }
  test("V2-low-re-diffusive")     { runVectorTest("v2-low-re-diffusive.json")     }
  test("V3-high-re-nonlinear")    { runVectorTest("v3-high-re-nonlinear.json")    }
  test("V4-double-resolution")    { runVectorTest("v4-double-resolution.json")    }
  test("V5-half-resolution-fast") { runVectorTest("v5-half-resolution-fast.json") }
  test("V6-long-time-evolution")  { runVectorTest("v6-long-time-evolution.json")  }
  test("V7-alt-mode-k1")          { runVectorTest("v7-alt-mode-k1.json")          }
  test("V8-dealiasing-edge-k4")   { runVectorTest("v8-dealiasing-edge-k4.json")   }
  test("V9-zero-perturbation")    { runVectorTest("v9-zero-perturbation.json")    }
