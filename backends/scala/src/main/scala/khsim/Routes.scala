package khsim

import cats.effect.IO
import io.circe.Json
import io.circe.syntax.*
import org.http4s.*
import org.http4s.circe.*
import org.http4s.dsl.io.*
import org.typelevel.log4cats.Logger

import physics.Solver

object Routes:
  private val BackendId = "kh-scala"

  private def requestError(bad: KhInvalid): IO[Response[IO]] =
    BadRequest(Json.obj("error" -> Json.fromString(bad.message), "field" -> Json.fromString(bad.field)))

  private def backendError(message: String): IO[Response[IO]] =
    InternalServerError(Json.obj("error" -> Json.fromString(message), "backend" -> Json.fromString(BackendId)))

  /** The date the backend's classes were packaged (UTC, YYYY-MM-DD). In the image
    * that is the build date, which lets a deploy tell a fresh build from a stale
    * one; "unknown" if it cannot be read.
    */
  private val buildDate: String =
    scala.util.Try {
      val loc = classOf[InfoResponse].getProtectionDomain.getCodeSource.getLocation.toURI
      val ms  = java.nio.file.Files.getLastModifiedTime(java.nio.file.Paths.get(loc)).toMillis
      java.time.Instant.ofEpochMilli(ms).atZone(java.time.ZoneOffset.UTC).toLocalDate.toString
    }.getOrElse("unknown")

  def apply()(using logger: Logger[IO]): HttpRoutes[IO] =
    HttpRoutes.of[IO]:

      // -- POST /simulate -----------------------------------------------------
      // 2026-10-03 (finding v0.5, S-1..S-3): the body is validated by Contract
      // (KH-VAL-1) before the solver runs, and every failure answers in the
      // openapi.yaml error shape instead of an empty 500.

      case req @ POST -> Root / "simulate" =>
        req.as[String].flatMap { text =>
          Contract.parse(text) match
            case Left(bad) => requestError(bad)
            case Right(parsed) =>
              // forcing_amplitude 0 means no forcing (kh_physics.py); do not let a
              // zero-amplitude band request reach buildBand's square-domain check.
              val zeroBand = parsed.forceType == "band" && parsed.forceAmp == 0.0
              val body     = if zeroBand then parsed.copy(forcingType = Some("none")) else parsed
              val run =
                logger.info(s"simulate nx=${body.nx} ny=${body.ny} steps=${body.nSteps} Re=${body.re}") >>
                IO.blocking(Solver.simulate(body))   // off the CE thread-pool
              run.attempt.flatMap {
                case Left(e) =>
                  logger.warn(s"simulate failed: $e") >>
                  backendError(s"computation failed: ${e.getClass.getSimpleName}: ${e.getMessage}")
                case Right(result) if !Contract.isFinite(result) =>
                  backendError("computation produced a non-finite value (dt beyond the stability bound?)")
                case Right(result) =>
                  logger.info(
                    s"done t=${result.tFinal} ke=${result.diagnostics.kineticEnergy} " +
                    s"ms=${result.computeTimeMs}") >>
                  Ok(result.copy(warnings = Contract.warnings(body, result.diagnostics.cflPeak)).asJson)
              }
        }

      // -- GET /health --------------------------------------------------------

      case GET -> Root / "health" =>
        Ok(HealthResponse("ok", "scala-http4s", 8002).asJson)

      // -- GET /info ----------------------------------------------------------

      case GET -> Root / "info" =>
        Ok(InfoResponse(
          backend        = "scala-http4s",
          language       = "Scala",
          framework      = "http4s 0.23 (Ember)",
          fftLibrary     = "JTransforms 3.1",
          gridConstraint = "any length >= 2",
          buildDate      = buildDate,
          port           = 8002,
          openapiSpec    = "kh-sim/shared/api/openapi.yaml",
        ).asJson)
