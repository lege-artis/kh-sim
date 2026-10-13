package khsim

import cats.effect.{IO, IOApp}
import com.comcast.ip4s.*
import org.http4s.ember.server.EmberServerBuilder
import org.http4s.server.middleware.{CORS, Logger as HttpLogger}
import org.typelevel.log4cats.slf4j.Slf4jLogger
import org.typelevel.log4cats.Logger

/** KH-SIM Scala/http4s backend — Task KH-004
  *
  * Binds on 0.0.0.0:8002.
  * Spec: kh-sim/shared/api/openapi.yaml
  */
object Main extends IOApp.Simple:

  def run: IO[Unit] =
    given logger: Logger[IO] = Slf4jLogger.getLogger[IO]

    val routes = CORS.policy
      .withAllowOriginAll
      .withAllowHeadersAll
      .withAllowMethodsAll
      .apply(Routes())

    val loggedRoutes = HttpLogger.httpRoutes(
      logHeaders = false,
      logBody    = false,
    )(routes)

    EmberServerBuilder
      .default[IO]
      .withHost(ipv4"0.0.0.0")
      .withPort(port"8002")
      .withHttpApp(loggedRoutes.orNotFound)
      .build
      .use: server =>
        logger.info(s"kh-sim Scala backend listening on http://0.0.0.0:8002") >>
        IO.never
