// KH-SIM — Scala backend (KH-004)
// Framework: http4s 0.23 (Ember server) + cats-effect 3
// FFT:       JTransforms 3.1 (DoubleFFT_2D)
// JSON:      circe

val http4sVersion   = "0.23.29"
val circeVersion    = "0.14.10"
val catsEffectVer   = "3.5.7"
val jTransformsVer  = "3.1"
val scalaTestVer    = "3.2.19"

lazy val root = (project in file("."))
  .enablePlugins(JavaAppPackaging)
  .settings(
    name             := "kh-sim-scala",
    version          := "0.1.0",
    scalaVersion     := "3.3.4",
    organization     := "org.zemla.khsim",
    mainClass        := Some("khsim.Main"),

    scalacOptions ++= Seq(
      "-deprecation",
      "-feature",
      "-unchecked",
      "-Wunused:all",
    ),

    libraryDependencies ++= Seq(
      // HTTP server (cats-effect 3 based)
      "org.http4s"     %% "http4s-ember-server" % http4sVersion,
      "org.http4s"     %% "http4s-circe"        % http4sVersion,
      "org.http4s"     %% "http4s-dsl"          % http4sVersion,

      // JSON codec
      "io.circe"       %% "circe-generic"       % circeVersion,
      "io.circe"       %% "circe-parser"        % circeVersion,

      // Functional effect system
      "org.typelevel"  %% "cats-effect"         % catsEffectVer,

      // 2D FFT (Java library, JVM-accessible from Scala)
      "com.github.wendykierp" % "JTransforms"   % jTransformsVer,

      // Logging
      "org.typelevel"  %% "log4cats-slf4j"      % "2.7.0",
      "ch.qos.logback"  % "logback-classic"      % "1.5.13",

      // Test
      "org.scalatest"  %% "scalatest"           % scalaTestVer % Test,
      "io.circe"       %% "circe-parser"        % circeVersion % Test,
    ),
  )
