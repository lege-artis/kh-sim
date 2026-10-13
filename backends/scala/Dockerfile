# KH-SIM Scala backend — KH-004
# Port: 8002  |  Main: khsim.Main  |  Framework: http4s 0.23 / cats-effect 3
# Build tool: sbt 1.x + sbt-native-packager (staged deployment)
#
# Build:  docker build -t kh-sim-scala .
# Run:    docker run -p 8002:8002 kh-sim-scala
# Health: curl http://localhost:8002/health

# ── Stage 1: compile + stage ──────────────────────────────────────────────────
FROM eclipse-temurin:17-jdk-jammy AS builder

# Install sbt
RUN apt-get update && apt-get install -y --no-install-recommends curl gnupg && \
    echo "deb https://repo.scala-sbt.org/scalasbt/debian all main" \
      > /etc/apt/sources.list.d/sbt.list && \
    curl -sL "https://keyserver.ubuntu.com/pks/lookup?op=get&search=0x2EE0EA64E40A89B84B2DF73499E82A75642AC823" \
      | apt-key add - && \
    apt-get update && apt-get install -y --no-install-recommends sbt && \
    rm -rf /var/lib/apt/lists/*

WORKDIR /build
# Cache ivy/sbt dependencies before copying source
COPY project/build.properties project/plugins.sbt ./project/
COPY build.sbt ./
RUN sbt update

# Full source + compile + stage
COPY src ./src
RUN sbt stage

# ── Stage 2: minimal JRE runtime ─────────────────────────────────────────────
FROM eclipse-temurin:17-jre-jammy AS runtime

WORKDIR /app
COPY --from=builder /build/target/universal/stage ./

EXPOSE 8002

HEALTHCHECK --interval=15s --timeout=5s --start-period=20s --retries=3 \
    CMD curl -sf http://localhost:8002/health || exit 1

CMD ["./bin/kh-sim-scala"]
