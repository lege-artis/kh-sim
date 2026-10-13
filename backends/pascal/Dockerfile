# KH-SIM Pascal backend — KH-007
# Port: 8005  |  Binary: kh_server  |  Compiler: FPC 3.2.2 (Debian bookworm)
# Units: fphttpserver, fpjson, jsonparser, Math, SysUtils
#
# Build:  docker build -t kh-sim-pascal .
# Run:    docker run -p 8005:8005 kh-sim-pascal
# Health: curl http://localhost:8005/health
#
# NOTE: FPC unit availability on Linux is confirmed equivalent to Windows FPC 3.2.2
# for the units used: fphttpserver, fpjson, jsonparser, Math, SysUtils, DateUtils.
# The fphttpserver unit in Debian FPC uses the same TFPHttpServer / TFPHTTPConnectionRequest
# interface as the Windows build.

# ── Stage 1: compile ─────────────────────────────────────────────────────────
FROM debian:bookworm AS builder

RUN apt-get update && apt-get install -y --no-install-recommends \
    fpc \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /build
COPY src ./src

# Build physics unit first (kh_physics.pas), then HTTP server
# -O2: optimise  -Mobjfpc: ObjFPC dialect  -Fu./src: unit search path
# -FE./build: output dir
# NOTE: -Cs (stack size override) is a Windows-only FPC flag; ppcx64 on Linux
#       rejects it with "Illegal parameter". Linux stack is managed by ulimit.
#       Default 8 MB thread stack is sufficient for 64x32 FFT recursion depth.
RUN mkdir -p build && \
    fpc -O2 -Mobjfpc \
        -Fu./src \
        -FE./build \
        src/kh_server.pas

# ── Stage 2: minimal runtime ─────────────────────────────────────────────────
FROM debian:bookworm-slim AS runtime

# Pascal binaries are statically linked by default -- no runtime libs needed
# except libc which bookworm-slim provides.
# curl is required by the HEALTHCHECK CMD below.
RUN apt-get update && apt-get install -y --no-install-recommends \
    ca-certificates \
    curl \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /app
COPY --from=builder /build/build/kh_server ./kh_server

EXPOSE 8005

HEALTHCHECK --interval=15s --timeout=3s --start-period=5s --retries=3 \
    CMD curl -sf http://localhost:8005/health || exit 1

CMD ["./kh_server"]
