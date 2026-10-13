# KH-SIM Fortran backend — KH-006
# Port: 8004  |  Binary: kh-sim-fortran
# Language: Fortran 2008 (physics) + C++20 HTTP shim (cpp-httplib)
# Compiler: gfortran 12 + g++ 12  |  Build: CMake + Ninja + FetchContent
#
# Build:  docker build -t kh-sim-fortran .
# Run:    docker run -p 8004:8004 kh-sim-fortran
# Health: curl http://localhost:8004/health

# ── Stage 1: build ───────────────────────────────────────────────────────────
FROM debian:bookworm AS builder

RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    gfortran \
    cmake \
    ninja-build \
    git \
    ca-certificates \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /build
COPY CMakeLists.txt ./
COPY src ./src

# Configure + build (FetchContent pulls cpp-httplib + nlohmann/json)
RUN cmake -S . -B build \
      -G Ninja \
      -DCMAKE_BUILD_TYPE=Release \
      -DCMAKE_Fortran_COMPILER=gfortran \
      -DCMAKE_CXX_COMPILER=g++ && \
    cmake --build build --target kh-sim-fortran --parallel

# ── Stage 2: minimal runtime ─────────────────────────────────────────────────
FROM debian:bookworm-slim AS runtime

RUN apt-get update && apt-get install -y --no-install-recommends \
    libgfortran5 \
    libstdc++6 \
    libgomp1 \
    ca-certificates \
    curl \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /app
COPY --from=builder /build/build/kh-sim-fortran ./kh-sim-fortran

EXPOSE 8004

HEALTHCHECK --interval=15s --timeout=3s --start-period=5s --retries=3 \
    CMD curl -sf http://localhost:8004/health || exit 1

CMD ["./kh-sim-fortran"]
