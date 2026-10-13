# KH-SIM C++ backend — KH-005
# Port: 8003  |  Binary: kh-sim-cpp  |  Standard: C++20  |  HTTP: cpp-httplib 0.18
# Dependencies: fetched at build time via CMake FetchContent (internet required during build)
#
# Build:  docker build -t kh-sim-cpp .
# Run:    docker run -p 8003:8003 kh-sim-cpp
# Health: curl http://localhost:8003/health

# ── Stage 1: build ───────────────────────────────────────────────────────────
FROM debian:bookworm AS builder

RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    cmake \
    ninja-build \
    git \
    ca-certificates \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /build
COPY CMakeLists.txt ./
COPY src ./src

# Configure + build Release (FetchContent pulls httplib + nlohmann/json at this step)
RUN cmake -S . -B build \
      -G Ninja \
      -DCMAKE_BUILD_TYPE=Release \
      -DCMAKE_CXX_STANDARD=20 && \
    cmake --build build --target kh-sim-cpp --parallel

# ── Stage 2: minimal runtime ─────────────────────────────────────────────────
FROM debian:bookworm-slim AS runtime

RUN apt-get update && apt-get install -y --no-install-recommends \
    libstdc++6 \
    ca-certificates \
    curl \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /app
COPY --from=builder /build/build/kh-sim-cpp ./kh-sim-cpp

EXPOSE 8003

HEALTHCHECK --interval=15s --timeout=3s --start-period=5s --retries=3 \
    CMD curl -sf http://localhost:8003/health || exit 1

CMD ["./kh-sim-cpp"]
