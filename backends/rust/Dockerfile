# KH-SIM Rust backend — KH-003
# Port: 8001  |  Binary: kh-sim-rust  |  Framework: axum 0.7 / tokio
#
# Build:  docker build -t kh-sim-rust .
# Run:    docker run -p 8001:8001 kh-sim-rust
# Health: curl http://localhost:8001/health

# ── Stage 1: compile ─────────────────────────────────────────────────────────
# rust 1.86: Cargo.lock pins dev-dependencies (icu_* 2.2, idna_adapter via reqwest) that
# need rustc 1.86, so no test built in a 1.85 builder (2026-10-09, ref-v1 run 1).
FROM rust:1.86-bookworm AS builder

WORKDIR /build
COPY Cargo.toml Cargo.lock ./
# Cache dependency compilation layer separately
RUN mkdir src && echo "fn main(){}" > src/main.rs && \
    echo "pub fn lib(){}" > src/lib.rs && \
    cargo build --release && \
    rm -rf src

COPY src ./src
# Touch src to invalidate the cached empty-source layer
RUN touch src/main.rs src/lib.rs && cargo build --release

# ── Stage 2: minimal runtime ─────────────────────────────────────────────────
FROM debian:bookworm-slim AS runtime

RUN apt-get update && apt-get install -y --no-install-recommends \
    ca-certificates \
    curl \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /app
COPY --from=builder /build/target/release/kh-sim-rust ./kh-sim-rust

ENV RUST_LOG=info
EXPOSE 8001

HEALTHCHECK --interval=15s --timeout=3s --start-period=5s --retries=3 \
    CMD curl -sf http://localhost:8001/health || exit 1

CMD ["./kh-sim-rust"]
