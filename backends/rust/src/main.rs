/// main.rs — kh-sim Rust/axum backend entry point
///
/// Starts HTTP server on 0.0.0.0:8001 (port per kh-sim.config.yaml).
/// Routes:
///   POST /simulate
///   GET  /health
///   GET  /info
///
/// Task: KH-003
/// Spec: kh-sim/shared/api/openapi.yaml
use kh_sim_rust::{contract, handlers};

use axum::{extract::DefaultBodyLimit, routing::{get, post}, Router};
use tower_http::{
    cors::{Any, CorsLayer},
    trace::TraceLayer,
};
use tracing_subscriber::{layer::SubscriberExt, util::SubscriberInitExt, EnvFilter};

#[tokio::main]
async fn main() {
    // Structured logging — honour RUST_LOG env var, default=info
    tracing_subscriber::registry()
        .with(EnvFilter::try_from_default_env().unwrap_or_else(|_| "info".into()))
        .with(tracing_subscriber::fmt::layer())
        .init();

    let cors = CorsLayer::new()
        .allow_origin(Any)
        .allow_headers(Any)
        .allow_methods(Any);

    let app = Router::new()
        .route("/simulate", post(handlers::simulate))
        .route("/health",   get(handlers::health))
        .route("/info",     get(handlers::info))
        // Every request the contract allows fits (contract::MAX_BODY_BYTES); a larger
        // body is refused by the handler as a 400 naming "body".
        .layer(DefaultBodyLimit::max(contract::MAX_BODY_BYTES))
        .layer(TraceLayer::new_for_http())
        .layer(cors);

    let addr = "0.0.0.0:8001";
    tracing::info!("kh-sim Rust backend listening on http://{}", addr);

    let listener = tokio::net::TcpListener::bind(addr).await.unwrap();
    axum::serve(listener, app).await.unwrap();
}
