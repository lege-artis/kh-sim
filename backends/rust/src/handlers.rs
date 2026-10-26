/// handlers.rs -- axum route handlers
///
/// Endpoints (per openapi.yaml):
///   POST /simulate   -> run KH simulation, return SimulationResult JSON
///   GET  /health     -> liveness probe
///   GET  /info       -> backend metadata
///
/// 2026-10-03 (finding v0.5, R-1..R-3): the body is read as bytes and validated by
/// contract.rs (KH-VAL-1) instead of axum's Json extractor. Until then a refused
/// request came back as 422 with a plain-text body, grid_nx = 0 panicked the
/// physics task and dropped the connection, and a blown-up run returned 200 with
/// null numbers.
use axum::{body::Bytes, extract::rejection::BytesRejection, http::StatusCode, response::IntoResponse, Json};
use serde_json::json;
use tracing::{info, warn};

use crate::contract;
use crate::models::{HealthResponse, InfoResponse, SimulationResult};
use crate::physics;

const BACKEND_ID: &str = "kh-rust";

fn error(status: StatusCode, message: &str, key: &str, value: &str) -> axum::response::Response {
    (status, Json(json!({ "error": message, key: value }))).into_response()
}

fn is_finite(r: &SimulationResult) -> bool {
    let d = &r.diagnostics;
    [d.kinetic_energy, d.enstrophy, d.max_vorticity, d.divergence_rms, d.cfl_peak]
        .iter()
        .chain(r.u_velocity.iter())
        .chain(r.v_velocity.iter())
        .chain(r.vorticity.iter())
        .chain(r.pressure.iter())
        .all(|x| x.is_finite())
}

// -- POST /simulate -------------------------------------------------------------

pub async fn simulate(body: Result<Bytes, BytesRejection>) -> axum::response::Response {
    // A body above contract::MAX_BODY_BYTES (no request the contract allows) or one that
    // could not be read: 400 in the contract's form, not axum's plain-text 413.
    let body = match body {
        Ok(b) => b,
        Err(e) => return error(StatusCode::BAD_REQUEST,
                               &format!("request body not read: {e} (limit {} bytes: a request at \
                                         the {}-cell limit with initial_omega)",
                                        contract::MAX_BODY_BYTES, contract::MAX_CELLS),
                               "field", "body"),
    };
    let req = match contract::parse(&body) {
        Ok(r) => r,
        Err(e) => return error(StatusCode::BAD_REQUEST, &e.message, "field", &e.field),
    };

    info!(
        nx = req.grid_nx,
        ny = req.grid_ny,
        steps = req.steps,
        re = req.reynolds_number,
        "simulate request received"
    );

    let wargs = (req.velocity_shear, req.reynolds_number, req.dt,
                 req.domain_lx / req.grid_nx as f64, req.domain_ly / req.grid_ny as f64);

    // Run physics on a blocking thread so the async executor is not starved.
    // A panic there is the backend's failure: 500, and the server keeps serving.
    let mut result = match tokio::task::spawn_blocking(move || physics::simulate(&req)).await {
        Ok(r) => r,
        Err(e) => {
            warn!(error = %e, "physics task failed");
            return error(StatusCode::INTERNAL_SERVER_ERROR,
                         &format!("computation failed: {e}"), "backend", BACKEND_ID);
        }
    };
    if !is_finite(&result) {
        return error(StatusCode::INTERNAL_SERVER_ERROR,
                     "computation produced a non-finite value (dt beyond the stability bound?)",
                     "backend", BACKEND_ID);
    }

    info!(
        t_final         = result.t_final,
        ke              = result.diagnostics.kinetic_energy,
        enstrophy       = result.diagnostics.enstrophy,
        max_vorticity   = result.diagnostics.max_vorticity,
        divergence_rms  = result.diagnostics.divergence_rms,
        compute_ms      = result.compute_time_ms,
        "simulate complete"
    );

    result.warnings = contract::warnings(wargs.0, wargs.1, wargs.2, wargs.3, wargs.4,
                                         result.diagnostics.cfl_peak);
    (StatusCode::OK, Json(result)).into_response()
}

// -- GET /health ----------------------------------------------------------------

pub async fn health() -> impl IntoResponse {
    Json(HealthResponse {
        status:  "ok",
        backend: "rust-axum",
        port:    8001,
    })
}

// -- GET /info ------------------------------------------------------------------

/// The executable's modification date (UTC, YYYY-MM-DD). In the image that is the
/// build date, which lets a deploy tell a fresh build from a stale one without a
/// build script; "unknown" if the file cannot be read.
fn build_date() -> String {
    let secs = std::env::current_exe()
        .and_then(|p| p.metadata())
        .and_then(|m| m.modified())
        .ok()
        .and_then(|t| t.duration_since(std::time::UNIX_EPOCH).ok())
        .map(|d| d.as_secs() as i64);
    let Some(secs) = secs else { return "unknown".to_string() };
    // Civil date from days since 1970-01-01 (H. Hinnant's days_from_civil inverse).
    let z = secs.div_euclid(86_400) + 719_468;
    let era = z.div_euclid(146_097);
    let doe = z - era * 146_097;
    let yoe = (doe - doe / 1_460 + doe / 36_524 - doe / 146_096) / 365;
    let doy = doe - (365 * yoe + yoe / 4 - yoe / 100);
    let mp = (5 * doy + 2) / 153;
    let day = doy - (153 * mp + 2) / 5 + 1;
    let month = if mp < 10 { mp + 3 } else { mp - 9 };
    let year = yoe + era * 400 + i64::from(month <= 2);
    format!("{year:04}-{month:02}-{day:02}")
}

pub async fn info() -> impl IntoResponse {
    Json(InfoResponse {
        backend:         "rust-axum",
        language:        "Rust",
        framework:       "axum 0.7",
        fft_library:     "rustfft 6",
        grid_constraint: "any",
        build_date:      build_date(),
        port:            8001,
        openapi_spec:    "kh-sim/shared/api/openapi.yaml",
    })
}
