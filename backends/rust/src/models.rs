/// models.rs — REST API data structures (mirrors openapi.yaml / kh_physics.py dataclasses)
use serde::{Deserialize, Serialize};

// ── Request ──────────────────────────────────────────────────────────────────

#[derive(Debug, Deserialize)]
pub struct SimulationRequest {
    #[serde(default = "default_grid_nx")]
    pub grid_nx: usize,

    #[serde(default = "default_grid_ny")]
    pub grid_ny: usize,

    #[serde(default = "default_domain_lx")]
    pub domain_lx: f64,

    #[serde(default = "default_domain_ly")]
    pub domain_ly: f64,

    #[serde(default = "default_dt")]
    pub dt: f64,

    #[serde(default = "default_steps")]
    pub steps: usize,

    #[serde(default = "default_reynolds_number")]
    pub reynolds_number: f64,

    #[serde(default = "default_velocity_shear")]
    pub velocity_shear: f64,

    #[serde(default = "default_perturbation_amplitude")]
    pub perturbation_amplitude: f64,

    #[serde(default = "default_perturbation_mode")]
    pub perturbation_mode: usize,

    /// Optional: pre-initialised vorticity (row-major flattened, length nx*ny).
    /// If None, the analytic IC is used.
    pub initial_omega: Option<Vec<f64>>,

    // -- Phase 7f forcing extension --------------------------------------------
    // PHASE-7F-FORCING-API-DESIGN-v0.1.md (sec.1 v0.1 + sec.7 v0.2). Identical
    // field names + formulas across all six backends. Defaults reproduce pre-7f
    // behaviour EXACTLY (no forcing, no drag) -> fully backward compatible:
    // phase-6c golden vectors and phase-7d testbeds are unaffected.
    #[serde(default = "default_forcing_type")]
    pub forcing_type: String, // "none" | "mode" | "band"

    #[serde(default)]
    pub forcing_amplitude: f64, // A

    #[serde(default = "default_forcing_k")]
    pub forcing_kx: i64, // "mode" only: integer mode along x

    #[serde(default = "default_forcing_k")]
    pub forcing_ky: i64, // "mode" only: integer mode along y

    #[serde(default)]
    pub forcing_kmin: f64, // "band" only: inner integer-mode radius

    #[serde(default)]
    pub forcing_kmax: f64, // "band" only: outer integer-mode radius

    #[serde(default)]
    pub forcing_omega0: f64, // v0.2: phase-rotation rate scale (band only)

    #[serde(default)]
    pub linear_drag: f64, // v0.2: Ekman-type linear drag coefficient alpha

    /// Time of the initial field (2026-10-10, CP-7F-B/D harness): time-dependent
    /// forcing runs from t_start, and t_final = t_start + steps*dt.
    #[serde(default)]
    pub t_start: f64,
}

fn default_forcing_type() -> String { "none".to_string() }
fn default_forcing_k()    -> i64    { 1 }

fn default_grid_nx()                -> usize { 128 }
fn default_grid_ny()                -> usize { 64  }
fn default_domain_lx()              -> f64   { 1.0 }
fn default_domain_ly()              -> f64   { 0.5 }
fn default_dt()                     -> f64   { 0.001 }
fn default_steps()                  -> usize { 100 }
fn default_reynolds_number()        -> f64   { 1000.0 }
fn default_velocity_shear()         -> f64   { 1.0 }
fn default_perturbation_amplitude() -> f64   { 0.01 }
fn default_perturbation_mode()      -> usize { 2 }

// ── Response ─────────────────────────────────────────────────────────────────

#[derive(Debug, Serialize)]
pub struct Diagnostics {
    pub kinetic_energy:  f64,
    pub enstrophy:       f64,
    pub max_vorticity:   f64,
    pub divergence_rms:  f64,
    pub cfl_peak:        f64,   // max over all steps of dt*(max|u|/dx + max|v|/dy)
}

#[derive(Debug, Serialize)]
pub struct SimulationResult {
    pub backend:         &'static str,
    pub language:        &'static str,
    pub steps_completed: usize,
    pub t_final:         f64,
    pub grid_nx:         usize,
    pub grid_ny:         usize,
    pub u_velocity:      Vec<f64>,
    pub v_velocity:      Vec<f64>,
    pub vorticity:       Vec<f64>,
    pub pressure:        Vec<f64>,   // streamfunction psi
    pub diagnostics:     Diagnostics,
    pub compute_time_ms: f64,
    pub warnings:        Vec<String>,   // KH-VAL-2 K5 (filled by the HTTP layer)
}

// ── Health / Info ─────────────────────────────────────────────────────────────

#[derive(Debug, Serialize)]
pub struct HealthResponse {
    pub status:  &'static str,
    pub backend: &'static str,
    pub port:    u16,
}

#[derive(Debug, Serialize)]
pub struct InfoResponse {
    pub backend:         &'static str,
    pub language:        &'static str,
    pub framework:       &'static str,
    pub fft_library:     &'static str,
    pub grid_constraint: &'static str,  // "any": rustfft transforms every length
    pub build_date:      String,        // YYYY-MM-DD, the binary's own timestamp
    pub port:            u16,
    pub openapi_spec:    &'static str,
}
