/// golden_vectors.rs — per-backend golden-vector regression tests, Rust backend.
///
/// Mode 1: HTTP POST /simulate to running server on port 8001.
/// Compares backend response against oracle diagnostics + field_summaries using
/// Higham §4.2 hybrid tolerance with backend_factor = 1.0 (rust has no override entry).
///
/// Run (requires server running on :8001):
///   cd backends/rust
///   cargo run --release &
///   sleep 2
///   cargo test --test golden_vectors
///   kill %1
///
/// Each of the 9 test functions covers all three check types for that vector:
///   1. scalar diagnostics (KE, enstrophy, max_vort, div_rms; cfl_peak skipped if absent)
///   2. field summary statistics (min/max/mean per field)
///   3. field element comparison (first_10 / last_10 per field)
///
/// If the server is unreachable the test is skipped (not failed) with a message
/// printed to stdout — use RUST_TEST_NOCAPTURE=1 to see it.

use reqwest::blocking::Client;
use serde_json::Value;
use std::path::PathBuf;

const BACKEND_URL:  &str = "http://localhost:8001/simulate";
const BACKEND_NAME: &str = "kh-rust";
const MACHINE_EPS:  f64  = 2.22e-16;

// ── Vector loading ─────────────────────────────────────────────────────────────

fn golden_vectors_dir() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("..")
        .join("..")
        .join("shared")
        .join("golden-vectors")
}

fn load_vector(filename: &str) -> Value {
    let path = golden_vectors_dir().join(filename);
    let content = std::fs::read_to_string(&path)
        .unwrap_or_else(|e| panic!("cannot read {}: {e}", path.display()));
    serde_json::from_str(&content)
        .unwrap_or_else(|e| panic!("cannot parse {}: {e}", path.display()))
}

// ── Tolerance helpers (Higham §4.2 hybrid) ────────────────────────────────────

fn backend_factor(tols: &Value, kind: &str) -> f64 {
    tols["backend_overrides"][BACKEND_NAME]
        .get(format!("{kind}_factor"))
        .and_then(Value::as_f64)
        .unwrap_or(1.0)
}

fn scalar_tol(oracle: f64, abs_: f64, rel_: f64, factor: f64) -> f64 {
    (abs_ * factor).max(rel_ * factor * oracle.abs())
}

fn field_tol(oracle: f64, abs_: f64, rel_: f64, factor: f64, sqrt_n: f64, use_sqn: bool) -> f64 {
    let eps_term = if use_sqn { sqrt_n * MACHINE_EPS * oracle.abs() } else { 0.0 };
    (abs_ * factor).max(rel_ * factor * oracle.abs() + eps_term)
}

// ── Request builder ───────────────────────────────────────────────────────────

fn params_to_request(params: &Value) -> Value {
    serde_json::json!({
        "grid_nx":                params["grid_nx"],
        "grid_ny":                params["grid_ny"],
        "domain_lx":              params["domain_lx"],
        "domain_ly":              params["domain_ly"],
        "dt":                     params["dt"],
        "steps":                  params["nsteps"],   // REST field name differs from JSON key
        "reynolds_number":        params["reynolds_number"],
        "velocity_shear":         params["velocity_shear"],
        "perturbation_amplitude": params["perturbation_amplitude"],
        "perturbation_mode":      params["perturbation_mode"],
    })
}

// ── Core comparison logic ──────────────────────────────────────────────────────

struct CheckCtx<'a> {
    vector:  &'a Value,
    body:    &'a Value,
    name:    &'a str,
}

impl<'a> CheckCtx<'a> {
    fn check_scalars(&self) -> Vec<String> {
        let tols   = &self.vector["tolerances"];
        let oracle = &self.vector["diagnostics"];
        let sf     = backend_factor(tols, "scalar");
        let s_abs  = tols["scalar_abs"].as_f64().unwrap();
        let s_rel  = tols["scalar_rel"].as_f64().unwrap();
        let sd_abs = tols["scalar_div_rms_abs"].as_f64().unwrap();
        let sd_rel = tols["scalar_div_rms_rel"].as_f64().unwrap();

        let b_diag = &self.body["diagnostics"];
        let mut failures = Vec::new();

        let checks: &[(&str, f64, f64)] = &[
            ("kinetic_energy", s_abs,  s_rel),
            ("enstrophy",      s_abs,  s_rel),
            ("max_vorticity",  s_abs,  s_rel),
            ("divergence_rms", sd_abs, sd_rel),
        ];
        for &(key, abs_, rel_) in checks {
            let bv = match b_diag[key].as_f64() {
                Some(v) => v,
                None    => { failures.push(format!("{key}: missing from response")); continue; }
            };
            let ov  = oracle[key].as_f64().unwrap();
            let tol = scalar_tol(ov, abs_, rel_, sf);
            let err = (bv - ov).abs();
            if err >= tol {
                failures.push(format!(
                    "{key}: |diff|={err:.3e} >= tol={tol:.3e} (actual={bv:.6e}, oracle={ov:.6e})"
                ));
            }
        }
        // cfl_peak: skip gracefully if backend does not return it
        if let Some(bv) = b_diag["cfl_peak"].as_f64() {
            let ov  = oracle["cfl_peak"].as_f64().unwrap();
            let tol = scalar_tol(ov, s_abs, s_rel, sf);
            let err = (bv - ov).abs();
            if err >= tol {
                failures.push(format!(
                    "cfl_peak: |diff|={err:.3e} >= tol={tol:.3e} (actual={bv:.6e}, oracle={ov:.6e})"
                ));
            }
        }
        failures
    }

    fn check_field_summaries_and_elements(&self) -> Vec<String> {
        let tols     = &self.vector["tolerances"];
        let oracle   = &self.vector["field_summaries"];
        let ff       = backend_factor(tols, "field");
        let fs_abs   = tols["field_summary_abs"].as_f64().unwrap();
        let fs_rel   = tols["field_summary_rel"].as_f64().unwrap();
        let fe_abs   = tols["field_element_abs"].as_f64().unwrap();
        let fe_rel   = tols["field_element_rel"].as_f64().unwrap();
        let nx       = self.vector["parameters"]["grid_nx"].as_u64().unwrap() as f64;
        let ny       = self.vector["parameters"]["grid_ny"].as_u64().unwrap() as f64;
        let sqrt_n   = (nx * ny).sqrt();
        let use_sqn  = tols["sqrt_n_scaling"].as_bool().unwrap_or(false);

        let fields = [
            ("u_velocity", "u_velocity"),
            ("v_velocity", "v_velocity"),
            ("vorticity",  "vorticity"),
            ("pressure",   "pressure"),
        ];
        let mut failures = Vec::new();

        for (oracle_key, body_key) in &fields {
            let o_fs = &oracle[oracle_key];
            let b_field = &self.body[body_key];

            // Backend returns full Vec<f64> — compute summary on-the-fly
            let flat: Vec<f64> = match serde_json::from_value(b_field.clone()) {
                Ok(v) => v,
                Err(_) => { failures.push(format!("{body_key}: not a float array")); continue; }
            };
            if flat.is_empty() {
                failures.push(format!("{body_key}: empty array"));
                continue;
            }

            // Summary: min/max/mean
            let b_min  = flat.iter().cloned().fold(f64::INFINITY, f64::min);
            let b_max  = flat.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
            let b_mean = flat.iter().sum::<f64>() / flat.len() as f64;

            for (stat, bv, ov_key) in [
                ("min",  b_min,  "min"),
                ("max",  b_max,  "max"),
                ("mean", b_mean, "mean"),
            ] {
                let ov  = o_fs[ov_key].as_f64().unwrap();
                let tol = field_tol(ov, fs_abs, fs_rel, ff, sqrt_n, use_sqn);
                let err = (bv - ov).abs();
                if err >= tol {
                    failures.push(format!(
                        "{oracle_key}.{stat}: |diff|={err:.3e} >= tol={tol:.3e} \
                         (actual={bv:.6e}, oracle={ov:.6e})"
                    ));
                }
            }

            // Elements: first_10 / last_10
            let o_first: Vec<f64> = serde_json::from_value(o_fs["first_10"].clone()).unwrap();
            let o_last:  Vec<f64> = serde_json::from_value(o_fs["last_10"].clone()).unwrap();
            let b_first = &flat[..o_first.len().min(flat.len())];
            let b_last  = &flat[flat.len().saturating_sub(o_last.len())..];

            for (slot, b_elems, o_elems) in [
                ("first_10", b_first, o_first.as_slice()),
                ("last_10",  b_last,  o_last.as_slice()),
            ] {
                for (i, (&bv, &ov)) in b_elems.iter().zip(o_elems.iter()).enumerate() {
                    let tol = field_tol(ov, fe_abs, fe_rel, ff, sqrt_n, use_sqn);
                    let err = (bv - ov).abs();
                    if err >= tol {
                        failures.push(format!(
                            "{oracle_key}.{slot}[{i}]: |diff|={err:.3e} >= tol={tol:.3e} \
                             (actual={bv:.6e}, oracle={ov:.6e})"
                        ));
                    }
                }
            }
        }
        failures
    }
}

// ── HTTP test helper ───────────────────────────────────────────────────────────

/// Post the vector's parameters to the backend, run all checks.
/// Returns Ok(()) on pass, Err(msg) on failure, Err("SKIP: ...") if server unreachable.
fn run_vector(filename: &str) -> Result<(), String> {
    let vector  = load_vector(filename);
    let payload = params_to_request(&vector["parameters"]);

    let client = Client::builder()
        .timeout(std::time::Duration::from_secs(120))
        .build()
        .expect("http client build");

    let resp = client
        .post(BACKEND_URL)
        .json(&payload)
        .send();

    let resp = match resp {
        Ok(r) => r,
        Err(e) if e.is_connect() => {
            return Err(format!("SKIP: server not reachable at {BACKEND_URL}: {e}"));
        }
        Err(e) => return Err(format!("HTTP error: {e}")),
    };

    if !resp.status().is_success() {
        return Err(format!("HTTP {}: {}", resp.status(), resp.text().unwrap_or_default()));
    }

    let body: Value = resp.json().map_err(|e| format!("response parse: {e}"))?;
    let ctx = CheckCtx { vector: &vector, body: &body, name: filename };

    let mut all_failures = ctx.check_scalars();
    all_failures.extend(ctx.check_field_summaries_and_elements());

    if all_failures.is_empty() {
        Ok(())
    } else {
        Err(format!(
            "{filename}: {} failure(s):\n{}",
            all_failures.len(),
            all_failures.iter().map(|f| format!("  {f}")).collect::<Vec<_>>().join("\n")
        ))
    }
}

/// When KH_GOLDEN_REQUIRE_SERVER=1 a missing server is a hard panic.
/// When unset (default) the test prints the skip message and passes —
/// safe for CI environments without the Rust backend running.
fn require_server() -> bool {
    std::env::var("KH_GOLDEN_REQUIRE_SERVER").map(|v| v == "1").unwrap_or(false)
}

macro_rules! golden_test {
    ($name:ident, $file:expr) => {
        #[test]
        fn $name() {
            match run_vector($file) {
                Ok(()) => {}
                Err(msg) if msg.starts_with("SKIP:") => {
                    if require_server() {
                        panic!("{msg}");
                    }
                    println!("{msg}");
                }
                Err(msg) => panic!("{msg}"),
            }
        }
    };
}

// ── 9 parametrised test cases (one per vector) ────────────────────────────────

golden_test!(test_v1_canonical_reference,  "v1-canonical-reference.json");
golden_test!(test_v2_low_re_diffusive,     "v2-low-re-diffusive.json");
golden_test!(test_v3_high_re_nonlinear,    "v3-high-re-nonlinear.json");
golden_test!(test_v4_double_resolution,    "v4-double-resolution.json");
golden_test!(test_v5_half_resolution_fast, "v5-half-resolution-fast.json");
golden_test!(test_v6_long_time_evolution,  "v6-long-time-evolution.json");
golden_test!(test_v7_alt_mode_k1,          "v7-alt-mode-k1.json");
golden_test!(test_v8_dealiasing_edge_k4,  "v8-dealiasing-edge-k4.json");
golden_test!(test_v9_zero_perturbation,    "v9-zero-perturbation.json");
