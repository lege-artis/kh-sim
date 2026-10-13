/// contract.rs -- POST /simulate request contract KH-VAL-1 (Rust implementation)
///
/// Same rules as backends/cpp/src/kh_http_contract.hpp and the Pascal server
/// (finding v0.5, kh-sim/_audit/FINDING-PASCAL-HTTP-PARAMETER-DROP-2026-09-16.md):
///   - body must be a JSON object                          -> 400 field "body"
///   - a key present with the wrong JSON type               -> 400 field <key>
///   - grid_nx, grid_ny in [1, 65536] (rustfft takes any length, so no
///     power-of-two rule); grid_nx*grid_ny <= 4194304      -> 400 field grid_nx|grid_ny
///   - domain_lx, domain_ly, reynolds_number non-zero      -> 400 field <key>
///   - initial_omega (legacy alias init_omega, not both)
///     has exactly grid_nx*grid_ny entries                 -> 400 field initial_omega
///   - forcing_type in {none, mode, band}; band needs
///     domain_lx == domain_ly                              -> 400 field forcing_type
/// Rust-only, because the solver's types are unsigned: a negative
/// perturbation_mode is refused (400) before the grid checks.
///
/// Contract KH-VAL-2 (kh-sim/_audit/KH-SIM-MODELLING-RULES-v0.1.md, ruling
/// 2026-10-07), after the KH-VAL-1 checks, same rules as
/// backends/cpp/src/kh_http_contract.hpp:
///   K3 outside the model: domain_lx, domain_ly, reynolds_number < 0; dt < 0;
///      steps < 0; linear_drag < 0; band forcing with kmin > kmax -> 400 field
///   K4 not representable (analytic field, amplitude != 0): perturbation_mode 0
///      or beyond the 2/3-rule band (3*mode > 2*(grid_nx/2))      -> 400 field
///   K4b (ruling 2026-10-09) a non-zero perturbation_amplitude sent with
///      initial_omega (it would be silently ignored)               -> 400 field
///   K4 forcing (amplitude != 0, 2026-10-10): mode forcing_kx/ky beyond the
///      2/3-rule band, or a band with 3*floor(kmax) > 2*(min(nx, ny)/2) -> 400 field
///   K3 a band with no integer mode kmin <= |k| <= kmax             -> 400 forcing_kmax
///   K5 warnings in the response (warnings()); K7 unknown keys logged.
/// t_start (2026-10-10): the time of the initial field (forcing phase, t_final).
use serde_json::{Map, Value};

use crate::models::SimulationRequest;

pub const MAX_AXIS: i64 = 65536;
pub const MAX_CELLS: i64 = 4_194_304;

/// A request the backend refuses: HTTP 400 {"error": message, "field": field}.
#[derive(Debug)]
pub struct Invalid {
    pub field: String,
    pub message: String,
}

fn invalid(field: &str, message: impl Into<String>) -> Invalid {
    Invalid { field: field.to_string(), message: message.into() }
}

fn type_name(v: &Value) -> &'static str {
    match v {
        Value::Null => "null",
        Value::Bool(_) => "boolean",
        Value::Number(_) => "number",
        Value::String(_) => "string",
        Value::Array(_) => "array",
        Value::Object(_) => "object",
    }
}

fn present<'a>(o: &'a Map<String, Value>, key: &str) -> Option<&'a Value> {
    match o.get(key) {
        None | Some(Value::Null) => None,
        Some(v) => Some(v),
    }
}

fn read_int(o: &Map<String, Value>, key: &str) -> Result<Option<i64>, Invalid> {
    let Some(v) = present(o, key) else { return Ok(None) };
    let n = if let Some(i) = v.as_i64() {
        i
    } else if let Some(f) = v.as_f64() {
        if !f.is_finite() || f != f.floor() || f.abs() > 2_147_483_647.0 {
            return Err(invalid(key, format!("{key} must be an integer")));
        }
        f as i64
    } else {
        return Err(invalid(key, format!("{key} must be an integer, got {}", type_name(v))));
    };
    if !(-2_147_483_647..=2_147_483_647).contains(&n) {
        return Err(invalid(key, format!("{key} is outside the 32-bit integer range")));
    }
    Ok(Some(n))
}

fn read_f64(o: &Map<String, Value>, key: &str) -> Result<Option<f64>, Invalid> {
    let Some(v) = present(o, key) else { return Ok(None) };
    v.as_f64()
        .map(Some)
        .ok_or_else(|| invalid(key, format!("{key} must be a number, got {}", type_name(v))))
}

fn read_string(o: &Map<String, Value>, key: &str) -> Result<Option<String>, Invalid> {
    let Some(v) = present(o, key) else { return Ok(None) };
    v.as_str()
        .map(|s| Some(s.to_string()))
        .ok_or_else(|| invalid(key, format!("{key} must be a string, got {}", type_name(v))))
}

fn read_field(o: &Map<String, Value>, key: &str) -> Result<Option<Vec<f64>>, Invalid> {
    let Some(v) = present(o, key) else { return Ok(None) };
    let arr = v
        .as_array()
        .ok_or_else(|| invalid("initial_omega", format!("{key} must be an array of numbers")))?;
    if arr.len() as i64 > MAX_CELLS {
        return Err(invalid("initial_omega",
            format!("{key} has more entries than the {MAX_CELLS}-cell limit")));
    }
    arr.iter()
        .map(|x| x.as_f64().ok_or_else(|| invalid("initial_omega",
            format!("{key} must contain numbers only"))))
        .collect::<Result<Vec<f64>, Invalid>>()
        .map(Some)
}

fn check_axis(key: &str, n: i64) -> Result<usize, Invalid> {
    if !(1..=MAX_AXIS).contains(&n) {
        return Err(invalid(key, format!("{key} = {n} is out of range [1, {MAX_AXIS}]")));
    }
    Ok(n as usize)
}

/// Parse and validate a POST /simulate body into the solver's request.
pub fn parse(body: &[u8]) -> Result<SimulationRequest, Invalid> {
    let v: Value = serde_json::from_slice(body)
        .map_err(|_| invalid("body", "request body is not valid JSON"))?;
    let o = v
        .as_object()
        .ok_or_else(|| invalid("body", "request body must be a JSON object"))?;

    // Start from the serde defaults, so a key that is absent means exactly what
    // it meant before this module existed.
    let mut r: SimulationRequest = serde_json::from_value(Value::Object(Map::new()))
        .expect("SimulationRequest has a default for every field");

    let nx = read_int(o, "grid_nx")?.unwrap_or(r.grid_nx as i64);
    let ny = read_int(o, "grid_ny")?.unwrap_or(r.grid_ny as i64);
    if let Some(x) = read_f64(o, "domain_lx")? { r.domain_lx = x; }
    if let Some(x) = read_f64(o, "domain_ly")? { r.domain_ly = x; }
    if let Some(x) = read_f64(o, "dt")? { r.dt = x; }
    let steps_raw = read_int(o, "steps")?;
    if let Some(n) = steps_raw { r.steps = n.max(0) as usize; }
    if let Some(x) = read_f64(o, "reynolds_number")? { r.reynolds_number = x; }
    if let Some(x) = read_f64(o, "velocity_shear")? { r.velocity_shear = x; }
    if let Some(x) = read_f64(o, "perturbation_amplitude")? { r.perturbation_amplitude = x; }
    if let Some(n) = read_int(o, "perturbation_mode")? {
        if n < 0 {
            return Err(invalid("perturbation_mode",
                "perturbation_mode must be >= 0 on this backend"));
        }
        r.perturbation_mode = n as usize;
    }
    if let Some(s) = read_string(o, "forcing_type")? { r.forcing_type = s; }
    if let Some(x) = read_f64(o, "forcing_amplitude")? { r.forcing_amplitude = x; }
    if let Some(n) = read_int(o, "forcing_kx")? { r.forcing_kx = n; }
    if let Some(n) = read_int(o, "forcing_ky")? { r.forcing_ky = n; }
    if let Some(x) = read_f64(o, "forcing_kmin")? { r.forcing_kmin = x; }
    if let Some(x) = read_f64(o, "forcing_kmax")? { r.forcing_kmax = x; }
    if let Some(x) = read_f64(o, "forcing_omega0")? { r.forcing_omega0 = x; }
    if let Some(x) = read_f64(o, "linear_drag")? { r.linear_drag = x; }
    if let Some(x) = read_f64(o, "t_start")? { r.t_start = x; }

    let canonical = read_field(o, "initial_omega")?;
    let legacy = read_field(o, "init_omega")?;
    if canonical.is_some() && legacy.is_some() {
        return Err(invalid("initial_omega",
            "send initial_omega or the legacy init_omega, not both"));
    }
    r.initial_omega = canonical.or(legacy);

    r.grid_nx = check_axis("grid_nx", nx)?;
    r.grid_ny = check_axis("grid_ny", ny)?;
    let cells = nx * ny;
    if cells > MAX_CELLS {
        return Err(invalid("grid_nx", format!(
            "grid {nx}x{ny} has {cells} cells, more than the {MAX_CELLS}-cell limit")));
    }
    if r.domain_lx == 0.0 { return Err(invalid("domain_lx", "domain_lx must be non-zero")); }
    if r.domain_ly == 0.0 { return Err(invalid("domain_ly", "domain_ly must be non-zero")); }
    if r.reynolds_number == 0.0 {
        return Err(invalid("reynolds_number", "reynolds_number must be non-zero"));
    }
    if let Some(w) = &r.initial_omega {
        if w.len() as i64 != cells {
            return Err(invalid("initial_omega", format!(
                "initial_omega has {} entries; grid {nx}x{ny} needs {cells}", w.len())));
        }
    }
    match r.forcing_type.as_str() {
        "none" | "mode" => {}
        "band" => {
            // Same test as the solver's own assert (physics/solver.rs, band arm),
            // so a request the solver would abort on is refused here instead.
            if r.forcing_amplitude != 0.0 && (r.domain_lx - r.domain_ly).abs() >= 1e-12 {
                return Err(invalid("forcing_type", "band forcing requires domain_lx == domain_ly"));
            }
        }
        other => {
            return Err(invalid("forcing_type", format!(
                "unknown forcing_type '{other}'; expected none, mode or band")));
        }
    }

    // KH-VAL-2, K3: outside the model.
    if r.domain_lx < 0.0 { return Err(invalid("domain_lx", "domain_lx must be positive")); }
    if r.domain_ly < 0.0 { return Err(invalid("domain_ly", "domain_ly must be positive")); }
    if r.reynolds_number < 0.0 {
        return Err(invalid("reynolds_number",
            "reynolds_number must be positive (a negative value is negative viscosity)"));
    }
    if r.dt < 0.0 {
        return Err(invalid("dt", "dt must be >= 0 (backward integration of a viscous flow is ill-posed)"));
    }
    if matches!(steps_raw, Some(n) if n < 0) { return Err(invalid("steps", "steps must be >= 0")); }
    if r.linear_drag < 0.0 {
        return Err(invalid("linear_drag", "linear_drag must be >= 0 (negative drag injects energy)"));
    }
    if r.forcing_type == "band" && r.forcing_amplitude != 0.0 && r.forcing_kmin > r.forcing_kmax {
        return Err(invalid("forcing_kmin", "forcing_kmin > forcing_kmax: the forcing band is empty"));
    }
    // KH-VAL-2, K4 (forcing): the forced modes must lie inside the 2/3-rule band.
    if r.forcing_type == "mode" && r.forcing_amplitude != 0.0 {
        if 3 * r.forcing_kx.unsigned_abs() as i64 > 2 * (nx / 2) {
            return Err(invalid("forcing_kx", format!(
                "forcing_kx {} lies beyond the 2/3-rule band of grid_nx {nx}", r.forcing_kx)));
        }
        if 3 * r.forcing_ky.unsigned_abs() as i64 > 2 * (ny / 2) {
            return Err(invalid("forcing_ky", format!(
                "forcing_ky {} lies beyond the 2/3-rule band of grid_ny {ny}", r.forcing_ky)));
        }
    }
    if r.forcing_type == "band" && r.forcing_amplitude != 0.0 {
        let kfl = r.forcing_kmax.floor();
        if 3.0 * kfl > 2.0 * (nx.min(ny) / 2) as f64 {
            return Err(invalid("forcing_kmax", format!(
                "forcing_kmax reaches beyond the 2/3-rule band of the {nx}x{ny} grid")));
        }
        // K3: the band must contain an integer mode (half-plane, kmin <= |k| <= kmax).
        let km = kfl as i64;
        let any = (0..=km.max(-1)).any(|ix| (-km..=km).any(|iy| {
            if ix == 0 && iy <= 0 { return false; }
            let kr = ((ix * ix + iy * iy) as f64).sqrt();
            kr >= r.forcing_kmin && kr <= r.forcing_kmax
        }));
        if !any {
            return Err(invalid("forcing_kmax",
                "the forcing band contains no integer mode: nothing would be forced"));
        }
    }
    // KH-VAL-2, K4: the analytic perturbation must exist on the grid, inside the 2/3-rule band.
    if r.initial_omega.is_none() && r.perturbation_amplitude != 0.0 {
        if r.perturbation_mode == 0 {
            return Err(invalid("perturbation_mode",
                "perturbation_mode 0 with a non-zero amplitude perturbs nothing"));
        }
        if 3 * r.perturbation_mode > 2 * (r.grid_nx / 2) {
            return Err(invalid("perturbation_mode", format!(
                "perturbation_mode {} lies beyond the 2/3-rule band of grid_nx {} (largest resolved mode {})",
                r.perturbation_mode, r.grid_nx, 2 * (r.grid_nx / 2) / 3)));
        }
    }
    // KH-VAL-2, K4b: initial_omega replaces the analytic field; an amplitude sent with it would be ignored.
    if r.initial_omega.is_some() && present(o, "perturbation_amplitude").is_some()
        && r.perturbation_amplitude != 0.0 {
        return Err(invalid("perturbation_amplitude",
            "perturbation_amplitude has no effect with initial_omega (send 0 or omit it)"));
    }

    // KH-VAL-2, K7: unknown keys are logged now and refused once no client sends one.
    const KNOWN: [&str; 21] = [
        "grid_nx", "grid_ny", "domain_lx", "domain_ly", "dt", "steps", "reynolds_number",
        "velocity_shear", "perturbation_amplitude", "perturbation_mode", "initial_omega",
        "init_omega", "forcing_type", "forcing_amplitude", "forcing_kx", "forcing_ky",
        "forcing_kmin", "forcing_kmax", "forcing_omega0", "linear_drag", "t_start"];
    for k in o.keys() {
        if !KNOWN.contains(&k.as_str()) {
            eprintln!("KH-VAL-2 K7: unknown key '{k}' ignored");
        }
    }
    Ok(r)
}

/// KH-VAL-2, K5: computed, but possibly not resolved. Same conditions, order and
/// prefixes in every backend and the oracle (kh_physics.request_warnings).
pub fn warnings(velocity_shear: f64, reynolds_number: f64, dt: f64,
                dx: f64, dy: f64, cfl_peak: f64) -> Vec<String> {
    let mut out = Vec::new();
    let nu = velocity_shear.abs() / reynolds_number;
    if cfl_peak > 1.0 {
        out.push(format!("cfl_advective: peak Courant number {} exceeds 1", g3(cfl_peak)));
    }
    let dmin = dx.min(dy);
    let diff = 2.0 * nu * dt / (dmin * dmin);
    if diff > 1.0 {
        out.push(format!("cfl_diffusive: 2*nu*dt/min(dx,dy)^2 = {} exceeds 1", g3(diff)));
    }
    if nu == 0.0 {
        out.push("inviscid: nu = 0, no physical dissipation; enstrophy accumulates at the 2/3 cutoff".to_string());
    }
    out
}

/// Up to three decimals. The prefix before ':' is the contract; the number's
/// formatting may differ between backends.
fn g3(x: f64) -> String {
    format!("{:.3}", x).trim_end_matches('0').trim_end_matches('.').to_string()
}
