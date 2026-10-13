/// validation_test.rs — KH-003 acceptance criteria
///
/// Validates Rust diagnostics against the oracle golden vector V1
/// (kh-sim/shared/golden-vectors/v1-canonical-reference.json, class A2), which
/// replaced the retired kh-sim/shared/physics/kh_reference_output.json on
/// 2026-10-08. The request is read from V1's "parameters"; the tolerance is
/// V1's A2 band: |got - ref| <= max(scalar_abs * f, scalar_rel * f * |ref|),
/// f = the "kh-rust" scalar_factor override (1 when absent).
/// Incompressibility: divergence_rms < 1e-10.
///
/// Run: cargo test --test validation_test
use std::path::PathBuf;

use serde_json::Value;

const BACKEND: &str = "kh-rust";

// ── Load golden vector V1 ─────────────────────────────────────────────────────

fn load_v1() -> Value {
    // Path relative to the crate root (kh-sim/backends/rust/)
    let mut path = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    path.push("../../shared/golden-vectors/v1-canonical-reference.json");

    let raw = std::fs::read_to_string(&path)
        .unwrap_or_else(|_| panic!("Cannot open golden vector at {:?}", path));
    serde_json::from_str(&raw).expect("Failed to parse golden vector JSON")
}

fn f64_at(v: &Value, key: &str) -> f64 {
    v[key].as_f64().unwrap_or_else(|| panic!("golden vector: {} is not a number", key))
}

fn usize_at(v: &Value, key: &str) -> usize {
    v[key].as_u64().unwrap_or_else(|| panic!("golden vector: {} is not an integer", key)) as usize
}

// ── Tests ─────────────────────────────────────────────────────────────────────

#[test]
fn test_diagnostics_within_tolerance() {
    let v1 = load_v1();
    let p = &v1["parameters"];
    let ref_d = &v1["diagnostics"];
    let tols = &v1["tolerances"];
    let factor = tols["backend_overrides"][BACKEND]["scalar_factor"].as_f64().unwrap_or(1.0);
    let s_abs = f64_at(tols, "scalar_abs") * factor;
    let s_rel = f64_at(tols, "scalar_rel") * factor;

    use kh_sim_rust::physics::simulate;
    use kh_sim_rust::models::SimulationRequest;

    let req = SimulationRequest {
        grid_nx:               usize_at(p, "grid_nx"),
        grid_ny:               usize_at(p, "grid_ny"),
        domain_lx:             f64_at(p, "domain_lx"),
        domain_ly:             f64_at(p, "domain_ly"),
        dt:                    f64_at(p, "dt"),
        steps:                 usize_at(p, "nsteps"),
        reynolds_number:       f64_at(p, "reynolds_number"),
        velocity_shear:        f64_at(p, "velocity_shear"),
        perturbation_amplitude:f64_at(p, "perturbation_amplitude"),
        perturbation_mode:     usize_at(p, "perturbation_mode"),
        initial_omega:         None,
        // Phase-7f fields; values = forcing off.
        forcing_type:          "none".to_string(),
        forcing_amplitude:     0.0,
        forcing_kx:            1,
        forcing_ky:            1,
        forcing_kmin:          0.0,
        forcing_kmax:          0.0,
        forcing_omega0:        0.0,
        linear_drag:           0.0,
        t_start:               0.0,
    };

    let result = simulate(&req);
    let got = &result.diagnostics;

    println!("--- Rust vs golden vector V1 (A2 band, factor {}) ---", factor);
    let mut failures = Vec::new();
    for (name, value) in [
        ("kinetic_energy", got.kinetic_energy),
        ("enstrophy",      got.enstrophy),
        ("max_vorticity",  got.max_vorticity),
    ] {
        let reference = f64_at(ref_d, name);
        let diff = (value - reference).abs();
        let band = s_abs.max(s_rel * reference.abs());
        println!("{:<15}: got={:.17}  ref={:.17}  |diff|={:.3e}  band={:.3e}",
                 name, value, reference, diff, band);
        if !(diff <= band) {
            failures.push(format!("{}: |diff| {:.3e} > band {:.3e}", name, diff, band));
        }
    }
    println!("divergence_rms : got={:.2e}", got.divergence_rms);
    if !(got.divergence_rms < 1e-10) {
        failures.push(format!("divergence_rms {:.2e} >= 1e-10", got.divergence_rms));
    }
    assert!(failures.is_empty(), "outside the A2 band of V1: {}", failures.join("; "));
}

#[test]
fn test_field_shapes() {
    use kh_sim_rust::physics::simulate;
    use kh_sim_rust::models::SimulationRequest;

    let req = SimulationRequest {
        grid_nx: 32, grid_ny: 16,
        domain_lx: 1.0, domain_ly: 0.5,
        dt: 0.001, steps: 10,
        reynolds_number: 1000.0, velocity_shear: 1.0,
        perturbation_amplitude: 0.01, perturbation_mode: 2,
        initial_omega: None,
        forcing_type: "none".to_string(), forcing_amplitude: 0.0,
        forcing_kx: 1, forcing_ky: 1, forcing_kmin: 0.0, forcing_kmax: 0.0,
        forcing_omega0: 0.0, linear_drag: 0.0, t_start: 0.0,
    };

    let result = simulate(&req);
    let n = req.grid_nx * req.grid_ny;
    assert_eq!(result.u_velocity.len(), n, "u_velocity shape mismatch");
    assert_eq!(result.v_velocity.len(), n, "v_velocity shape mismatch");
    assert_eq!(result.vorticity.len(),  n, "vorticity shape mismatch");
    assert_eq!(result.pressure.len(),   n, "pressure shape mismatch");
    assert!((result.t_final - 0.01).abs() < 1e-9, "t_final mismatch");
}
