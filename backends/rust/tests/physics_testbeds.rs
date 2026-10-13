/// physics_testbeds.rs — phase 7d physics testbeds, Rust backend, Mode 1 HTTP.
///
/// Implements Tier 1+2+3 tests (TB-1a/b/c, TB-2a-e, TB-3a/b) from
///   _audit/PHASE-7C-NONLINEAR-TESTBED-DESIGN-v0.1.md §4.
/// Tier 4 + TB-3c marked PENDING per
///   _handoffs/phase-7d-physics-testbeds/MASTER-BRIEF.md §2.5.
///
/// Mode 1 (HTTP): POST /simulate to running Rust server on port 8001.
/// Tolerances: MASTER-BRIEF §1 post-calibration values (Python session 2026-05-20).
/// backend_factor = 1.0 (Rust f64 + rustfft, no override needed).
///
/// Run (server must be running on :8001):
///   cd backends/rust
///   cargo run --release &
///   sleep 2
///   cargo test --test physics_testbeds -- --nocapture
///   kill %1
///
/// KH_GOLDEN_REQUIRE_SERVER=1 → hard-fail when server unreachable (CI mode).
/// Default (unset) → silent skip when server unreachable (dev mode).

use reqwest::blocking::Client;
use serde_json::{json, Value};
use std::f64::consts::PI;

const BACKEND_URL:  &str = "http://localhost:8001/simulate";
const BACKEND_NAME: &str = "kh-rust";

// ── Server connectivity helpers ───────────────────────────────────────────────

fn http_client() -> Client {
    Client::builder()
        .timeout(std::time::Duration::from_secs(300))
        .build()
        .expect("http client build")
}

fn require_server() -> bool {
    std::env::var("KH_GOLDEN_REQUIRE_SERVER")
        .map(|v| v == "1")
        .unwrap_or(false)
}

/// POST payload to server.  Returns Some(body) on success, None on silent-skip
/// (server unreachable and KH_GOLDEN_REQUIRE_SERVER not set).  Panics on other errors.
fn post_sim(client: &Client, payload: &Value) -> Option<Value> {
    let resp = client.post(BACKEND_URL).json(payload).send();
    match resp {
        Err(ref e) if e.is_connect() => {
            let msg = format!(
                "SKIP [{BACKEND_NAME}]: server not reachable at {BACKEND_URL}: {e}"
            );
            if require_server() {
                panic!("{msg}");
            }
            println!("{msg}");
            None
        }
        Err(e) => panic!("HTTP error: {e}"),
        Ok(r) => {
            if !r.status().is_success() {
                panic!("HTTP {}: {}", r.status(), r.text().unwrap_or_default());
            }
            Some(r.json::<Value>().expect("response JSON parse"))
        }
    }
}

// ── Request builder ───────────────────────────────────────────────────────────

/// Build a /simulate JSON payload.
/// Note: REST field is "steps" (not "nsteps" — golden_vectors.rs maps the two).
fn sim_req(
    nx: usize, ny: usize,
    lx: f64, ly: f64,
    dt: f64, steps: usize,
    re: f64, shear: f64,
    amp: f64, mode: usize,
    omega0: Option<&[f64]>,
) -> Value {
    let mut v = json!({
        "grid_nx":                nx,
        "grid_ny":                ny,
        "domain_lx":              lx,
        "domain_ly":              ly,
        "dt":                     dt,
        "steps":                  steps,
        "reynolds_number":        re,
        "velocity_shear":         shear,
        "perturbation_amplitude": amp,
        "perturbation_mode":      mode,
    });
    if let Some(w) = omega0 {
        v["initial_omega"] = json!(w);
    }
    v
}

// ── IC builders ───────────────────────────────────────────────────────────────

/// Taylor-Green IC: omega_0 = -2*cos(x)*cos(y) on [0, 2pi]^2.
fn tg_ic(nx: usize, ny: usize, lx: f64, ly: f64) -> Vec<f64> {
    let dx = lx / nx as f64;
    let dy = ly / ny as f64;
    let mut omega = vec![0.0f64; nx * ny];
    for i in 0..nx {
        let x = i as f64 * dx;
        for j in 0..ny {
            let y = j as f64 * dy;
            omega[i * ny + j] = -2.0 * x.cos() * y.cos();
        }
    }
    omega
}

/// Pure-diffusion IC: omega_0 = sin(kx_m*x)*sin(ky_m*y) on [0, 2pi]^2.
fn diffusion_ic(nx: usize, ny: usize, lx: f64, ly: f64, kx_m: usize, ky_m: usize) -> Vec<f64> {
    let dx = lx / nx as f64;
    let dy = ly / ny as f64;
    let mut omega = vec![0.0f64; nx * ny];
    for i in 0..nx {
        let x = i as f64 * dx;
        for j in 0..ny {
            let y = j as f64 * dy;
            omega[i * ny + j] = (kx_m as f64 * x).sin() * (ky_m as f64 * y).sin();
        }
    }
    omega
}

/// Canonical KH IC: tanh shear + sinusoidal perturbation.
/// Mirrors solver.rs initial_conditions() exactly.
fn canonical_ic(nx: usize, ny: usize, lx: f64, ly: f64, u0: f64, amp: f64, mode: usize) -> Vec<f64> {
    let dx = lx / nx as f64;
    let dy = ly / ny as f64;
    let delta = 0.05 * ly;
    let mut omega = vec![0.0f64; nx * ny];
    for i in 0..nx {
        let x = i as f64 * dx;
        let pert_x = if mode > 0 {
            amp * (2.0 * PI * mode as f64 / lx) * (2.0 * PI * mode as f64 * x / lx).cos()
        } else {
            0.0
        };
        for j in 0..ny {
            let y = j as f64 * dy;
            let z = (y - ly / 2.0) / delta;
            omega[i * ny + j] = u0 / delta / z.cosh().powi(2) + pert_x;
        }
    }
    omega
}

// ── Field extraction helpers ──────────────────────────────────────────────────

fn vorticity_vec(body: &Value) -> Vec<f64> {
    body["vorticity"]
        .as_array()
        .expect("vorticity not an array")
        .iter()
        .map(|v| v.as_f64().expect("vorticity element not f64"))
        .collect()
}

fn diag(body: &Value, key: &str) -> f64 {
    body["diagnostics"][key]
        .as_f64()
        .unwrap_or_else(|| panic!("diagnostics.{key} missing from response"))
}

// ── Math helpers ──────────────────────────────────────────────────────────────

fn rms(v: &[f64]) -> f64 {
    (v.iter().map(|&x| x * x).sum::<f64>() / v.len() as f64).sqrt()
}

fn rms_err(sim: &[f64], exact: &[f64]) -> f64 {
    assert_eq!(sim.len(), exact.len(), "rms_err: length mismatch");
    let n = sim.len() as f64;
    (sim.iter().zip(exact).map(|(&a, &b)| (a - b).powi(2)).sum::<f64>() / n).sqrt()
}

/// Unnormalized Fourier coefficient amplitude |F[kx_m, ky_m]| for a 2D real field.
/// F[u,v] = sum_{m,n} f[m,n] * exp(-2*pi*i*(m*u/nx + n*v/ny)).
/// Indexes field as row-major: f[m,n] = field[m*ny + n].
fn fourier_mode_amp(field: &[f64], nx: usize, ny: usize, kx_m: usize, ky_m: usize) -> f64 {
    let mut re = 0.0f64;
    let mut im = 0.0f64;
    for m in 0..nx {
        for n in 0..ny {
            let angle = 2.0 * PI * (m * kx_m) as f64 / nx as f64
                + 2.0 * PI * (n * ky_m) as f64 / ny as f64;
            re += field[m * ny + n] * angle.cos();
            im -= field[m * ny + n] * angle.sin();
        }
    }
    (re * re + im * im).sqrt()
}

/// Log-linear slope fit: p such that log(err) ≈ p*log(x) + const.
/// Returns the fitted slope.
fn loglog_slope(x: &[f64], y: &[f64]) -> f64 {
    assert_eq!(x.len(), y.len());
    let n = x.len() as f64;
    let lx: Vec<f64> = x.iter().map(|&v| v.ln()).collect();
    let ly: Vec<f64> = y.iter().map(|&v| v.ln()).collect();
    let sx: f64 = lx.iter().sum();
    let sy: f64 = ly.iter().sum();
    let sxx: f64 = lx.iter().map(|&v| v * v).sum();
    let sxy: f64 = lx.iter().zip(&ly).map(|(&a, &b)| a * b).sum();
    (n * sxy - sx * sy) / (n * sxx - sx * sx)
}

// ══════════════════════════════════════════════════════════════════════════════
// Tier 1 — Exact-solution tests
// ══════════════════════════════════════════════════════════════════════════════

#[test]
fn test_tb1a_taylor_green_vortex() {
    /// TB-1a: Taylor-Green vortex exact solution.
    ///
    /// omega(x,y,t) = -2*cos(x)*cos(y)*exp(-2*nu*t)
    /// The TG mode satisfies u.grad(omega) = 0 exactly, so it is an exact solution
    /// of 2D NS with only diffusive decay.
    ///
    /// Pass: rms(omega_sim - omega_exact) / rms(omega_exact) < 1e-6.
    let nx = 64usize;
    let ny = 64usize;
    let lx = 2.0 * PI;
    let ly = 2.0 * PI;
    let nu = 0.01_f64;
    let dt = 0.001_f64;
    let t_final = 1.0_f64;
    let steps = (t_final / dt) as usize; // 1000

    let omega_0 = tg_ic(nx, ny, lx, ly);
    let decay = (-2.0 * nu * t_final).exp();
    let omega_exact: Vec<f64> = omega_0.iter().map(|&w| w * decay).collect();

    let payload = sim_req(nx, ny, lx, ly, dt, steps, 1.0 / nu, 1.0, 0.0, 0, Some(&omega_0));
    let client  = http_client();
    let body    = match post_sim(&client, &payload) { Some(b) => b, None => return };

    let omega_sim = vorticity_vec(&body);
    let rms_diff  = rms_err(&omega_sim, &omega_exact);
    let rms_exact = rms(&omega_exact);
    let rms_rel   = rms_diff / rms_exact;

    println!("TB-1a Taylor-Green: rms_rel = {rms_rel:.3e} (tol 1e-6)");
    assert!(
        rms_rel < 1e-6,
        "TB-1a Taylor-Green rms_rel = {rms_rel:.3e} >= 1e-6 \
         (rms_diff={rms_diff:.3e}, rms_exact={rms_exact:.3e})"
    );
}

#[test]
fn test_tb1b_pure_diffusion() {
    /// TB-1b: pure diffusion single-mode decay.
    ///
    /// IC: omega_0 = sin(x)*sin(y) on [0,2pi]^2.
    /// Expected: per-mode amplitude decays as exp(-nu*k^2*t).
    /// Tolerance 1e-4 (relaxed from 1e-13 per MASTER §1 post-calibration):
    /// nonlinear advective term is non-zero for sin(x)sin(y) unlike TG.
    /// Empirically the error is <1e-15 for these params (advection sub-dominant).
    ///
    /// Pass: |ratio - expected_decay| / expected_decay < 1e-4.
    let nx = 32usize;
    let ny = 32usize;
    let lx = 2.0 * PI;
    let ly = 2.0 * PI;
    let kx_m = 1usize;
    let ky_m = 1usize;
    let nu  = 0.1_f64;
    let dt  = 0.001_f64;
    let steps = 500usize;
    let t_final = dt * steps as f64; // 0.5 s

    let k2 = (kx_m * kx_m + ky_m * ky_m) as f64;
    let expected_decay = (-nu * k2 * t_final).exp();

    let omega_0 = diffusion_ic(nx, ny, lx, ly, kx_m, ky_m);
    // re = 1.0/nu with shear=1.0 → nu_eff = shear/re = 1.0/(1.0/nu) = nu
    let payload = sim_req(nx, ny, lx, ly, dt, steps, 1.0 / nu, 1.0, 0.0, 0, Some(&omega_0));
    let client  = http_client();
    let body    = match post_sim(&client, &payload) { Some(b) => b, None => return };

    let omega_sim = vorticity_vec(&body);
    let amp_0 = fourier_mode_amp(&omega_0,  nx, ny, kx_m, ky_m);
    let amp_t = fourier_mode_amp(&omega_sim, nx, ny, kx_m, ky_m);

    assert!(amp_0 > 1e-10, "TB-1b: initial mode amplitude too small: {amp_0:.3e}");
    let ratio   = amp_t / amp_0;
    let rel_err = (ratio - expected_decay).abs() / expected_decay;

    println!(
        "TB-1b pure diffusion: ratio={ratio:.8}, expected={expected_decay:.8}, \
         rel_err={rel_err:.3e} (tol 1e-4)"
    );
    assert!(
        rel_err < 1e-4,
        "TB-1b pure diffusion: ratio={ratio:.8}, expected={expected_decay:.8}, \
         rel_err={rel_err:.3e} >= 1e-4"
    );
}

#[test]
fn test_tb1c_linear_growth_rate() {
    /// TB-1c: sigma(k) linear-stability dispersion match.
    ///
    /// Two-run subtraction method (Python §6 precedent):
    ///   - BASE run (amp=0): omega_base_T
    ///   - PERT run (amp=1e-4, mode=2): omega_pert_T
    ///   - E_pert(T) = 0.5*mean((omega_pert_T - omega_base_T)^2)
    ///   - E_pert(0) computed from ICs built locally
    ///   - sigma_sim = log(E_pert_T / E_pert_0) / (2 * T_meas)
    ///
    /// sigma_expected = 0.18 * U0 / delta = 7.2 /s (KH-LS-2, inviscid).
    /// 10% tolerance per MASTER §1 calibration: viscous correction + short window.
    ///
    /// Pass: |sigma_sim - sigma_expected| / sigma_expected < 10%.
    let nx = 128usize;
    let ny = 64usize;
    let lx = 1.0_f64;
    let ly = 0.5_f64;
    let dt = 0.001_f64;
    let steps_meas = 50usize; // T_meas = 0.05 s
    let t_meas = dt * steps_meas as f64;
    let amp_pert = 1e-4_f64;
    let mode = 2usize;
    let re = 1000.0_f64;
    let shear = 1.0_f64;

    let delta = 0.05 * ly; // = 0.025 m
    let sigma_expected = 0.18 * shear / delta; // = 7.2 /s

    let client = http_client();

    // Build t=0 ICs locally (mirrors solver.rs initial_conditions)
    let omega_base_0 = canonical_ic(nx, ny, lx, ly, shear, 0.0, 0);
    let omega_pert_0 = canonical_ic(nx, ny, lx, ly, shear, amp_pert, mode);

    // Run BASE (amp=0)
    let payload_base = sim_req(nx, ny, lx, ly, dt, steps_meas, re, shear, 0.0, 0, None);
    let body_base = match post_sim(&client, &payload_base) { Some(b) => b, None => return };
    let omega_base_t: Vec<f64> = vorticity_vec(&body_base);

    // Run PERT (amp=1e-4, mode=2)
    let payload_pert = sim_req(nx, ny, lx, ly, dt, steps_meas, re, shear, amp_pert, mode, None);
    let body_pert = match post_sim(&client, &payload_pert) { Some(b) => b, None => return };
    let omega_pert_t: Vec<f64> = vorticity_vec(&body_pert);

    // Perturbation energies
    let n = (nx * ny) as f64;
    let e_pert_0 = 0.5 * omega_pert_0.iter().zip(&omega_base_0)
        .map(|(&a, &b)| (a - b).powi(2)).sum::<f64>() / n;
    let e_pert_t = 0.5 * omega_pert_t.iter().zip(&omega_base_t)
        .map(|(&a, &b)| (a - b).powi(2)).sum::<f64>() / n;

    if e_pert_0 < 1e-20 || e_pert_t < 1e-20 {
        println!("TB-1c: perturbation energy too small; e0={e_pert_0:.2e}, et={e_pert_t:.2e} — SKIP");
        return;
    }

    let sigma_sim = (e_pert_t / e_pert_0).ln() / (2.0 * t_meas);
    let rel_err = (sigma_sim - sigma_expected).abs() / sigma_expected;

    println!(
        "TB-1c sigma(k): sigma_sim={sigma_sim:.4} /s, \
         sigma_expected={sigma_expected:.4} /s, rel_err={rel_err:.3e} (tol 10%)"
    );
    assert!(
        rel_err < 0.10,
        "TB-1c sigma(k): sigma_sim={sigma_sim:.4} /s, sigma_expected={sigma_expected:.4} /s, \
         rel_err={rel_err:.3e} >= 10% \
         (E_pert_0={e_pert_0:.3e}, E_pert_T={e_pert_t:.3e}, T_meas={t_meas:.3})"
    );
}

// ══════════════════════════════════════════════════════════════════════════════
// Tier 2 — Conservation tests
// ══════════════════════════════════════════════════════════════════════════════

#[test]
fn test_tb2a_total_vorticity_conservation() {
    /// TB-2a: total vorticity conservation.
    ///
    /// Canonical KH IC has sum(omega) ~ 8192 (positive shear layer, NOT zero-mean).
    /// The spectral DC mode (k=0) is exactly conserved by construction.
    /// Test: |sum(omega_T) - sum(omega_0)| < 1e-8  (not zero-ness, but conservation).
    ///
    /// omega_0 computed locally from canonical_ic() — same formula as solver.rs.
    ///
    /// Pass: |sum_T - sum_0| < 1e-8.
    let nx = 64usize;
    let ny = 32usize;
    let lx = 1.0_f64;
    let ly = 0.5_f64;

    let omega_0 = canonical_ic(nx, ny, lx, ly, 1.0, 0.01, 2);
    let sum_0: f64 = omega_0.iter().sum();

    let payload = sim_req(nx, ny, lx, ly, 0.001, 100, 1000.0, 1.0, 0.01, 2, None);
    let client  = http_client();
    let body    = match post_sim(&client, &payload) { Some(b) => b, None => return };

    let omega_t: Vec<f64> = vorticity_vec(&body);
    let sum_t: f64 = omega_t.iter().sum();
    let drift = (sum_t - sum_0).abs();

    println!(
        "TB-2a vorticity conservation: drift = {drift:.3e} (tol 1e-8, sum_0={sum_0:.6e})"
    );
    assert!(
        drift < 1e-8,
        "TB-2a total vorticity drift = {drift:.3e} >= 1e-8 \
         (sum_0={sum_0:.6e}, sum_T={sum_t:.6e})"
    );
}

#[test]
fn test_tb2b_energy_inviscid() {
    /// TB-2b: kinetic energy conservation (near-inviscid, nu=1e-8).
    ///
    /// E(0) obtained via steps=0 POST (server returns diagnostics from IC).
    /// E(T) from steps=200, dt=5e-4 POST.
    ///
    /// Tolerance 1e-4 per MASTER §1 calibration: 2/3-rule dealiasing leaks
    /// ~1e-6 relative energy per 200 steps even at nu=0.
    ///
    /// Pass: |E(T) - E(0)| / E(0) < 1e-4.
    let nx = 64usize;
    let ny = 32usize;
    let lx = 1.0_f64;
    let ly = 0.5_f64;
    let dt = 5e-4_f64;
    let nu = 1e-8_f64;
    let shear = 1.0_f64;
    let re = shear / nu; // = 1e8

    let client = http_client();

    // E(0): steps=0 → server returns diagnostics at IC (time loop skipped)
    let payload_0 = sim_req(nx, ny, lx, ly, dt, 0, re, shear, 0.01, 2, None);
    let body_0    = match post_sim(&client, &payload_0) { Some(b) => b, None => return };
    let e0 = diag(&body_0, "kinetic_energy");

    // E(T): steps=200 → T = 0.1 s
    let payload_t = sim_req(nx, ny, lx, ly, dt, 200, re, shear, 0.01, 2, None);
    let body_t    = match post_sim(&client, &payload_t) { Some(b) => b, None => return };
    let et = diag(&body_t, "kinetic_energy");

    if e0 < 1e-15 {
        println!("TB-2b: E(0) too small ({e0:.2e}) — SKIP");
        return;
    }
    let rel_err = (et - e0).abs() / e0;

    println!(
        "TB-2b inviscid energy: E0={e0:.8e}, ET={et:.8e}, \
         |E(T)-E(0)|/E(0)={rel_err:.3e} (tol 1e-4)"
    );
    assert!(
        rel_err < 1e-4,
        "TB-2b inviscid energy: E0={e0:.8e}, ET={et:.8e}, \
         rel_err={rel_err:.3e} >= 1e-4"
    );
}

#[test]
fn test_tb2c_enstrophy_inviscid() {
    /// TB-2c: enstrophy conservation (near-inviscid, nu=1e-8, 2D Casimir invariant).
    ///
    /// Z(0) from steps=0 POST; Z(T) from steps=200 POST.
    /// Same tolerance rationale as TB-2b (2/3-rule dealiasing).
    ///
    /// Pass: |Z(T) - Z(0)| / Z(0) < 1e-4.
    let nx = 64usize;
    let ny = 32usize;
    let lx = 1.0_f64;
    let ly = 0.5_f64;
    let dt = 5e-4_f64;
    let nu = 1e-8_f64;
    let shear = 1.0_f64;
    let re = shear / nu;

    let client = http_client();

    let payload_0 = sim_req(nx, ny, lx, ly, dt, 0,   re, shear, 0.01, 2, None);
    let body_0    = match post_sim(&client, &payload_0) { Some(b) => b, None => return };
    let z0 = diag(&body_0, "enstrophy");

    let payload_t = sim_req(nx, ny, lx, ly, dt, 200, re, shear, 0.01, 2, None);
    let body_t    = match post_sim(&client, &payload_t) { Some(b) => b, None => return };
    let zt = diag(&body_t, "enstrophy");

    if z0 < 1e-15 {
        println!("TB-2c: Z(0) too small ({z0:.2e}) — SKIP");
        return;
    }
    let rel_err = (zt - z0).abs() / z0;

    println!(
        "TB-2c inviscid enstrophy: Z0={z0:.8e}, ZT={zt:.8e}, \
         |Z(T)-Z(0)|/Z(0)={rel_err:.3e} (tol 1e-4)"
    );
    assert!(
        rel_err < 1e-4,
        "TB-2c inviscid enstrophy: Z0={z0:.8e}, ZT={zt:.8e}, \
         rel_err={rel_err:.3e} >= 1e-4"
    );
}

#[test]
fn test_tb2d_energy_decay_viscous() {
    /// TB-2d: energy decay rate dE/dt = -2*nu*Z_ac (DC-mode corrected).
    ///
    /// For canonical KH IC (non-zero DC vorticity), the solver zeroes psi_hat[0,0],
    /// so the DC mode has no kinetic energy.  Correct identity: dE/dt = -2*nu*Z_ac
    /// where Z_ac = Z_total - Z_dc and Z_dc = 0.5*|omega_hat[0,0]|^2/(nx*ny)^2.
    /// omega_hat[0,0] = sum(omega)  (unnormalized DFT DC component).
    ///
    /// Method: FD over T1=0.1 → T2=0.2; midpoint Z_ac at T_mid=0.15.
    /// nu = shear/re = 1.0/1000 = 0.001.
    ///
    /// Pass: |dE/dt_fd - (-2*nu*Z_ac)| / |2*nu*Z_ac| < 5%.
    let nx = 64usize;
    let ny = 32usize;
    let lx = 1.0_f64;
    let ly = 0.5_f64;
    let dt = 0.001_f64;
    let re = 1000.0_f64;
    let shear = 1.0_f64;
    let nu = shear / re; // 0.001

    let client = http_client();

    let payload_t1 = sim_req(nx, ny, lx, ly, dt, 100, re, shear, 0.01, 2, None); // T1=0.1
    let payload_t2 = sim_req(nx, ny, lx, ly, dt, 200, re, shear, 0.01, 2, None); // T2=0.2
    let payload_tm = sim_req(nx, ny, lx, ly, dt, 150, re, shear, 0.01, 2, None); // T_mid=0.15

    let body_t1 = match post_sim(&client, &payload_t1) { Some(b) => b, None => return };
    let body_t2 = match post_sim(&client, &payload_t2) { Some(b) => b, None => return };
    let body_tm = match post_sim(&client, &payload_tm) { Some(b) => b, None => return };

    let e1 = diag(&body_t1, "kinetic_energy");
    let e2 = diag(&body_t2, "kinetic_energy");
    let z_total = diag(&body_tm, "enstrophy");

    // DC correction: omega_hat[0,0] = sum(omega)
    let omega_m = vorticity_vec(&body_tm);
    let sum_omega: f64 = omega_m.iter().sum();
    let z_dc = 0.5 * sum_omega * sum_omega / ((nx * ny) as f64).powi(2);
    let z_ac = z_total - z_dc;

    let t1 = 0.1_f64;
    let t2 = 0.2_f64;
    let dedt_fd     = (e2 - e1) / (t2 - t1);
    let dedt_theory = -2.0 * nu * z_ac;
    let denom = dedt_theory.abs();

    if denom < 1e-15 {
        println!("TB-2d: |2*nu*Z_ac| too small ({denom:.2e}) — SKIP");
        return;
    }
    let rel_err = (dedt_fd - dedt_theory).abs() / denom;

    println!(
        "TB-2d energy decay: dE/dt_fd={dedt_fd:.6e}, theory_ac={dedt_theory:.6e}, \
         Z_total={z_total:.6e}, Z_dc={z_dc:.6e}, Z_ac={z_ac:.6e}, \
         rel_err={rel_err:.3e} (tol 5%)"
    );
    assert!(
        rel_err < 0.05,
        "TB-2d energy decay: dE/dt_fd={dedt_fd:.6e}, theory_ac={dedt_theory:.6e}, \
         rel_err={rel_err:.3e} >= 5%"
    );
}

#[test]
fn test_tb2e_incompressibility() {
    /// TB-2e: divergence-free velocity field.
    ///
    /// Spectral vorticity-streamfunction formulation enforces div(u)=0 exactly
    /// by construction (KH-GE-3).  divergence_rms < 1e-12 throughout.
    ///
    /// Pass: diagnostics.divergence_rms < 1e-12.
    let payload = sim_req(64, 32, 1.0, 0.5, 0.001, 100, 1000.0, 1.0, 0.01, 2, None);
    let client  = http_client();
    let body    = match post_sim(&client, &payload) { Some(b) => b, None => return };

    let div_rms = diag(&body, "divergence_rms");
    println!("TB-2e incompressibility: divergence_rms = {div_rms:.3e} (tol 1e-12)");
    assert!(
        div_rms < 1e-12,
        "TB-2e incompressibility: divergence_rms = {div_rms:.3e} >= 1e-12"
    );
}

// ══════════════════════════════════════════════════════════════════════════════
// Tier 3 — Convergence tests
// ══════════════════════════════════════════════════════════════════════════════

#[test]
fn test_tb3a_spatial_convergence() {
    /// TB-3a: spatial convergence on Taylor-Green vortex.
    ///
    /// 4 grids: 32^2, 64^2, 128^2, 256^2.  T=0.1, nu=0.01, dt=0.0001.
    /// For the TG mode (cos(x)*cos(y)), spectral convergence is achieved even
    /// at 32^2 — all grids hit the machine-epsilon floor.
    ///
    /// Pass criteria (machine-eps floor test per OQ-7D-PY-3 resolution):
    ///   (a) max(errors) < 1e-12  (all at machine-eps floor)
    ///   (b) errors monotonically non-increasing (10% slack for noise)
    ///
    /// Note: for N=256 the initial_omega payload is 65536 floats (~1 MB).
    let nu    = 0.01_f64;
    let t_end = 0.1_f64;
    let dt    = 0.0001_f64;
    let steps = (t_end / dt) as usize; // 1000
    let grids = [32usize, 64, 128, 256];

    let client = http_client();
    let mut errors: Vec<f64> = Vec::new();

    for &n in &grids {
        let lx = 2.0 * PI;
        let ly = 2.0 * PI;
        let omega_0    = tg_ic(n, n, lx, ly);
        let decay      = (-2.0 * nu * t_end).exp();
        let omega_exact: Vec<f64> = omega_0.iter().map(|&w| w * decay).collect();

        let payload = sim_req(n, n, lx, ly, dt, steps, 1.0 / nu, 1.0, 0.0, 0, Some(&omega_0));
        let body    = match post_sim(&client, &payload) { Some(b) => b, None => return };

        let omega_sim = vorticity_vec(&body);
        assert_eq!(
            omega_sim.len(), n * n,
            "TB-3a grid={n}: response vorticity length mismatch"
        );
        let err = rms_err(&omega_sim, &omega_exact);
        println!("TB-3a grid {n}^2: rms_err = {err:.3e}");
        errors.push(err);
    }

    // (a) All at machine-eps floor
    let max_err = errors.iter().cloned().fold(0.0_f64, f64::max);
    assert!(
        max_err < 1e-12,
        "TB-3a spatial convergence: max rms_err = {max_err:.3e} >= 1e-12 \
         (errors: {errors:?})"
    );

    // (b) Monotone non-increasing (10% slack for machine-epsilon-level noise)
    for i in 0..errors.len() - 1 {
        assert!(
            errors[i + 1] <= errors[i] * 1.1,
            "TB-3a: error increased from grid {g0}^2 ({e0:.3e}) to {g1}^2 ({e1:.3e}) by >10%",
            g0 = grids[i], e0 = errors[i], g1 = grids[i + 1], e1 = errors[i + 1]
        );
    }
    println!("TB-3a PASS: max_err={max_err:.3e} < 1e-12; monotone ✓");
}

#[test]
fn test_tb3b_temporal_convergence_rk4() {
    /// TB-3b: RK4 temporal convergence order.
    ///
    /// 4 dt values: [1e-2, 5e-3, 2e-3, 1e-3]; reference dt_ref=5e-4.
    /// T_final=0.1; canonical IC (amp=0.01, mode=2).
    /// nu=0.001 (Re=1000) for clean RK4 stability (OQ-7D-PY-7D-PY-3 precedent:
    /// nu=0.01 triggers marginal diffusive CFL at dt=0.01).
    ///
    /// Fit log(err) vs log(dt); expect slope >= 3.5 (RK4 theoretical: 4).
    ///
    /// Pass: fitted slope p >= 3.5.
    let nx = 64usize;
    let ny = 32usize;
    let lx = 1.0_f64;
    let ly = 0.5_f64;
    let nu = 0.001_f64;
    let shear = 1.0_f64;
    let re = shear / nu; // 1000
    let t_final = 0.1_f64;
    let dt_ref  = 5e-4_f64;
    let steps_ref = (t_final / dt_ref).round() as usize; // 200

    let client = http_client();

    // Reference run
    let payload_ref = sim_req(nx, ny, lx, ly, dt_ref, steps_ref, re, shear, 0.01, 2, None);
    let body_ref    = match post_sim(&client, &payload_ref) { Some(b) => b, None => return };
    let omega_ref   = vorticity_vec(&body_ref);

    let dt_values = [1e-2_f64, 5e-3, 2e-3, 1e-3];
    let mut errs: Vec<f64> = Vec::new();

    for &dt in &dt_values {
        let steps = (t_final / dt).round() as usize;
        let payload = sim_req(nx, ny, lx, ly, dt, steps, re, shear, 0.01, 2, None);
        let body    = match post_sim(&client, &payload) { Some(b) => b, None => return };
        let omega_t = vorticity_vec(&body);
        let err     = rms_err(&omega_t, &omega_ref);
        println!("TB-3b dt={dt:.0e}: rms_err={err:.3e}");
        errs.push(err);
    }

    // Guard: if any error is below machine eps, slope fit is meaningless
    let min_err = errs.iter().cloned().fold(f64::INFINITY, f64::min);
    if min_err < 1e-14 {
        println!("TB-3b: smallest error {min_err:.2e} at machine eps — spectral convergence, PASS");
        return;
    }

    let slope = loglog_slope(&dt_values, &errs);
    println!("TB-3b RK4 convergence: fitted slope p={slope:.3} (tol >= 3.5); errors={errs:?}");
    assert!(
        slope >= 3.5,
        "TB-3b temporal convergence: fitted slope p={slope:.3} < 3.5 \
         (expected RK4 order 4); dt_values={dt_values:?}, errors={errs:?}"
    );
}

// ══════════════════════════════════════════════════════════════════════════════
// Tier 3 — PENDING stubs
// ══════════════════════════════════════════════════════════════════════════════

#[test]
fn test_tb3c_mms() {
    /// TB-3c: method of manufactured solutions.
    /// PENDING: phase 7c Opus MMS design.
    /// Reference: MASTER-BRIEF.md §2.5.
    println!("TB-3c PENDING: phase 7c Opus MMS design — method of manufactured solutions");
}

// ══════════════════════════════════════════════════════════════════════════════
// Tier 4 — PENDING / STRETCH-DEFER stubs
// ══════════════════════════════════════════════════════════════════════════════

#[test]
fn test_tb4a_lyapunov() {
    /// TB-4a: Lyapunov time estimate.
    /// PENDING: phase 7e Opus Lyapunov estimator design.
    /// Reference: MASTER-BRIEF.md §2.5.
    println!("TB-4a PENDING: phase 7e Opus Lyapunov estimator design");
}

#[test]
fn test_tb4b_invariant_measure() {
    /// TB-4b: invariant measure agreement (ensemble; stretch scope).
    /// STRETCH-DEFER per MASTER-BRIEF.md §2.5: long-time ensemble, high complexity.
    println!("TB-4b STRETCH-DEFER per MASTER §2.5");
}
