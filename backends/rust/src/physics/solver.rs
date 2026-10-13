/// solver.rs — KH instability physics kernel (Rust port of kh_physics.py)
///
/// Physics: 2D incompressible Navier-Stokes, vorticity-streamfunction form.
/// Numerics: pseudo-spectral (full-complex 2D FFT), RK4 time integration.
///
/// Reference: kh-sim/shared/physics/KH-PHYSICS.md
/// Canonical: kh-sim/shared/physics/kh_physics.py
use std::f64::consts::PI;

use num_complex::Complex64;

use super::fft2d::{Fft2DPlans, make_wavenumbers};
use crate::models::{Diagnostics, SimulationRequest, SimulationResult};

// -- Phase 7f forcing (PHASE-7F-FORCING-API-DESIGN-v0.1.md sec.3 v0.1 + sec.7 v0.2) --
//
// Time-independent "mode" forcing and (optionally rotating-phase) "band" forcing,
// constructed by a direct double-loop with deterministic, RNG-free phases so the
// field is bit-near-identical across all six backends. Linear (Ekman) drag is a
// separate -alpha*omega RHS term, threaded alongside.
//
// Grid points: x_i = i*lx/nx (i=0..nx-1), y_j = j*ly/ny (j=0..ny-1).
// Arrays indexed (i,j) = (x,y), layout idx = i*ny + j (matches omega).

/// Forcing operator: precomputed once before the time loop, evaluated per RK4
/// stage at the stage time. None => no forcing (pre-7f RHS bit-identical).
enum Forcing {
    /// Steady Taylor-Green-mode field (omega0 ignored -- CP-7F-D needs steadiness).
    Static(Vec<f64>),
    /// Rotating-phase isotropic band: F(t) = (A/sqrt(M)) * sum_m
    ///   [ C_m*cos(Om_m*t) - S_m*sin(Om_m*t) ], with C_m/S_m the per-mode spatial
    /// cos/sin fields and Om_m the per-mode rotation rate (v0.2 sec.7).
    Band {
        c_modes: Vec<Vec<f64>>, // C_m = (A/sqrt(M)) * cos(arg0_m(x,y))
        s_modes: Vec<Vec<f64>>, // S_m = (A/sqrt(M)) * sin(arg0_m(x,y))
        omega_m: Vec<f64>,      // Om_m, per-mode rotation rate
        n: usize,
    },
}

impl Forcing {
    /// Evaluate the forcing field at time t. Returns a (nx*ny) vector.
    fn eval(&self, t: f64) -> Vec<f64> {
        match self {
            Forcing::Static(f) => f.clone(),
            Forcing::Band { c_modes, s_modes, omega_m, n } => {
                let mut out = vec![0.0f64; *n];
                for m in 0..omega_m.len() {
                    let ph = omega_m[m] * t;
                    let (c, s) = (ph.cos(), ph.sin());
                    let cm = &c_modes[m];
                    let sm = &s_modes[m];
                    for idx in 0..*n {
                        // C_m*cos(Om*t) - S_m*sin(Om*t)
                        out[idx] += cm[idx] * c - sm[idx] * s;
                    }
                }
                out
            }
        }
    }
}

/// Build the forcing operator from the request. Returns None when forcing is off
/// (type "none" or amplitude 0) so the no-forcing RHS path is taken unchanged.
fn build_forcing(req: &SimulationRequest, nx: usize, ny: usize, lx: f64, ly: f64) -> Option<Forcing> {
    let ftype = req.forcing_type.as_str();
    let amp = req.forcing_amplitude;
    if ftype == "none" || amp == 0.0 {
        return None;
    }

    let dx = lx / nx as f64;
    let dy = ly / ny as f64;

    match ftype {
        "mode" => {
            // F(x,y) = A * cos(2*pi*kx*x/lx) * cos(2*pi*ky*y/ly)   (steady)
            let kx = req.forcing_kx as f64;
            let ky = req.forcing_ky as f64;
            let mut f = vec![0.0f64; nx * ny];
            for i in 0..nx {
                let x = i as f64 * dx;
                let cx = (2.0 * PI * kx * x / lx).cos();
                for j in 0..ny {
                    let y = j as f64 * dy;
                    let cy = (2.0 * PI * ky * y / ly).cos();
                    f[i * ny + j] = amp * cx * cy;
                }
            }
            Some(Forcing::Static(f))
        }
        "band" => {
            // Isotropic shell of integer mode pairs (ix,iy) in the half-plane
            // (ix >= 0; skip ix==0 && iy<=0) with kmin <= sqrt(ix^2+iy^2) <= kmax.
            // Requires lx == ly (integer-mode-radius isotropy).
            assert!(
                (lx - ly).abs() < 1e-12,
                "band forcing requires lx == ly (got lx={lx}, ly={ly})"
            );
            let kmin = req.forcing_kmin;
            let kmax = req.forcing_kmax;
            let omega0 = req.forcing_omega0;

            // Enumerate contributing modes.
            let kmax_i = kmax.floor() as i64;
            let mut modes: Vec<(i64, i64)> = Vec::new();
            for ix in 0..=kmax_i {
                for iy in -kmax_i..=kmax_i {
                    if ix == 0 && iy <= 0 {
                        continue;
                    }
                    let kr = ((ix * ix + iy * iy) as f64).sqrt();
                    if kr >= kmin && kr <= kmax {
                        modes.push((ix, iy));
                    }
                }
            }
            let m_count = modes.len();
            if m_count == 0 {
                return None;
            }
            let scale = amp / (m_count as f64).sqrt();

            // Deterministic, RNG-free phases (sec.3 + sec.7):
            //   phi0(ix,iy) = 2*pi*frac(ix*sqrt(2) + iy*sqrt(3))
            //   Om(ix,iy)   = 2*pi*frac(ix*sqrt(11) + iy*sqrt(13)) * omega0
            let frac = |v: f64| v - v.floor();
            let s2 = 2.0_f64.sqrt();
            let s3 = 3.0_f64.sqrt();
            let s11 = 11.0_f64.sqrt();
            let s13 = 13.0_f64.sqrt();

            let mut c_modes: Vec<Vec<f64>> = Vec::with_capacity(m_count);
            let mut s_modes: Vec<Vec<f64>> = Vec::with_capacity(m_count);
            let mut omega_m: Vec<f64> = Vec::with_capacity(m_count);

            for &(ix, iy) in &modes {
                let phi0 = 2.0 * PI * frac(ix as f64 * s2 + iy as f64 * s3);
                let om = 2.0 * PI * frac(ix as f64 * s11 + iy as f64 * s13) * omega0;
                let mut cm = vec![0.0f64; nx * ny];
                let mut sm = vec![0.0f64; nx * ny];
                for i in 0..nx {
                    let x = i as f64 * dx;
                    for j in 0..ny {
                        let y = j as f64 * dy;
                        // arg0 = 2*pi*(ix*x/lx + iy*y/ly) + phi0
                        let arg0 = 2.0 * PI * (ix as f64 * x / lx + iy as f64 * y / ly) + phi0;
                        cm[i * ny + j] = scale * arg0.cos();
                        sm[i * ny + j] = scale * arg0.sin();
                    }
                }
                c_modes.push(cm);
                s_modes.push(sm);
                omega_m.push(om);
            }
            Some(Forcing::Band { c_modes, s_modes, omega_m, n: nx * ny })
        }
        other => panic!("unknown forcing_type: {other:?}"),
    }
}

// ── Initial conditions ────────────────────────────────────────────────────────

/// Analytic initial vorticity: omega = dv/dx - du/dy
///   u(x,y) = U0 * tanh((y - Ly/2) / delta)
///   v(x,y) = amp * sin(2*PI*mode*x / Lx)
/// => omega_shear = U0/delta * sech^2((y-Ly/2)/delta)
///    omega_pert  = amp * (2*PI*mode/Lx) * cos(2*PI*mode*x/Lx)
#[allow(clippy::too_many_arguments)]
fn initial_conditions(
    nx: usize,
    ny: usize,
    lx: f64,
    ly: f64,
    u0: f64,
    delta: f64,
    amp: f64,
    mode: usize,
) -> Vec<f64> {
    let dx = lx / nx as f64;
    let dy = ly / ny as f64;
    let mut omega = vec![0.0f64; nx * ny];

    for i in 0..nx {
        let x = i as f64 * dx;
        let pert_x = amp * (2.0 * PI * mode as f64 / lx) * (2.0 * PI * mode as f64 * x / lx).cos();

        for j in 0..ny {
            let y = j as f64 * dy;
            let z = (y - ly / 2.0) / delta;
            let shear = u0 / delta / z.cosh().powi(2); // sech^2 = 1/cosh^2
            omega[i * ny + j] = shear + pert_x;
        }
    }
    omega
}

// ── Spectral Poisson solve ────────────────────────────────────────────────────

/// Solve nabla^2 psi = -omega in Fourier space:
///   psi_hat = omega_hat / (kx^2 + ky^2)
///   psi_hat[0,0] = 0  (zero mean streamfunction)
fn solve_poisson(
    omega_hat: &[Complex64],
    kx: &[f64],
    ky: &[f64],
    nx: usize,
    ny: usize,
) -> Vec<Complex64> {
    let mut psi_hat = omega_hat.to_vec();

    for i in 0..nx {
        for j in 0..ny {
            let idx = i * ny + j;
            let k2 = kx[idx].powi(2) + ky[idx].powi(2);
            if i == 0 && j == 0 {
                psi_hat[idx] = Complex64::new(0.0, 0.0); // zero mean
            } else {
                psi_hat[idx] = omega_hat[idx] / k2;
            }
        }
    }
    psi_hat
}

// ── Velocity from streamfunction ──────────────────────────────────────────────

/// u = d(psi)/dy  =>  u_hat = i * ky * psi_hat
/// v = -d(psi)/dx =>  v_hat = -i * kx * psi_hat
fn velocity_from_psi_hat(
    psi_hat: &[Complex64],
    kx: &[f64],
    ky: &[f64],
    plans: &Fft2DPlans,
    nx: usize,
    ny: usize,
) -> (Vec<f64>, Vec<f64>) {
    let i_unit = Complex64::new(0.0, 1.0);

    let u_hat: Vec<Complex64> = (0..nx * ny)
        .map(|idx| i_unit * ky[idx] * psi_hat[idx])
        .collect();

    let v_hat: Vec<Complex64> = (0..nx * ny)
        .map(|idx| -i_unit * kx[idx] * psi_hat[idx])
        .collect();

    let u = plans.ifft2_real(&u_hat);
    let v = plans.ifft2_real(&v_hat);
    (u, v)
}

// ── 2/3-rule dealiasing + CFL bound ──────────────────────────────────────────

/// KH-NM-4: 2/3-rule dealiasing (Orszag 1971). Zeros modes outside 2N/3 cutoff in-place.
fn dealias_two_thirds(field_hat: &mut [Complex64], kx: &[f64], ky: &[f64], nx: usize, ny: usize) {
    let kx_max = kx.iter().map(|k| k.abs()).fold(0.0_f64, f64::max);
    let ky_max = ky.iter().map(|k| k.abs()).fold(0.0_f64, f64::max);
    let kx_cut = (2.0 / 3.0) * kx_max;
    let ky_cut = (2.0 / 3.0) * ky_max;
    for idx in 0..nx * ny {
        if kx[idx].abs() > kx_cut || ky[idx].abs() > ky_cut {
            field_hat[idx] = Complex64::new(0.0, 0.0);
        }
    }
}

/// KH-NM-5: CFL-limited time step. Returns min of advective and diffusive bounds.
fn cfl_dt_bound(u_max: f64, v_max: f64, dx: f64, dy: f64, nu: f64) -> f64 {
    let dt_adv  = (dx / (u_max + 1e-30)).min(dy / (v_max + 1e-30));
    // nu = 0 (inviscid) is a valid limit: no diffusive bound (modelling rules v0.1, D).
    let dt_diff = if nu == 0.0 { f64::INFINITY } else { (dx * dx / (2.0 * nu)).min(dy * dy / (2.0 * nu)) };
    dt_adv.min(dt_diff)
}

// ── Vorticity RHS ─────────────────────────────────────────────────────────────

/// f(omega) = -u * d(omega)/dx - v * d(omega)/dy + nu * laplacian(omega)
///            [- alpha*omega] [+ F(x,y,t)]   (phase 7f optional terms)
///
/// `drag` (alpha) is 0.0 in the pre-7f path; `forcing` is None in the pre-7f
/// path. When both are inert the accumulated RHS is byte-for-byte the original
/// `nonlinear + nu*lap_omega` expression (backward-compatibility guarantee).
#[allow(clippy::too_many_arguments)]
fn vorticity_rhs(
    omega: &[f64],
    u: &[f64],
    v: &[f64],
    nu: f64,
    kx: &[f64],
    ky: &[f64],
    plans: &Fft2DPlans,
    nx: usize,
    ny: usize,
    drag: f64,
    forcing: Option<&[f64]>,
) -> Vec<f64> {
    let i_unit = Complex64::new(0.0, 1.0);
    let omega_hat = plans.fft2(omega);

    let domega_dx: Vec<f64> = {
        let spec: Vec<Complex64> = (0..nx * ny)
            .map(|idx| i_unit * kx[idx] * omega_hat[idx])
            .collect();
        plans.ifft2_real(&spec)
    };

    let domega_dy: Vec<f64> = {
        let spec: Vec<Complex64> = (0..nx * ny)
            .map(|idx| i_unit * ky[idx] * omega_hat[idx])
            .collect();
        plans.ifft2_real(&spec)
    };

    let lap_omega: Vec<f64> = {
        let spec: Vec<Complex64> = (0..nx * ny)
            .map(|idx| {
                let k2 = kx[idx].powi(2) + ky[idx].powi(2);
                -k2 * omega_hat[idx]
            })
            .collect();
        plans.ifft2_real(&spec)
    };

    // Nonlinear term with 2/3-rule dealiasing (KH-NM-4, Orszag 1971)
    // Same pattern as Python oracle commit 65acf10 and Scala/C++ backends.
    let nonlinear_phys: Vec<f64> = (0..nx * ny)
        .map(|idx| -u[idx] * domega_dx[idx] - v[idx] * domega_dy[idx])
        .collect();
    let mut nl_hat = plans.fft2(&nonlinear_phys);
    dealias_two_thirds(&mut nl_hat, kx, ky, nx, ny);
    let nonlinear = plans.ifft2_real(&nl_hat);

    (0..nx * ny)
        .map(|idx| {
            let mut r = nonlinear[idx] + nu * lap_omega[idx];
            if drag != 0.0 {
                r -= drag * omega[idx];
            }
            if let Some(f) = forcing {
                r += f[idx];
            }
            r
        })
        .collect()
}

// ── Single RK4 step ───────────────────────────────────────────────────────────

fn add_scaled(a: &[f64], b: &[f64], scale: f64) -> Vec<f64> {
    a.iter().zip(b.iter()).map(|(&x, &y)| x + scale * y).collect()
}

#[allow(clippy::too_many_arguments)]
fn rk4_step(
    omega: &[f64],
    nu: f64,
    kx: &[f64],
    ky: &[f64],
    plans: &Fft2DPlans,
    nx: usize,
    ny: usize,
    dt: f64,
    t: f64,
    drag: f64,
    forcing: Option<&Forcing>,
) -> Vec<f64> {
    // Stage forcing fields evaluated at RK4 stage times (k1@t, k2/k3@t+dt/2,
    // k4@t+dt). None when forcing is off (pre-7f path bit-identical).
    let f_at = |ts: f64| -> Option<Vec<f64>> { forcing.map(|fc| fc.eval(ts)) };
    let f1 = f_at(t);
    let f23 = f_at(t + 0.5 * dt);
    let f4 = f_at(t + dt);

    let step = |w: &[f64], f: &Option<Vec<f64>>| -> Vec<f64> {
        let w_hat = plans.fft2(w);
        let psi_hat = solve_poisson(&w_hat, kx, ky, nx, ny);
        let (u, v) = velocity_from_psi_hat(&psi_hat, kx, ky, plans, nx, ny);
        vorticity_rhs(w, &u, &v, nu, kx, ky, plans, nx, ny, drag, f.as_deref())
    };

    let k1 = step(omega, &f1);
    let k2 = step(&add_scaled(omega, &k1, 0.5 * dt), &f23);
    let k3 = step(&add_scaled(omega, &k2, 0.5 * dt), &f23);
    let k4 = step(&add_scaled(omega, &k3, dt), &f4);

    (0..nx * ny)
        .map(|idx| {
            omega[idx] + (dt / 6.0) * (k1[idx] + 2.0 * k2[idx] + 2.0 * k3[idx] + k4[idx])
        })
        .collect()
}

// ── Diagnostics ───────────────────────────────────────────────────────────────

fn mean(v: &[f64]) -> f64 {
    v.iter().sum::<f64>() / v.len() as f64
}

#[allow(clippy::too_many_arguments)]
fn compute_diagnostics(
    omega: &[f64],
    u: &[f64],
    v: &[f64],
    kx: &[f64],
    ky: &[f64],
    plans: &Fft2DPlans,
    nx: usize,
    ny: usize,
) -> Diagnostics {
    let i_unit = Complex64::new(0.0, 1.0);

    let ke = 0.5 * mean(&u.iter().zip(v.iter()).map(|(&ui, &vi)| ui * ui + vi * vi).collect::<Vec<_>>());
    let enstrophy = 0.5 * mean(&omega.iter().map(|&w| w * w).collect::<Vec<_>>());
    let max_vorticity = omega.iter().map(|&w| w.abs()).fold(0.0_f64, f64::max);

    // Divergence: d(u)/dx + d(v)/dy should be machine-zero for spectral method
    let u_hat = plans.fft2(u);
    let v_hat = plans.fft2(v);
    let div_spec: Vec<Complex64> = (0..nx * ny)
        .map(|idx| i_unit * kx[idx] * u_hat[idx] + i_unit * ky[idx] * v_hat[idx])
        .collect();
    let div = plans.ifft2_real(&div_spec);
    let divergence_rms = mean(&div.iter().map(|&d| d * d).collect::<Vec<_>>()).sqrt();

    Diagnostics { kinetic_energy: ke, enstrophy, max_vorticity, divergence_rms, cfl_peak: 0.0 }
}

// ── Main entry point ──────────────────────────────────────────────────────────

pub fn simulate(req: &SimulationRequest) -> SimulationResult {
    let start = std::time::Instant::now();

    let nx    = req.grid_nx;
    let ny    = req.grid_ny;
    let lx    = req.domain_lx;
    let ly    = req.domain_ly;
    let dt    = req.dt;
    let steps = req.steps;
    // nu = |U| / Re: a reversed shear is the mirrored flow, not negative viscosity
    // (modelling rules v0.1, K3). Identical for U >= 0.
    let nu    = req.velocity_shear.abs() / req.reynolds_number;
    let delta = 0.05 * ly;

    let dx = lx / nx as f64;
    let dy = ly / ny as f64;

    let plans = Fft2DPlans::new(nx, ny);
    let (kx, ky) = make_wavenumbers(nx, ny, dx, dy);

    // Phase 7f: build the forcing operator once (None => pre-7f path) + drag coeff.
    let forcing = build_forcing(req, nx, ny, lx, ly);
    let drag = req.linear_drag;

    // Initialise vorticity
    let mut omega: Vec<f64> = match &req.initial_omega {
        Some(v) => {
            assert_eq!(v.len(), nx * ny, "initial_omega length must equal nx*ny");
            v.clone()
        }
        None => initial_conditions(nx, ny, lx, ly, req.velocity_shear, delta,
                                   req.perturbation_amplitude, req.perturbation_mode),
    };

    // CFL guard (KH-NM-5): warn once if user-supplied dt exceeds stability bound
    {
        let omega_hat0 = plans.fft2(&omega);
        let psi_hat0   = solve_poisson(&omega_hat0, &kx, &ky, nx, ny);
        let (u0, v0)   = velocity_from_psi_hat(&psi_hat0, &kx, &ky, &plans, nx, ny);
        let u_max = u0.iter().map(|&x| x.abs()).fold(0.0_f64, f64::max);
        let v_max = v0.iter().map(|&x| x.abs()).fold(0.0_f64, f64::max);
        let dt_bound = cfl_dt_bound(u_max, v_max, dx, dy, nu);
        if dt > dt_bound {
            eprintln!("KH-NM-5: dt={:.3e} exceeds CFL bound {:.3e}; simulation may be unstable",
                      dt, dt_bound);
        }
    }

    // Time integration — track peak CFL pre-step (matches Python oracle timing)
    // t starts at t_start: a continued run keeps the forcing phase (2026-10-10).
    let mut t = req.t_start;
    let mut cfl_peak = 0.0_f64;
    for _ in 0..steps {
        // CFL from pre-step omega (velocity at start of this step, i.e. at current t)
        let w_hat_pre   = plans.fft2(&omega);
        let psi_pre     = solve_poisson(&w_hat_pre, &kx, &ky, nx, ny);
        let (us, vs)    = velocity_from_psi_hat(&psi_pre, &kx, &ky, &plans, nx, ny);
        let u_max = us.iter().fold(0.0_f64, |a, &b| a.max(b.abs()));
        let v_max = vs.iter().fold(0.0_f64, |a, &b| a.max(b.abs()));
        let cfl_step = dt * (u_max / dx + v_max / dy);
        cfl_peak = cfl_peak.max(cfl_step);

        omega = rk4_step(&omega, nu, &kx, &ky, &plans, nx, ny, dt, t, drag, forcing.as_ref());
        t += dt;
    }

    // Final field recovery (reuse last-step psi_hat by re-computing from final omega)
    let omega_hat = plans.fft2(&omega);
    let psi_hat   = solve_poisson(&omega_hat, &kx, &ky, nx, ny);
    let psi       = plans.ifft2_real(&psi_hat);
    let (u, v)    = velocity_from_psi_hat(&psi_hat, &kx, &ky, &plans, nx, ny);

    let mut diagnostics = compute_diagnostics(&omega, &u, &v, &kx, &ky, &plans, nx, ny);
    diagnostics.cfl_peak = cfl_peak;

    let elapsed_ms = start.elapsed().as_secs_f64() * 1000.0;

    SimulationResult {
        backend:         "rust-axum",
        language:        "Rust",
        steps_completed: steps,
        t_final:         (t * 1e10).round() / 1e10,
        grid_nx:         nx,
        grid_ny:         ny,
        u_velocity:      u,
        v_velocity:      v,
        vorticity:       omega,
        pressure:        psi,
        diagnostics,
        compute_time_ms: (elapsed_ms * 100.0).round() / 100.0,
        warnings:        Vec::new(),
    }
}
