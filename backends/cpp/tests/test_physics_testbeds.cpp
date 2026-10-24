// test_physics_testbeds.cpp — Phase 7d physics testbeds, C++ Mode 1 HTTP.
//
// Ports tests/integration/test_physics_testbeds_python.py to C++.
// Mode 1: POST /simulate to http://localhost:8003
// backend_factor = 1.0 (post-FFT-sign-fix; matches Python ref at machine precision)
// Reference: development record, phase-7d brief

#include <algorithm>
#include <cmath>
#include <functional>
#include <iomanip>
#include <iostream>
#include <numbers>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include <httplib.h>
#include <nlohmann/json.hpp>

using json = nlohmann::json;

static constexpr double PI = std::numbers::pi;
static const std::string HOST = "localhost";
static constexpr int PORT = 8003;

// ── HTTP helper ───────────────────────────────────────────────────────────────

static json post_simulate(const json& payload, int timeout_s = 120) {
    httplib::Client cli(HOST, PORT);
    cli.set_connection_timeout(30, 0);
    cli.set_read_timeout(timeout_s, 0);
    auto res = cli.Post("/simulate", payload.dump(), "application/json");
    if (!res)
        throw std::runtime_error("HTTP request failed (no response)");
    if (res->status != 200)
        throw std::runtime_error("HTTP " + std::to_string(res->status) + ": " + res->body);
    return json::parse(res->body);
}

// ── Payload builder ───────────────────────────────────────────────────────────

static json make_payload(int nx, int ny, double lx, double ly,
                          double dt, int steps,
                          double re, double shear,
                          double amp, int mode,
                          const std::vector<double>* omega0 = nullptr) {
    json j;
    j["grid_nx"]                 = nx;
    j["grid_ny"]                 = ny;
    j["domain_lx"]               = lx;
    j["domain_ly"]               = ly;
    j["dt"]                      = dt;
    j["steps"]                   = steps;
    j["reynolds_number"]         = re;
    j["velocity_shear"]          = shear;
    j["perturbation_amplitude"]  = amp;
    j["perturbation_mode"]       = mode;
    if (omega0)
        j["initial_omega"] = *omega0;
    return j;
}

// ── IC helper — mirrors physics.cpp initial_conditions() exactly ──────────────

static std::vector<double> make_ic(int nx, int ny, double lx, double ly,
                                    double u0, double delta,
                                    double amp, int mode) {
    double dx = lx / nx, dy = ly / ny;
    std::vector<double> omega(nx * ny);
    for (int i = 0; i < nx; ++i) {
        double x     = i * dx;
        double pert_x = amp * (2.0 * PI * mode / lx)
                             * std::cos(2.0 * PI * mode * x / lx);
        for (int j = 0; j < ny; ++j) {
            double y  = j * dy;
            double z  = (y - ly / 2.0) / delta;
            double ch = std::cosh(z);
            omega[i * ny + j] = u0 / delta / (ch * ch) + pert_x;
        }
    }
    return omega;
}

// ── Stats helpers ─────────────────────────────────────────────────────────────

static double rms_diff(const std::vector<double>& a, const std::vector<double>& b) {
    double s = 0.0;
    for (size_t k = 0; k < a.size(); ++k) {
        double d = a[k] - b[k];
        s += d * d;
    }
    return std::sqrt(s / static_cast<double>(a.size()));
}

static double rms_of(const std::vector<double>& a) {
    double s = 0.0;
    for (double x : a) s += x * x;
    return std::sqrt(s / static_cast<double>(a.size()));
}

static double sum_of(const std::vector<double>& a) {
    return std::accumulate(a.begin(), a.end(), 0.0);
}

static double mean_of(const std::vector<double>& a) {
    return sum_of(a) / static_cast<double>(a.size());
}

// ── Log-log slope via ordinary least squares ─────────────────────────────────

static double log_slope(const std::vector<double>& xs, const std::vector<double>& ys) {
    int n = static_cast<int>(xs.size());
    double sx = 0.0, sy = 0.0, sxx = 0.0, sxy = 0.0;
    for (int i = 0; i < n; ++i) {
        double lx = std::log(xs[i]);
        double ly = std::log(ys[i]);
        sx  += lx;
        sy  += ly;
        sxx += lx * lx;
        sxy += lx * ly;
    }
    return (n * sxy - sx * sy) / (n * sxx - sx * sx);
}

// ── Test result ───────────────────────────────────────────────────────────────

struct TR {
    std::string id;
    std::string status;  // "PASS", "FAIL", "SKIP"
    std::string msg;
};

// ── TB-1a: Taylor-Green vortex ────────────────────────────────────────────────
//
// Exact solution: omega(x,y,t) = -2*cos(x)*cos(y)*exp(-2*nu*t).
// IC: omega_0 = -2*cos(x)*cos(y) on [0, 2pi]^2. The TG advective term
// vanishes identically, leaving pure diffusion and exact decay.
// Pass: rms(omega_sim - omega_exact) / rms(omega_exact) < 1e-6.

static TR tb1a() {
    const int    nx = 64,  ny = 64;
    const double lx = 2.0 * PI, ly = 2.0 * PI;
    const double nu = 0.01, T = 1.0, dt = 0.001;
    const int    nsteps = static_cast<int>(T / dt);  // 1000

    double dx = lx / nx, dy = ly / ny;
    std::vector<double> omega0(nx * ny), omega_exact(nx * ny);
    double decay = std::exp(-2.0 * nu * T);
    for (int i = 0; i < nx; ++i) {
        for (int j = 0; j < ny; ++j) {
            double v = -2.0 * std::cos(i * dx) * std::cos(j * dy);
            omega0[i * ny + j]     = v;
            omega_exact[i * ny + j] = v * decay;
        }
    }

    auto resp  = post_simulate(
        make_payload(nx, ny, lx, ly, dt, nsteps, 1.0 / nu, 1.0, 0.0, 0, &omega0));
    auto vort  = resp["vorticity"].get<std::vector<double>>();

    double rms_d   = rms_diff(vort, omega_exact);
    double rms_ex  = rms_of(omega_exact);
    double rms_rel = rms_d / rms_ex;

    std::ostringstream os;
    os << std::scientific << std::setprecision(3)
       << "rms_rel=" << rms_rel << " tol=1e-6";
    if (rms_rel < 1e-6) return {"TB-1a", "PASS", os.str()};
    return {"TB-1a", "FAIL", os.str()};
}

// ── TB-1b: pure diffusion ─────────────────────────────────────────────────────
//
// Single-mode IC: omega_0 = sin(kx*x)*sin(ky*y).
// Expected decay: amp(T)/amp(0) = exp(-nu*(kx^2+ky^2)*T).
// Mode amplitude extracted via direct DFT sum (NumPy-compatible convention).
// Pass: |ratio - expected_decay| / expected_decay < 1e-4.

static TR tb1b() {
    const int    nx = 32, ny = 32;
    const double lx = 2.0 * PI, ly = 2.0 * PI;
    const int    kxt = 1, kyt = 1;
    const double nu = 0.1, dt = 0.001;
    const int    nsteps = 500;
    const double T  = dt * nsteps;

    double dx = lx / nx, dy = ly / ny;
    std::vector<double> omega0(nx * ny);
    for (int i = 0; i < nx; ++i)
        for (int j = 0; j < ny; ++j)
            omega0[i * ny + j] = std::sin(kxt * i * dx) * std::sin(kyt * j * dy);

    double k2             = kxt * kxt + kyt * kyt;
    double expected_decay = std::exp(-nu * k2 * T);

    auto resp = post_simulate(
        make_payload(nx, ny, lx, ly, dt, nsteps, 1.0 / nu, 1.0, 0.0, 0, &omega0));
    auto vort = resp["vorticity"].get<std::vector<double>>();

    // DFT at (kxt, kyt): X[k] = sum_ij x[i,j] * exp(-2pi*i*(kxt*i/nx + kyt*j/ny))
    auto dft_amp = [&](const std::vector<double>& field) {
        double re = 0.0, im = 0.0;
        for (int i = 0; i < nx; ++i) {
            for (int j = 0; j < ny; ++j) {
                double ang = -2.0 * PI * (kxt * i / static_cast<double>(nx)
                                         + kyt * j / static_cast<double>(ny));
                re += field[i * ny + j] * std::cos(ang);
                im += field[i * ny + j] * std::sin(ang);
            }
        }
        return std::sqrt(re * re + im * im);
    };

    double amp0 = dft_amp(omega0);
    double ampt = dft_amp(vort);
    if (amp0 < 1e-10)
        return {"TB-1b", "SKIP", "initial mode amplitude too small"};

    double ratio   = ampt / amp0;
    double rel_err = std::abs(ratio - expected_decay) / expected_decay;

    std::ostringstream os;
    os << std::fixed << std::setprecision(8)
       << "ratio=" << ratio << " exp=" << expected_decay
       << std::scientific << std::setprecision(3)
       << " rel_err=" << rel_err << " tol=1e-4";
    if (rel_err < 1e-4) return {"TB-1b", "PASS", os.str()};
    return {"TB-1b", "FAIL", os.str()};
}

// ── TB-1c: linear growth rate (two-run subtraction) ───────────────────────────
//
// Small-amplitude perturbation; measure energy growth rate of the difference
// field (perturbed - base) vs analytic sigma_expected = 0.18*U0/delta = 7.2/s.
// Pass: |sigma_sim - sigma_expected| / sigma_expected < 10%.

static TR tb1c() {
    const int    nx = 128, ny = 64;
    const double lx = 1.0, ly = 0.5;
    const double dt = 0.001;
    const int    nsteps = 50;
    const double amp_pert = 1e-4;
    const int    mode  = 2;
    const double re    = 1000.0, shear = 1.0;
    const double delta = 0.05 * ly;                      // 0.025
    const double sigma_expected = 0.18 * shear / delta;  // 7.2 /s

    // Client-side IC — matches server physics.cpp initial_conditions exactly
    auto omega_base_0 = make_ic(nx, ny, lx, ly, shear, delta, 0.0,      0);
    auto omega_pert_0 = make_ic(nx, ny, lx, ly, shear, delta, amp_pert, mode);

    // Base run (no perturbation)
    auto resp_base = post_simulate(
        make_payload(nx, ny, lx, ly, dt, nsteps, re, shear, 0.0, 0));
    auto omega_base_T = resp_base["vorticity"].get<std::vector<double>>();

    // Perturbed run
    auto resp_pert = post_simulate(
        make_payload(nx, ny, lx, ly, dt, nsteps, re, shear, amp_pert, mode));
    auto omega_pert_T = resp_pert["vorticity"].get<std::vector<double>>();

    // Perturbation energy at t=0
    double E_pert_0 = 0.0;
    for (size_t k = 0; k < omega_base_0.size(); ++k) {
        double d = omega_pert_0[k] - omega_base_0[k];
        E_pert_0 += 0.5 * d * d;
    }
    E_pert_0 /= static_cast<double>(omega_base_0.size());

    // Perturbation energy at T_meas
    double E_pert_T = 0.0;
    for (size_t k = 0; k < omega_base_T.size(); ++k) {
        double d = omega_pert_T[k] - omega_base_T[k];
        E_pert_T += 0.5 * d * d;
    }
    E_pert_T /= static_cast<double>(omega_base_T.size());

    if (E_pert_0 < 1e-20 || E_pert_T < 1e-20) {
        std::ostringstream os;
        os << std::scientific << std::setprecision(3)
           << "E0=" << E_pert_0 << " ET=" << E_pert_T;
        return {"TB-1c", "SKIP", "perturbation energy too small; " + os.str()};
    }

    double T_meas   = dt * nsteps;
    double sigma_sim = std::log(E_pert_T / E_pert_0) / (2.0 * T_meas);
    double rel_err   = std::abs(sigma_sim - sigma_expected) / sigma_expected;

    std::ostringstream os;
    os << std::fixed << std::setprecision(4)
       << "sigma_sim=" << sigma_sim << " sigma_exp=" << sigma_expected
       << std::scientific << std::setprecision(3)
       << " rel_err=" << rel_err << " tol=10%";
    if (rel_err < 0.10) return {"TB-1c", "PASS", os.str()};
    return {"TB-1c", "FAIL", os.str()};
}

// ── TB-2a: total vorticity conservation ──────────────────────────────────────
//
// Canonical KH IC has sum(omega) ~ 8192 (positive shear layer).
// Spectral DC mode is exactly conserved: |sum_T - sum_0| < 1e-8.

static TR tb2a() {
    const int    nx = 64, ny = 32;
    const double lx = 1.0, ly = 0.5;
    const double delta = 0.05 * ly;

    auto omega0 = make_ic(nx, ny, lx, ly, 1.0, delta, 0.01, 2);
    double sum_0 = sum_of(omega0);

    auto resp = post_simulate(
        make_payload(nx, ny, lx, ly, 0.001, 100, 1000.0, 1.0, 0.01, 2));
    auto vort  = resp["vorticity"].get<std::vector<double>>();
    double sum_T = sum_of(vort);

    double drift = std::abs(sum_T - sum_0);

    std::ostringstream os;
    os << std::scientific << std::setprecision(3)
       << "drift=" << drift << " tol=1e-8 (sum_0=" << sum_0 << ")";
    if (drift < 1e-8) return {"TB-2a", "PASS", os.str()};
    return {"TB-2a", "FAIL", os.str()};
}

// ── TB-2b: inviscid energy conservation ──────────────────────────────────────
//
// Near-inviscid (nu=1e-8): |E(T)-E(0)| / E(0) < 1e-4.
// E(0) via steps=0 call (server returns diagnostics on IC, no time-stepping).
// Tolerance accounts for 2/3-rule dealiasing energy leak (~7e-7 measured Python).

static TR tb2b() {
    const int    nx = 64, ny = 32;
    const double lx = 1.0, ly = 0.5;
    const double dt = 0.0005, nu = 1e-8, shear = 1.0;
    const double re = shear / nu;

    auto resp0 = post_simulate(
        make_payload(nx, ny, lx, ly, dt, 0, re, shear, 0.01, 2));
    double E0 = resp0["diagnostics"]["kinetic_energy"].get<double>();
    if (E0 < 1e-15)
        return {"TB-2b", "SKIP", "E0 too small"};

    auto respT = post_simulate(
        make_payload(nx, ny, lx, ly, dt, 200, re, shear, 0.01, 2));
    double ET = respT["diagnostics"]["kinetic_energy"].get<double>();

    double rel_err = std::abs(ET - E0) / E0;

    std::ostringstream os;
    os << std::scientific << std::setprecision(3)
       << "E0=" << E0 << " ET=" << ET << " rel_err=" << rel_err << " tol=1e-4";
    if (rel_err < 1e-4) return {"TB-2b", "PASS", os.str()};
    return {"TB-2b", "FAIL", os.str()};
}

// ── TB-2c: inviscid enstrophy conservation ────────────────────────────────────
//
// Near-inviscid (nu=1e-8): |Z(T)-Z(0)| / Z(0) < 1e-4.
// Same rationale as TB-2b.

static TR tb2c() {
    const int    nx = 64, ny = 32;
    const double lx = 1.0, ly = 0.5;
    const double dt = 0.0005, nu = 1e-8, shear = 1.0;
    const double re = shear / nu;

    auto resp0 = post_simulate(
        make_payload(nx, ny, lx, ly, dt, 0, re, shear, 0.01, 2));
    double Z0 = resp0["diagnostics"]["enstrophy"].get<double>();
    if (Z0 < 1e-15)
        return {"TB-2c", "SKIP", "Z0 too small"};

    auto respT = post_simulate(
        make_payload(nx, ny, lx, ly, dt, 200, re, shear, 0.01, 2));
    double ZT = respT["diagnostics"]["enstrophy"].get<double>();

    double rel_err = std::abs(ZT - Z0) / Z0;

    std::ostringstream os;
    os << std::scientific << std::setprecision(3)
       << "Z0=" << Z0 << " ZT=" << ZT << " rel_err=" << rel_err << " tol=1e-4";
    if (rel_err < 1e-4) return {"TB-2c", "PASS", os.str()};
    return {"TB-2c", "FAIL", os.str()};
}

// ── TB-2d: energy decay rate (DC-mode corrected) ─────────────────────────────
//
// dE/dt = -2*nu*Z_ac where Z_ac = Z - Z_dc (AC enstrophy, DC mode excluded).
// Z_dc = 0.5 * mean(omega)^2 (DC mode = spatial mean of omega field).
// Pass: |dE/dt_fd - (-2*nu*Z_ac)| / |2*nu*Z_ac| < 5%.

static TR tb2d() {
    const int    nx = 64, ny = 32;
    const double lx = 1.0, ly = 0.5;
    const double dt = 0.001, re = 1000.0, shear = 1.0;
    const double nu = shear / re;
    const double T1 = 0.1, T2 = 0.2;

    auto resp1 = post_simulate(
        make_payload(nx, ny, lx, ly, dt, 100, re, shear, 0.01, 2));
    double E1 = resp1["diagnostics"]["kinetic_energy"].get<double>();

    auto resp2 = post_simulate(
        make_payload(nx, ny, lx, ly, dt, 200, re, shear, 0.01, 2));
    double E2 = resp2["diagnostics"]["kinetic_energy"].get<double>();

    // Z_ac at T_mid = 0.15 — need vorticity field for DC correction
    auto respM = post_simulate(
        make_payload(nx, ny, lx, ly, dt, 150, re, shear, 0.01, 2));
    double Z_total = respM["diagnostics"]["enstrophy"].get<double>();
    auto   vort_m  = respM["vorticity"].get<std::vector<double>>();

    // Z_dc = 0.5 * mean(omega)^2 = 0.5 * (FFT[0,0]/(nx*ny))^2
    double mean_w = mean_of(vort_m);
    double Z_dc   = 0.5 * mean_w * mean_w;
    double Z_ac   = Z_total - Z_dc;

    double dEdt_fd     = (E2 - E1) / (T2 - T1);
    double dEdt_theory = -2.0 * nu * Z_ac;
    double denom       = std::abs(dEdt_theory);
    if (denom < 1e-15)
        return {"TB-2d", "SKIP", "denominator too small"};

    double rel_err = std::abs(dEdt_fd - dEdt_theory) / denom;

    std::ostringstream os;
    os << std::scientific << std::setprecision(3)
       << "dE/dt_fd=" << dEdt_fd << " theory=" << dEdt_theory
       << " Z_ac=" << Z_ac << " Z_dc=" << Z_dc
       << " rel_err=" << rel_err << " tol=5%";
    if (rel_err < 0.05) return {"TB-2d", "PASS", os.str()};
    return {"TB-2d", "FAIL", os.str()};
}

// ── TB-2e: incompressibility ──────────────────────────────────────────────────
//
// Spectral vorticity-streamfunction enforces div(u)=0 to machine eps.
// Pass: diagnostics.divergence_rms < 1e-12.

static TR tb2e() {
    auto resp = post_simulate(
        make_payload(64, 32, 1.0, 0.5, 0.001, 100, 1000.0, 1.0, 0.01, 2));
    double div = resp["diagnostics"]["divergence_rms"].get<double>();

    std::ostringstream os;
    os << std::scientific << std::setprecision(3)
       << "div_rms=" << div << " tol=1e-12";
    if (div < 1e-12) return {"TB-2e", "PASS", os.str()};
    return {"TB-2e", "FAIL", os.str()};
}

// ── TB-3a: spatial convergence (Taylor-Green) ─────────────────────────────────
//
// 4 grids: 32^2, 64^2, 128^2, 256^2 on [0, 2pi]^2.
// T=0.1, nu=0.01, dt=0.0001, 1000 steps per grid.
// For pseudo-spectral, TG is exactly representable from 32^2 upward:
// all errors land at machine-epsilon floor.
// Pass: (a) max(errors) < 1e-12  (b) monotone non-increasing (10% slack).

static TR tb3a() {
    const double nu = 0.01, T = 0.1, dt = 0.0001;
    const int    nsteps = static_cast<int>(T / dt);  // 1000
    const std::vector<int> grids = {32, 64, 128, 256};
    std::vector<double> errors;

    for (int N : grids) {
        double lx = 2.0 * PI, ly = 2.0 * PI;
        double dx = lx / N, dy = ly / N;
        double decay = std::exp(-2.0 * nu * T);

        std::vector<double> omega0(N * N), omega_exact(N * N);
        for (int i = 0; i < N; ++i) {
            for (int j = 0; j < N; ++j) {
                double v = -2.0 * std::cos(i * dx) * std::cos(j * dy);
                omega0[i * N + j]      = v;
                omega_exact[i * N + j] = v * decay;
            }
        }

        // 300 s timeout: 256^2 x 1000 steps can take ~30s on a laptop
        auto resp = post_simulate(
            make_payload(N, N, lx, ly, dt, nsteps, 1.0 / nu, 1.0, 0.0, 0, &omega0),
            300);
        auto vort = resp["vorticity"].get<std::vector<double>>();

        errors.push_back(rms_diff(vort, omega_exact));
    }

    double max_err = *std::max_element(errors.begin(), errors.end());

    // Monotone check is only meaningful above machine-epsilon floor.
    // When all errors are < 1e-13 (noise-dominated), noise-level ratios can
    // exceed the 10% slack spuriously (see STATUS-python.md §4.6).
    // Primary gate is max_err < 1e-12; monotone is a secondary check.
    bool monotone = true;
    if (max_err > 1e-13) {
        for (size_t i = 0; i + 1 < errors.size(); ++i)
            if (errors[i + 1] > errors[i] * 1.1) { monotone = false; break; }
    }

    std::ostringstream os;
    os << std::scientific << std::setprecision(3)
       << "max_err=" << max_err << " tol=1e-12; monotone=" << (monotone ? "yes" : "no")
       << " [";
    for (size_t i = 0; i < errors.size(); ++i) {
        if (i) os << " ";
        os << errors[i];
    }
    os << "]";
    if (max_err < 1e-12 && monotone) return {"TB-3a", "PASS", os.str()};
    return {"TB-3a", "FAIL", os.str()};
}

// ── TB-3b: RK4 temporal convergence ─────────────────────────────────────────
//
// Reference: dt_ref=5e-4 (200 steps to T=0.1).
// Test: dt in {1e-2, 5e-3, 2e-3, 1e-3}; log-log slope of rms_error vs dt.
// nu=0.001 keeps all test dt well below diffusive CFL (~0.12).
// Pass: fitted slope >= 3.5 (RK4 theoretical: 4).

static TR tb3b() {
    const int    nx = 64, ny = 32;
    const double lx = 1.0, ly = 0.5;
    const double nu = 0.001, shear = 1.0;
    const double re = shear / nu;
    const double T  = 0.1, dt_ref = 5e-4;
    const int    nsteps_ref = static_cast<int>(T / dt_ref);  // 200

    auto resp_ref = post_simulate(
        make_payload(nx, ny, lx, ly, dt_ref, nsteps_ref, re, shear, 0.01, 2));
    auto omega_ref = resp_ref["vorticity"].get<std::vector<double>>();

    std::vector<double> dt_vals = {1e-2, 5e-3, 2e-3, 1e-3};
    std::vector<double> errs;
    for (double dt : dt_vals) {
        int nst = static_cast<int>(std::round(T / dt));
        auto resp = post_simulate(
            make_payload(nx, ny, lx, ly, dt, nst, re, shear, 0.01, 2));
        auto omega = resp["vorticity"].get<std::vector<double>>();
        errs.push_back(rms_diff(omega, omega_ref));
    }

    double slope = log_slope(dt_vals, errs);

    std::ostringstream os;
    os << std::fixed << std::setprecision(3) << "slope=" << slope << " tol>=3.5  [";
    for (size_t i = 0; i < errs.size(); ++i) {
        if (i) os << " ";
        os << std::scientific << std::setprecision(3) << errs[i];
    }
    os << "]";
    if (slope >= 3.5) return {"TB-3b", "PASS", os.str()};
    return {"TB-3b", "FAIL", os.str()};
}

// ── Pending / deferred stubs ──────────────────────────────────────────────────

static TR tb3c_pending() {
    return {"TB-3c", "SKIP",
            "PENDING-OPUS-7c: awaits MMS design (manufactured solution + source term)"};
}

static TR tb4a_pending() {
    return {"TB-4a", "SKIP",
            "PENDING-OPUS-7e: awaits Lyapunov estimator design"};
}

static TR tb4b_stretch() {
    return {"TB-4b", "SKIP",
            "STRETCH-DEFER per MASTER-BRIEF §2.5 (long-time ensemble; high complexity)"};
}

// ── Main ──────────────────────────────────────────────────────────────────────

int main() {
    std::cout << "=== kh-cpp phase 7d physics testbeds (port " << PORT << ") ===\n\n";

    std::vector<TR> results;

    auto run_test = [&](const std::string& id, std::function<TR()> fn) {
        std::cout << "Running " << id << "... " << std::flush;
        TR r;
        try {
            r = fn();
        } catch (const std::exception& e) {
            r = {id, "FAIL", std::string("exception: ") + e.what()};
        }
        std::cout << r.status;
        if (!r.msg.empty()) std::cout << "  [" << r.msg << "]";
        std::cout << "\n";
        results.push_back(r);
    };

    run_test("TB-1a", tb1a);
    run_test("TB-1b", tb1b);
    run_test("TB-1c", tb1c);
    run_test("TB-2a", tb2a);
    run_test("TB-2b", tb2b);
    run_test("TB-2c", tb2c);
    run_test("TB-2d", tb2d);
    run_test("TB-2e", tb2e);
    run_test("TB-3a", tb3a);
    run_test("TB-3b", tb3b);
    run_test("TB-3c", tb3c_pending);
    run_test("TB-4a", tb4a_pending);
    run_test("TB-4b", tb4b_stretch);

    std::cout << "\n--- Summary ---\n";
    int pass_cnt = 0, fail_cnt = 0, skip_cnt = 0;
    for (auto& r : results) {
        if      (r.status == "PASS") ++pass_cnt;
        else if (r.status == "FAIL") ++fail_cnt;
        else                         ++skip_cnt;
    }
    std::cout << "PASS=" << pass_cnt
              << "  FAIL=" << fail_cnt
              << "  SKIP=" << skip_cnt
              << "  total=" << results.size() << "\n";

    if (fail_cnt > 0) {
        std::cout << "\nFailed tests:\n";
        for (auto& r : results)
            if (r.status == "FAIL")
                std::cout << "  " << r.id << ": " << r.msg << "\n";
    }

    return (fail_cnt == 0) ? 0 : 1;
}
