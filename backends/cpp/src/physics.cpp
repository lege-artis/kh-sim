// physics.cpp -- KH instability physics kernel (C++ port of kh_physics.py)
//
// Physics: 2D incompressible Navier-Stokes, vorticity-streamfunction form.
// Numerics: pseudo-spectral (Cooley-Tukey in-place radix-2 FFT), RK4.
//
// Reference: kh-sim/shared/physics/KH-PHYSICS.md
// Canonical:  kh-sim/shared/physics/kh_physics.py
//
// Phase 7f extension: forcing API (v0.1 "mode"/"band" + v0.2 rotating band
//   + linear drag). Design: development record, forcing-API design v0.1.
// Citation (C3): Numerical Recipes 2007 3rd ed., Press et al., sec. 12.2 eq.
//   12.2.2 for the Cooley-Tukey FFT convention (forward: exp(-2*pi*i*k*n/N)).

#include "physics.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <complex>
#include <functional>
#include <limits>
#include <iostream>
#include <numbers>    // std::numbers::pi (C++20)
#include <numeric>
#include <stdexcept>
#include <vector>

namespace khsim {

using cd = std::complex<double>;
using vc = std::vector<cd>;
using vd = std::vector<double>;

static constexpr double PI = std::numbers::pi;

// ---- Cooley-Tukey in-place radix-2 FFT (requires n = power of 2) -----------
//
// C3 citation: Press et al. (2007) Numerical Recipes 3rd ed., sec. 12.2,
//   eq. 12.2.2. Forward sign convention: exp(-2*pi*i*k*n/N) (standard DFT).
//   Phase-6c sign fix (commit 118c345): forward uses ang=-2pi/len (minus sign).
// Complexity: O(N log N) -- log2(N) stages of N/2 butterflies.

static void fft1d(vc& a, bool inverse) {
    const int n = static_cast<int>(a.size());
    // A length that is not a power of two makes the butterflies below index past
    // the end of a (heap corruption, then a crash). The HTTP layer refuses such
    // grids with 400; this guard is the second line of defence (finding v0.5, C-1).
    if (n < 1 || (n & (n - 1)) != 0)
        throw std::invalid_argument("fft1d: length " + std::to_string(n) + " is not a power of two");
    // bit-reversal permutation
    for (int i = 1, j = 0; i < n; ++i) {
        int bit = n >> 1;
        for (; j & bit; bit >>= 1) j ^= bit;
        j ^= bit;
        if (i < j) std::swap(a[i], a[j]);
    }
    // butterfly stages
    for (int len = 2; len <= n; len <<= 1) {
        double ang = 2.0 * PI / len * (inverse ? 1.0 : -1.0);  // standard: forward -sign
        cd wlen(std::cos(ang), std::sin(ang));
        for (int i = 0; i < n; i += len) {
            cd w(1.0, 0.0);
            for (int j = 0; j < len / 2; ++j) {
                cd u = a[i + j];
                cd v = a[i + j + len / 2] * w;
                a[i + j]           = u + v;
                a[i + j + len / 2] = u - v;
                w *= wlen;
            }
        }
    }
    if (inverse)
        for (auto& x : a) x /= static_cast<double>(n);
}

// ---- 2D FFT: row-then-column (row-major layout, nx rows of ny elements) -----

static vc fft2(const vd& real, int nx, int ny) {
    vc buf(nx * ny);
    for (int k = 0; k < nx * ny; ++k)
        buf[k] = {real[k], 0.0};

    // FFT along axis-1 (rows)
    vc row(ny);
    for (int i = 0; i < nx; ++i) {
        for (int j = 0; j < ny; ++j) row[j] = buf[i * ny + j];
        fft1d(row, false);
        for (int j = 0; j < ny; ++j) buf[i * ny + j] = row[j];
    }

    // FFT along axis-0 (columns)
    vc col(nx);
    for (int j = 0; j < ny; ++j) {
        for (int i = 0; i < nx; ++i) col[i] = buf[i * ny + j];
        fft1d(col, false);
        for (int i = 0; i < nx; ++i) buf[i * ny + j] = col[i];
    }
    return buf;
}

static vd ifft2_real(vc buf, int nx, int ny) {
    // IFFT along axis-0 (columns)
    vc col(nx);
    for (int j = 0; j < ny; ++j) {
        for (int i = 0; i < nx; ++i) col[i] = buf[i * ny + j];
        fft1d(col, true);
        for (int i = 0; i < nx; ++i) buf[i * ny + j] = col[i];
    }
    // IFFT along axis-1 (rows)
    vc row(ny);
    for (int i = 0; i < nx; ++i) {
        for (int j = 0; j < ny; ++j) row[j] = buf[i * ny + j];
        fft1d(row, true);
        for (int j = 0; j < ny; ++j) buf[i * ny + j] = row[j];
    }
    vd out(nx * ny);
    for (int k = 0; k < nx * ny; ++k) out[k] = buf[k].real();
    return out;
}

// ---- Wavenumbers ------------------------------------------------------------

static vd angular_fftfreq(int n, double d) {
    vd out(n);
    int half = n / 2;
    for (int i = 0; i <= half; ++i)
        out[i] = 2.0 * PI * i / (n * d);
    for (int i = half + 1; i < n; ++i)
        out[i] = 2.0 * PI * (i - n) / (n * d);
    return out;
}

// ---- Initial conditions -----------------------------------------------------

static vd initial_conditions(int nx, int ny,
                              double lx, double ly,
                              double u0, double delta,
                              double amp, int mode) {
    double dx = lx / nx, dy = ly / ny;
    vd omega(nx * ny);
    for (int i = 0; i < nx; ++i) {
        double x = i * dx;
        double pert_x = amp * (2.0 * PI * mode / lx)
                            * std::cos(2.0 * PI * mode * x / lx);
        for (int j = 0; j < ny; ++j) {
            double y = j * dy;
            double z = (y - ly / 2.0) / delta;
            double ch = std::cosh(z);
            omega[i * ny + j] = u0 / delta / (ch * ch) + pert_x;
        }
    }
    return omega;
}

// ---- Poisson solve ----------------------------------------------------------

static vc solve_poisson(const vc& omega_hat,
                        const vd& kx, const vd& ky,
                        int nx, int ny) {
    vc psi_hat(omega_hat);
    for (int i = 0; i < nx; ++i) {
        for (int j = 0; j < ny; ++j) {
            int idx = i * ny + j;
            if (i == 0 && j == 0) {
                psi_hat[0] = {0.0, 0.0};
            } else {
                double k2 = kx[idx] * kx[idx] + ky[idx] * ky[idx];
                psi_hat[idx] = omega_hat[idx] / k2;
            }
        }
    }
    return psi_hat;
}

// ---- Velocity from psi ------------------------------------------------------

static std::pair<vd, vd> velocity_from_psi(const vc& psi_hat,
                                            const vd& kx, const vd& ky,
                                            int nx, int ny) {
    int n = nx * ny;
    vc u_hat(n), v_hat(n);
    cd I{0.0, 1.0};
    for (int k = 0; k < n; ++k) {
        u_hat[k] =  I * ky[k] * psi_hat[k];
        v_hat[k] = -I * kx[k] * psi_hat[k];
    }
    return {ifft2_real(std::move(u_hat), nx, ny),
            ifft2_real(std::move(v_hat), nx, ny)};
}

// ---- 2/3-rule dealiasing + CFL bound ----------------------------------------

/// KH-NM-4: 2/3-rule dealiasing (Orszag 1971). Zeros modes outside 2N/3 cutoff in-place.
static void dealias_two_thirds(vc& field_hat, const vd& kx, const vd& ky, int nx, int ny) {
    double kx_max = 0.0, ky_max = 0.0;
    for (int k = 0; k < nx * ny; ++k) {
        kx_max = std::max(kx_max, std::abs(kx[k]));
        ky_max = std::max(ky_max, std::abs(ky[k]));
    }
    double kx_cut = (2.0 / 3.0) * kx_max;
    double ky_cut = (2.0 / 3.0) * ky_max;
    for (int k = 0; k < nx * ny; ++k)
        if (std::abs(kx[k]) > kx_cut || std::abs(ky[k]) > ky_cut)
            field_hat[k] = {0.0, 0.0};
}

/// KH-NM-5: CFL-limited time step. Returns min of advective and diffusive bounds.
static double cfl_dt_bound(double u_max, double v_max, double dx, double dy, double nu) {
    double dt_adv  = std::min(dx / (u_max + 1e-30), dy / (v_max + 1e-30));
    // nu = 0 (inviscid) is a valid limit: no diffusive bound (modelling rules v0.1, D).
    double dt_diff = (nu == 0.0) ? std::numeric_limits<double>::infinity()
                                 : std::min(dx * dx / (2.0 * nu), dy * dy / (2.0 * nu));
    return std::min(dt_adv, dt_diff);
}

// ---- Forcing field construction (PHASE-7F-FORCING-API-DESIGN-v0.1.md) -------
//
// ForcingOperator: wraps the per-stage forcing evaluation.
// For "none" or amplitude==0: nullptr (no overhead in the time loop).
// For "mode": static (time-independent) vd field.
// For "band" v0.1 (omega0==0): static vd field.
// For "band" v0.2 (omega0>0): callable F(t) = precomputed spatial bases.
//
// All formulas are bit-reproducible across backends (no RNG; deterministic
// pseudo-random phases from irrational-number trick per design doc §3/§7).

struct ForcingContext {
    // Precomputed per-mode spatial cosine and sine fields for v0.2 rotating band.
    // F(t) = (A/sqrt(M)) * sum_m [ C_m * cos(Omega_m * t) - S_m * sin(Omega_m * t) ]
    std::vector<vd> C_modes;    // cosine spatial fields, one per mode
    std::vector<vd> S_modes;    // sine  spatial fields, one per mode
    std::vector<double> Omega;  // per-mode rotation rate (rad/time)
    bool rotating = false;
    // For static case (no rotation): just a single precomputed field.
    vd static_field;
    bool is_static = false;
    double drag = 0.0;
};

// Evaluate forcing field at time t.
static vd eval_forcing(const ForcingContext& fc, double t, int n) {
    if (fc.is_static) return fc.static_field;
    // rotating band: sum over modes
    vd F(n, 0.0);
    for (size_t m = 0; m < fc.C_modes.size(); ++m) {
        double cos_t = std::cos(fc.Omega[m] * t);
        double sin_t = std::sin(fc.Omega[m] * t);
        const vd& C = fc.C_modes[m];
        const vd& S = fc.S_modes[m];
        for (int k = 0; k < n; ++k)
            F[k] += C[k] * cos_t - S[k] * sin_t;
    }
    return F;
}

// Fractional part helper (positive definite).
static double frac_pos(double x) {
    double f = x - std::floor(x);
    return f < 0.0 ? f + 1.0 : f;
}

// Build the ForcingContext for the given request.
// Returns a context with is_static=false AND empty C_modes if no forcing.
static ForcingContext build_forcing(const SimulationRequest& req,
                                    const vd& kx_1d, const vd& ky_1d) {
    ForcingContext fc;
    fc.drag = req.linear_drag;

    if (req.forcing_type == "none" || req.forcing_amplitude == 0.0) {
        fc.is_static = true;
        fc.static_field.assign(req.grid_nx * req.grid_ny, 0.0);
        return fc;
    }

    int nx = req.grid_nx, ny = req.grid_ny;
    double lx = req.domain_lx, ly = req.domain_ly;
    double dx = lx / nx, dy = ly / ny;
    double A = req.forcing_amplitude;

    if (req.forcing_type == "mode") {
        // CP-7F-D: steady Taylor-Green-mode forcing.
        // F(x,y) = A * cos(2*pi*kx*x/lx) * cos(2*pi*ky*y/ly)
        // Reference: PHASE-7F-FORCING-API-DESIGN-v0.1.md sec.3
        fc.is_static = true;
        fc.static_field.resize(nx * ny);
        for (int i = 0; i < nx; ++i) {
            double cx = std::cos(2.0 * PI * req.forcing_kx * i * dx / lx);
            for (int j = 0; j < ny; ++j) {
                double cy = std::cos(2.0 * PI * req.forcing_ky * j * dy / ly);
                fc.static_field[i * ny + j] = A * cx * cy;
            }
        }
        return fc;
    }

    if (req.forcing_type == "band") {
        // CP-7F-B: isotropic shell forcing.
        // Constraint: lx must equal ly for isotropic integer-mode radius.
        if (std::abs(lx - ly) > 1e-10 * lx)
            throw std::invalid_argument("band forcing requires lx == ly");

        double kmin = req.forcing_kmin;
        double kmax = req.forcing_kmax;

        // Enumerate contributing modes (half-plane, skip DC-row negative-k):
        // ix >= 0; skip ix==0 AND iy<=0.
        // phi0(ix,iy) = 2*pi * frac(ix*sqrt(2) + iy*sqrt(3))
        // Omega(ix,iy) = 2*pi * frac(ix*sqrt(11) + iy*sqrt(13)) * omega0 (v0.2)
        std::vector<std::pair<int,int>> modes;
        int kmax_int = static_cast<int>(std::ceil(kmax)) + 1;
        for (int ix = 0; ix <= kmax_int; ++ix) {
            int iy_start = (ix == 0) ? 1 : -kmax_int;
            for (int iy = iy_start; iy <= kmax_int; ++iy) {
                double kr = std::sqrt((double)(ix*ix + iy*iy));
                if (kr >= kmin && kr <= kmax)
                    modes.push_back({ix, iy});
            }
        }

        int M = static_cast<int>(modes.size());
        if (M == 0)
            throw std::invalid_argument("band forcing: no modes in [kmin, kmax]");

        double scale = A / std::sqrt((double)M);

        bool rotating = (req.forcing_omega0 != 0.0);

        if (!rotating) {
            // v0.1 static band
            fc.is_static = true;
            fc.static_field.resize(nx * ny, 0.0);
            for (auto [ix, iy] : modes) {
                double phi0 = 2.0 * PI * frac_pos(ix * std::sqrt(2.0) + iy * std::sqrt(3.0));
                for (int i = 0; i < nx; ++i) {
                    double arg_x = 2.0 * PI * ix * i * dx / lx;
                    for (int j = 0; j < ny; ++j) {
                        double arg_y = 2.0 * PI * iy * j * dy / ly;
                        fc.static_field[i * ny + j] +=
                            scale * std::cos(arg_x + arg_y + phi0);
                    }
                }
            }
            return fc;
        }

        // v0.2 rotating band: precompute C_m and S_m spatial fields.
        // F(x,y,t) = sum_m [ C_m(x,y) * cos(Omega_m*t) - S_m(x,y) * sin(Omega_m*t) ]
        // where C_m(x,y) = scale * cos(k_m . x + phi0_m)
        //       S_m(x,y) = scale * sin(k_m . x + phi0_m)
        fc.rotating = true;
        fc.is_static = false;
        fc.C_modes.resize(M, vd(nx * ny));
        fc.S_modes.resize(M, vd(nx * ny));
        fc.Omega.resize(M);

        for (int m = 0; m < M; ++m) {
            auto [ix, iy] = modes[m];
            double phi0 = 2.0 * PI * frac_pos(ix * std::sqrt(2.0) + iy * std::sqrt(3.0));
            fc.Omega[m]  = 2.0 * PI * frac_pos(ix * std::sqrt(11.0) + iy * std::sqrt(13.0))
                           * req.forcing_omega0;
            for (int i = 0; i < nx; ++i) {
                double arg_x = 2.0 * PI * ix * i * dx / lx;
                for (int j = 0; j < ny; ++j) {
                    double arg_y = 2.0 * PI * iy * j * dy / ly;
                    double phase = arg_x + arg_y + phi0;
                    fc.C_modes[m][i * ny + j] = scale * std::cos(phase);
                    fc.S_modes[m][i * ny + j] = scale * std::sin(phase);
                }
            }
        }
        return fc;
    }

    throw std::invalid_argument("Unknown forcing_type: " + req.forcing_type);
}

// ---- Vorticity RHS ----------------------------------------------------------

static vd vorticity_rhs(const vd& omega, const vd& u, const vd& v,
                         double nu, double drag,
                         const vd& kx, const vd& ky,
                         int nx, int ny,
                         const vd& forcing) {
    int n = nx * ny;
    vc omega_hat = fft2(omega, nx, ny);
    cd I{0.0, 1.0};

    vc dx_spec(n), dy_spec(n), lap_spec(n);
    for (int k = 0; k < n; ++k) {
        dx_spec[k]  = I * kx[k] * omega_hat[k];
        dy_spec[k]  = I * ky[k] * omega_hat[k];
        double k2   = kx[k]*kx[k] + ky[k]*ky[k];
        lap_spec[k] = -k2 * omega_hat[k];
    }
    vd dox = ifft2_real(std::move(dx_spec), nx, ny);
    vd doy = ifft2_real(std::move(dy_spec), nx, ny);
    vd lap = ifft2_real(std::move(lap_spec), nx, ny);

    // Nonlinear term with 2/3-rule dealiasing (KH-NM-4, Orszag 1971)
    vd nonlinear_phys(n);
    for (int k = 0; k < n; ++k)
        nonlinear_phys[k] = -u[k] * dox[k] - v[k] * doy[k];
    vc nl_hat = fft2(nonlinear_phys, nx, ny);
    dealias_two_thirds(nl_hat, kx, ky, nx, ny);
    vd nonlinear = ifft2_real(std::move(nl_hat), nx, ny);

    vd rhs(n);
    for (int k = 0; k < n; ++k)
        rhs[k] = nonlinear[k] + nu * lap[k] - drag * omega[k] + forcing[k];
    return rhs;
}

// ---- RK4 step with forcing --------------------------------------------------

static vd add_scaled(const vd& a, const vd& b, double s) {
    vd out(a.size());
    for (size_t k = 0; k < a.size(); ++k) out[k] = a[k] + s * b[k];
    return out;
}

static vd rk4_step_forced(const vd& omega, double nu, double drag,
                           const vd& kx, const vd& ky,
                           int nx, int ny, double dt, double t,
                           const ForcingContext& fc) {
    int n = nx * ny;
    auto f = [&](const vd& w, double t_stage) -> vd {
        vc wh = fft2(w, nx, ny);
        vc ph = solve_poisson(wh, kx, ky, nx, ny);
        auto [u, v] = velocity_from_psi(ph, kx, ky, nx, ny);
        vd F = eval_forcing(fc, t_stage, n);
        return vorticity_rhs(w, u, v, nu, drag, kx, ky, nx, ny, F);
    };
    vd k1 = f(omega, t);
    vd k2 = f(add_scaled(omega, k1, 0.5 * dt), t + 0.5 * dt);
    vd k3 = f(add_scaled(omega, k2, 0.5 * dt), t + 0.5 * dt);
    vd k4 = f(add_scaled(omega, k3, dt),        t + dt);

    vd out(omega.size());
    for (size_t k = 0; k < omega.size(); ++k)
        out[k] = omega[k] + (dt / 6.0) * (k1[k] + 2*k2[k] + 2*k3[k] + k4[k]);
    return out;
}

// ---- Diagnostics ------------------------------------------------------------

static double mean(const vd& v) {
    return std::accumulate(v.begin(), v.end(), 0.0) / v.size();
}

static Diagnostics compute_diagnostics(const vd& omega,
                                       const vd& u, const vd& v,
                                       const vd& kx, const vd& ky,
                                       int nx, int ny) {
    int n = nx * ny;
    vd ke_arr(n), ens_arr(n);
    double max_w = 0.0;
    for (int k = 0; k < n; ++k) {
        ke_arr[k]  = 0.5 * (u[k]*u[k] + v[k]*v[k]);
        ens_arr[k] = 0.5 * omega[k] * omega[k];
        max_w = std::max(max_w, std::abs(omega[k]));
    }

    cd I{0.0, 1.0};
    vc uh = fft2(u, nx, ny);
    vc vh = fft2(v, nx, ny);
    vc div_spec(n);
    for (int k = 0; k < n; ++k)
        div_spec[k] = I * kx[k] * uh[k] + I * ky[k] * vh[k];
    vd div = ifft2_real(std::move(div_spec), nx, ny);
    vd div2(n);
    for (int k = 0; k < n; ++k) div2[k] = div[k] * div[k];

    return {mean(ke_arr), mean(ens_arr), max_w, std::sqrt(mean(div2))};
}

// ---- Main entry point -------------------------------------------------------

SimulationResult simulate(const SimulationRequest& req) {
    using clock = std::chrono::high_resolution_clock;
    auto t0 = clock::now();

    int nx = req.grid_nx, ny = req.grid_ny;
    double lx = req.domain_lx, ly = req.domain_ly;
    double dt = req.dt;
    int nsteps = req.steps;
    // nu = |U| / Re: a reversed shear is the mirrored flow, not negative viscosity
    // (modelling rules v0.1, K3). Identical for U >= 0.
    double nu  = std::abs(req.velocity_shear) / req.reynolds_number;
    double delta = 0.05 * ly;
    double dx = lx / nx, dy = ly / ny;

    vd kx_1d = angular_fftfreq(nx, dx);
    vd ky_1d = angular_fftfreq(ny, dy);
    vd kx(nx * ny), ky(nx * ny);
    for (int i = 0; i < nx; ++i)
        for (int j = 0; j < ny; ++j) {
            kx[i*ny+j] = kx_1d[i];
            ky[i*ny+j] = ky_1d[j];
        }

    vd omega;
    if (req.initial_omega) {
        omega = *req.initial_omega;
        if ((int)omega.size() != nx * ny)
            throw std::invalid_argument("initial_omega size mismatch");
    } else {
        omega = initial_conditions(nx, ny, lx, ly,
                                   req.velocity_shear, delta,
                                   req.perturbation_amplitude,
                                   req.perturbation_mode);
    }

    // Build forcing context (once before the time loop)
    ForcingContext fc = build_forcing(req, kx_1d, ky_1d);

    // CFL guard (KH-NM-5): warn once if user-supplied dt exceeds stability bound
    {
        vc oh0 = fft2(omega, nx, ny);
        vc ph0 = solve_poisson(oh0, kx, ky, nx, ny);
        auto [u0, v0] = velocity_from_psi(ph0, kx, ky, nx, ny);
        double u_max = *std::max_element(u0.begin(), u0.end(),
            [](double a, double b) { return std::abs(a) < std::abs(b); });
        double v_max = *std::max_element(v0.begin(), v0.end(),
            [](double a, double b) { return std::abs(a) < std::abs(b); });
        double dt_bound = cfl_dt_bound(std::abs(u_max), std::abs(v_max), dx, dy, nu);
        if (dt > dt_bound)
            std::cerr << "KH-NM-5: dt=" << dt << " exceeds CFL bound " << dt_bound
                      << "; simulation may be unstable\n";
    }

    // Time loop. t starts at t_start: a continued run keeps the phase of
    // time-dependent forcing (CP-7F-B/D harness, 2026-10-10).
    double cfl_peak = 0.0;
    double t = req.t_start;
    int steps_done = 0;
    for (int s = 0; s < nsteps; ++s) {
        // CFL diagnostic
        vc wh_step = fft2(omega, nx, ny);
        vc ph_step = solve_poisson(wh_step, kx, ky, nx, ny);
        auto [u_step, v_step] = velocity_from_psi(ph_step, kx, ky, nx, ny);
        double u_max = 0.0, v_max = 0.0;
        for (int k = 0; k < nx * ny; ++k) {
            u_max = std::max(u_max, std::abs(u_step[k]));
            v_max = std::max(v_max, std::abs(v_step[k]));
        }
        double cfl_step = dt * (u_max / dx + v_max / dy);
        cfl_peak = std::max(cfl_peak, cfl_step);

        omega = rk4_step_forced(omega, nu, fc.drag, kx, ky, nx, ny, dt, t, fc);
        t += dt;
        ++steps_done;
    }

    vc oh  = fft2(omega, nx, ny);
    vc ph  = solve_poisson(oh, kx, ky, nx, ny);
    vd psi = ifft2_real(ph, nx, ny);
    auto [u, v] = velocity_from_psi(ph, kx, ky, nx, ny);
    Diagnostics diag = compute_diagnostics(omega, u, v, kx, ky, nx, ny);
    diag.cfl_peak = cfl_peak;

    double elapsed_ms = std::chrono::duration<double, std::milli>(
        clock::now() - t0).count();

    // round t to avoid floating-point drift in JSON
    t = std::round(t * 1e10) / 1e10;

    return {
        "cpp-httplib", "C++",
        steps_done, t, nx, ny,
        std::move(u), std::move(v), std::move(omega), std::move(psi),
        diag,
        std::round(elapsed_ms * 100.0) / 100.0,
    };
}

} // namespace khsim
