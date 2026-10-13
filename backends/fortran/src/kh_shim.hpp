#pragma once
// kh_shim.hpp — C++ interface to the Fortran physics kernel
// Declares the C-linkage symbol exported by kh_physics.f90 bind(c)
// and a high-level simulate() function that wraps it.

#include <vector>
#include <optional>
#include <string>

extern "C" {
    void kh_simulate_c(
        int nx, int ny,
        double lx, double ly,
        double dt, int steps,
        double re, double u0, double amp, int mode,
        const double* init_omega,   // may be nullptr -> analytic IC
        int ftype, double famp,     // forcing (2026-10-10): 0 none, 1 mode, 2 band
        int fkx, int fky,
        double fkmin, double fkmax,
        double fomega0, double drag, double t_start,
        double* out_u,
        double* out_v,
        double* out_omega,
        double* out_psi,
        double* out_ke,
        double* out_enstrophy,
        double* out_max_vort,
        double* out_div_rms,
        double* out_cfl_peak
    );
}

namespace khsim {

struct SimRequest {
    int    nx = 128, ny = 64;
    double lx = 1.0, ly = 0.5;
    double dt = 0.001;
    int    steps = 100;
    double re = 1000.0, u0 = 1.0, amp = 0.01;
    int    mode = 2;
    std::optional<std::vector<double>> init_omega;
    // Forcing and continuation (2026-10-10, ruling 8): computed by the kernel.
    std::string forcing_type = "none";
    double forcing_amplitude = 0.0;
    int    forcing_kx = 1, forcing_ky = 1;
    double forcing_kmin = 0.0, forcing_kmax = 0.0, forcing_omega0 = 0.0;
    double linear_drag = 0.0;
    double t_start = 0.0;
};

struct Diagnostics {
    double kinetic_energy, enstrophy, max_vorticity, divergence_rms, cfl_peak;
};

struct SimResult {
    std::string backend = "fortran-fphttpapp";
    std::string language = "Fortran";
    int    steps_completed;
    double t_final;
    int    nx, ny;
    std::vector<double> u, v, omega, psi;
    Diagnostics diag;
    double compute_ms;
};

SimResult simulate(const SimRequest& req);

} // namespace khsim
