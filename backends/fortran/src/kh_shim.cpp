// kh_shim.cpp — C++ wrapper around Fortran kh_simulate_c()
// Allocates output buffers, calls Fortran, packages SimResult.

#include "kh_shim.hpp"
#include <chrono>
#include <cmath>
#include <stdexcept>

namespace khsim {

SimResult simulate(const SimRequest& req) {
    using clock = std::chrono::high_resolution_clock;
    auto t0 = clock::now();

    int n = req.nx * req.ny;
    std::vector<double> u(n), v(n), omega(n), psi(n);
    double ke, enstrophy, max_vort, div_rms, cfl_peak;

    const double* init_ptr = req.init_omega
        ? req.init_omega->data()
        : nullptr;

    // Forcing kind for the kernel: none when the amplitude is 0 (as the oracle).
    int ftype = 0;
    if (req.forcing_amplitude != 0.0) {
        if (req.forcing_type == "mode") ftype = 1;
        else if (req.forcing_type == "band") ftype = 2;
    }

    kh_simulate_c(
        req.nx, req.ny,
        req.lx, req.ly,
        req.dt, req.steps,
        req.re, req.u0, req.amp, req.mode,
        init_ptr,
        ftype, req.forcing_amplitude,
        req.forcing_kx, req.forcing_ky,
        req.forcing_kmin, req.forcing_kmax,
        req.forcing_omega0, req.linear_drag, req.t_start,
        u.data(), v.data(), omega.data(), psi.data(),
        &ke, &enstrophy, &max_vort, &div_rms, &cfl_peak
    );

    double elapsed = std::chrono::duration<double, std::milli>(
        clock::now() - t0).count();

    // The kernel's loop runs max(steps, 0) times; report what was computed.
    const int steps_done = req.steps > 0 ? req.steps : 0;
    double t_final = std::round((req.t_start + steps_done * req.dt) * 1e10) / 1e10;

    return SimResult{
        "fortran-httplib", "Fortran",
        steps_done, t_final, req.nx, req.ny,
        std::move(u), std::move(v), std::move(omega), std::move(psi),
        {ke, enstrophy, max_vort, div_rms, cfl_peak},
        std::round(elapsed * 100.0) / 100.0,
    };
}

} // namespace khsim
