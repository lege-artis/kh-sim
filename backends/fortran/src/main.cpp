// main.cpp -- KH-SIM Fortran backend HTTP server (KH-006)
// cpp-httplib on port 8004; physics delegated to Fortran via kh_shim
// Request contract KH-VAL-1: kh_http_contract.hpp (byte-identical copy of backends/cpp/src's)
//
// 2026-10-03 (finding v0.5, F-1..F-4): every request is validated before the
// Fortran kernel runs. A Fortran runtime error (the kernel is built with
// -fcheck=bounds) terminates the whole process, so nothing that can trip one may
// reach it: until then a zero or non-power-of-two grid, and ANY initial_omega,
// killed the server. Forcing parameters were silently dropped; this server does
// not implement forcing, so it now refuses a forcing request instead of computing
// an unforced flow and returning it as if forced.
//
// 2026-10-10 (ruling 8, with the shared CP-7F-B/D harness): the kernel computes
// forcing ("mode", "band" with rotating phases) and linear_drag, and continues a run
// from t_start; the refusal is gone. kh_forcing.f90 holds at most MAX_BAND_MODES
// (4096) band modes, so a larger band is a 400 here, never a kernel stop.

#include "kh_http_contract.hpp"
#include "kh_shim.hpp"

#include <httplib.h>
#include <nlohmann/json.hpp>

#include <cmath>
#include <iostream>
#include <stdexcept>

using json = nlohmann::json;

static const char* BACKEND_ID = "kh-fortran";

static void cors(httplib::Response& res) {
    res.set_header("Access-Control-Allow-Origin",  "*");
    res.set_header("Access-Control-Allow-Methods", "GET, POST, OPTIONS");
    res.set_header("Access-Control-Allow-Headers", "Content-Type");
}

static void send_error(httplib::Response& res, int status, const std::string& msg,
                       const char* key, const std::string& value) {
    res.status = status;
    res.set_content(json{{"error", msg}, {key, value}}.dump(), "application/json");
}

static khsim::SimRequest to_solver(const khcontract::Request& c) {
    khsim::SimRequest r;
    r.nx    = c.grid_nx;
    r.ny    = c.grid_ny;
    r.lx    = c.domain_lx;
    r.ly    = c.domain_ly;
    r.dt    = c.dt;
    r.steps = c.steps;
    r.re    = c.reynolds_number;
    r.u0    = c.velocity_shear;
    r.amp   = c.perturbation_amplitude;
    r.mode  = c.perturbation_mode;
    r.init_omega = c.initial_omega;
    r.forcing_type      = c.forcing_type;
    r.forcing_amplitude = c.forcing_amplitude;
    r.forcing_kx        = c.forcing_kx;
    r.forcing_ky        = c.forcing_ky;
    r.forcing_kmin      = c.forcing_kmin;
    r.forcing_kmax      = c.forcing_kmax;
    r.forcing_omega0    = c.forcing_omega0;
    r.linear_drag       = c.linear_drag;
    r.t_start           = c.t_start;
    return r;
}

// Band modes kh_forcing.f90 would build (half-plane, kmin <= |k| <= kmax).
static long band_modes(double kmin, double kmax) {
    long n = 0;
    const int km = static_cast<int>(std::floor(kmax));
    for (int ix = 0; ix <= km; ++ix)
        for (int iy = -km; iy <= km; ++iy) {
            if (ix == 0 && iy <= 0) continue;
            const double kr = std::sqrt(static_cast<double>(ix * ix + iy * iy));
            if (kr >= kmin && kr <= kmax) ++n;
        }
    return n;
}

static bool result_is_finite(const khsim::SimResult& r) {
    const auto& d = r.diag;
    for (double x : {d.kinetic_energy, d.enstrophy, d.max_vorticity, d.divergence_rms, d.cfl_peak})
        if (!std::isfinite(x)) return false;
    return khcontract::all_finite(r.u) && khcontract::all_finite(r.v)
        && khcontract::all_finite(r.omega) && khcontract::all_finite(r.psi);
}

int main() {
    httplib::Server svr;

    svr.Options(".*", [](const httplib::Request&, httplib::Response& res) {
        cors(res); res.status = 204;
    });

    svr.Post("/simulate", [](const httplib::Request& req, httplib::Response& res) {
        cors(res);
        khsim::SimRequest sr;
        khcontract::Request creq;
        try {
            creq = khcontract::parse(req.body, /*radix2=*/true, /*forcing_supported=*/true);
            if (creq.forcing_type == "band" && creq.forcing_amplitude != 0.0
                && band_modes(creq.forcing_kmin, creq.forcing_kmax) > 4096)
                throw khcontract::Invalid("forcing_kmax", "this backend's band forcing holds at most 4096 modes");
            sr = to_solver(creq);
        } catch (const khcontract::Invalid& e) {
            send_error(res, 400, e.what(), "field", e.field);
            return;
        }
        try {
            std::cout << "[simulate] nx=" << sr.nx << " ny=" << sr.ny
                      << " steps=" << sr.steps << " Re=" << sr.re << "\n";
            auto r = khsim::simulate(sr);
            if (!result_is_finite(r)) {
                send_error(res, 500, "computation produced a non-finite value "
                           "(dt beyond the stability bound?)", "backend", BACKEND_ID);
                return;
            }
            json out;
            out["backend"]         = r.backend;
            out["language"]        = r.language;
            out["steps_completed"] = r.steps_completed;
            out["t_final"]         = r.t_final;
            out["grid_nx"]         = r.nx;
            out["grid_ny"]         = r.ny;
            out["u_velocity"]      = r.u;
            out["v_velocity"]      = r.v;
            out["vorticity"]       = r.omega;
            out["pressure"]        = r.psi;
            out["diagnostics"] = {
                {"kinetic_energy",  r.diag.kinetic_energy},
                {"enstrophy",       r.diag.enstrophy},
                {"max_vorticity",   r.diag.max_vorticity},
                {"divergence_rms",  r.diag.divergence_rms},
                {"cfl_peak",        r.diag.cfl_peak},
            };
            out["compute_time_ms"] = r.compute_ms;
            out["warnings"]        = khcontract::warnings(creq, r.diag.cfl_peak);
            std::cout << "[simulate] done t=" << r.t_final
                      << " ke=" << r.diag.kinetic_energy
                      << " ms=" << r.compute_ms << "\n";
            res.set_content(out.dump(), "application/json");
        } catch (const std::exception& e) {
            send_error(res, 500, e.what(), "backend", BACKEND_ID);
        } catch (...) {
            send_error(res, 500, "unknown error", "backend", BACKEND_ID);
        }
    });

    svr.Get("/health", [](const httplib::Request&, httplib::Response& res) {
        cors(res);
        res.set_content(
            json{{"status","ok"},{"backend","fortran-httplib"},{"port",8004}}.dump(),
            "application/json");
    });

    svr.Get("/info", [](const httplib::Request&, httplib::Response& res) {
        cors(res);
        res.set_content(json{
            {"backend",         "fortran-httplib"},
            {"language",        "Fortran"},
            {"framework",       "cpp-httplib 0.18 (C-interop shim)"},
            {"fft_library",     "Cooley-Tukey (native Fortran)"},
            {"grid_constraint", "power-of-two"},
            {"forcing",         "mode, band (rotating phases), linear_drag; at most 4096 band modes"},
            {"build_date",      khcontract::build_date_iso()},
            {"port",            8004},
            {"openapi_spec",    "kh-sim/shared/api/openapi.yaml"},
        }.dump(), "application/json");
    });

    constexpr int PORT = 8004;
    std::cout << "kh-sim Fortran backend listening on http://0.0.0.0:" << PORT << "\n";
    svr.listen("0.0.0.0", PORT);
}
