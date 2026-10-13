// main.cpp -- KH-SIM C++ backend (KH-005)
// HTTP server: cpp-httplib 0.18 on port 8003
// Spec: kh-sim/shared/api/openapi.yaml; request contract KH-VAL-1 in kh_http_contract.hpp
//
// 2026-10-03 (finding v0.5, C-1..C-4): every request is validated before the solver
// runs. Until then a non-power-of-two or zero grid reached the radix-2 FFT, wrote out
// of bounds and killed the process; a JSON array body computed the defaults; every
// failure, the caller's or the backend's, was reported as 400 without a field; and a
// blown-up run returned 200 with null numbers.

#include "kh_http_contract.hpp"
#include "physics.hpp"

#include <httplib.h>
#include <nlohmann/json.hpp>

#include <iostream>
#include <stdexcept>

using json = nlohmann::json;

static const char* BACKEND_ID = "kh-cpp";

static void set_cors(httplib::Response& res) {
    res.set_header("Access-Control-Allow-Origin",  "*");
    res.set_header("Access-Control-Allow-Methods", "GET, POST, OPTIONS");
    res.set_header("Access-Control-Allow-Headers", "Content-Type");
}

static void send_error(httplib::Response& res, int status, const std::string& msg,
                       const char* key, const std::string& value) {
    res.status = status;
    res.set_content(json{{"error", msg}, {key, value}}.dump(), "application/json");
}

static khsim::SimulationRequest to_solver(const khcontract::Request& c) {
    khsim::SimulationRequest r;
    r.grid_nx                = c.grid_nx;
    r.grid_ny                = c.grid_ny;
    r.domain_lx              = c.domain_lx;
    r.domain_ly              = c.domain_ly;
    r.dt                     = c.dt;
    r.steps                  = c.steps;
    r.reynolds_number        = c.reynolds_number;
    r.velocity_shear         = c.velocity_shear;
    r.perturbation_amplitude = c.perturbation_amplitude;
    r.perturbation_mode      = c.perturbation_mode;
    r.initial_omega          = c.initial_omega;
    r.forcing_type           = c.forcing_type;
    r.forcing_amplitude      = c.forcing_amplitude;
    r.forcing_kx             = c.forcing_kx;
    r.forcing_ky             = c.forcing_ky;
    r.forcing_kmin           = c.forcing_kmin;
    r.forcing_kmax           = c.forcing_kmax;
    r.forcing_omega0         = c.forcing_omega0;
    r.linear_drag            = c.linear_drag;
    r.t_start                = c.t_start;
    return r;
}

static bool result_is_finite(const khsim::SimulationResult& r) {
    const auto& d = r.diagnostics;
    for (double x : {d.kinetic_energy, d.enstrophy, d.max_vorticity, d.divergence_rms, d.cfl_peak})
        if (!std::isfinite(x)) return false;
    return khcontract::all_finite(r.u_velocity) && khcontract::all_finite(r.v_velocity)
        && khcontract::all_finite(r.vorticity)  && khcontract::all_finite(r.pressure);
}

int main() {
    httplib::Server svr;

    // -- CORS preflight -------------------------------------------------------
    svr.Options(".*", [](const httplib::Request&, httplib::Response& res) {
        set_cors(res);
        res.status = 204;
    });

    // -- POST /simulate -------------------------------------------------------
    svr.Post("/simulate", [](const httplib::Request& req, httplib::Response& res) {
        set_cors(res);
        khsim::SimulationRequest sr;
        khcontract::Request creq;
        try {
            creq = khcontract::parse(req.body, /*radix2=*/true, /*forcing_supported=*/true);
            sr = to_solver(creq);
        } catch (const khcontract::Invalid& e) {
            send_error(res, 400, e.what(), "field", e.field);
            return;
        }
        try {
            std::cout << "[simulate] nx=" << sr.grid_nx
                      << " ny=" << sr.grid_ny
                      << " steps=" << sr.steps
                      << " Re=" << sr.reynolds_number << "\n";

            auto result = khsim::simulate(sr);
            if (!result_is_finite(result)) {
                send_error(res, 500, "computation produced a non-finite value "
                           "(dt beyond the stability bound?)", "backend", BACKEND_ID);
                return;
            }
            result.warnings = khcontract::warnings(creq, result.diagnostics.cfl_peak);
            json out = result;

            std::cout << "[simulate] done t=" << result.t_final
                      << " ke=" << result.diagnostics.kinetic_energy
                      << " ms=" << result.compute_time_ms << "\n";

            res.set_content(out.dump(), "application/json");
        } catch (const std::exception& e) {
            send_error(res, 500, e.what(), "backend", BACKEND_ID);
        } catch (...) {
            send_error(res, 500, "unknown error", "backend", BACKEND_ID);
        }
    });

    // -- GET /health ----------------------------------------------------------
    svr.Get("/health", [](const httplib::Request&, httplib::Response& res) {
        set_cors(res);
        json h{{"status", "ok"}, {"backend", "cpp-httplib"}, {"port", 8003}};
        res.set_content(h.dump(), "application/json");
    });

    // -- GET /info ------------------------------------------------------------
    svr.Get("/info", [](const httplib::Request&, httplib::Response& res) {
        set_cors(res);
        json info{
            {"backend",         "cpp-httplib"},
            {"language",        "C++"},
            {"framework",       "cpp-httplib 0.18"},
            {"fft_library",     "Cooley-Tukey (built-in)"},
            {"grid_constraint", "power-of-two"},
            {"build_date",      khcontract::build_date_iso()},
            {"port",            8003},
            {"openapi_spec",    "kh-sim/shared/api/openapi.yaml"},
        };
        res.set_content(info.dump(), "application/json");
    });

    constexpr int PORT = 8003;
    std::cout << "kh-sim C++ backend listening on http://0.0.0.0:" << PORT << "\n";
    svr.listen("0.0.0.0", PORT);
}
