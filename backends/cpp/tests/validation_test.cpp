// validation_test.cpp — KH-005 acceptance test
//
// Validates C++ diagnostics against the oracle golden vector V1
// (shared/golden-vectors/v1-canonical-reference.json, class A2), which replaced
// the retired shared/physics/kh_reference_output.json on 2026-10-08.
// The request is read from V1's "parameters"; the tolerance is V1's A2 band:
// |got - ref| <= max(scalar_abs * f, scalar_rel * f * |ref|), f = the
// "kh-cpp" scalar_factor override (1 when absent).
// Incompressibility: divergence_rms < 1e-10.
//
// Usage: kh-sim-cpp-validate <path/to/v1-canonical-reference.json>
// CTest: cmake --build build && ctest --test-dir build -V

#include "physics.hpp"

#include <nlohmann/json.hpp>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>

using json = nlohmann::json;

static const char* kBackend = "kh-cpp";

int main(int argc, char* argv[]) {
    std::string ref_path = "../../shared/golden-vectors/v1-canonical-reference.json";
    if (argc >= 2) ref_path = argv[1];

    std::ifstream f(ref_path);
    if (!f) {
        std::cerr << "Cannot open golden vector: " << ref_path << "\n";
        return 1;
    }
    json v = json::parse(f);
    const json& p    = v["parameters"];
    const json& rd   = v["diagnostics"];
    const json& tols = v["tolerances"];

    double factor = 1.0;
    if (tols.contains("backend_overrides") && tols["backend_overrides"].contains(kBackend))
        factor = tols["backend_overrides"][kBackend].value("scalar_factor", 1.0);
    const double s_abs = tols["scalar_abs"].get<double>() * factor;
    const double s_rel = tols["scalar_rel"].get<double>() * factor;

    khsim::SimulationRequest req;
    req.grid_nx                = p["grid_nx"];
    req.grid_ny                = p["grid_ny"];
    req.domain_lx              = p["domain_lx"];
    req.domain_ly              = p["domain_ly"];
    req.dt                     = p["dt"];
    req.steps                  = p["nsteps"];
    req.reynolds_number        = p["reynolds_number"];
    req.velocity_shear         = p["velocity_shear"];
    req.perturbation_amplitude = p["perturbation_amplitude"];
    req.perturbation_mode      = p["perturbation_mode"];

    auto result = khsim::simulate(req);
    const auto& d = result.diagnostics;

    bool pass = true;
    std::cout << "--- C++ vs golden vector V1 (A2 band, factor " << factor << ") ---\n";
    auto check = [&](const char* name, double got) {
        const double ref  = rd[name];
        const double diff = std::abs(got - ref);
        const double tol  = std::max(s_abs, s_rel * std::abs(ref));
        std::cout << std::left << std::setw(15) << name << ": got=" << std::setprecision(17) << got
                  << "  ref=" << ref << std::scientific << std::setprecision(3)
                  << "  |diff|=" << diff << "  band=" << tol << std::defaultfloat << "\n";
        if (!(diff <= tol)) {
            std::cerr << "FAIL " << name << ": |diff| " << diff << " > band " << tol << "\n";
            pass = false;
        }
    };
    check("kinetic_energy", d.kinetic_energy);
    check("enstrophy",      d.enstrophy);
    check("max_vorticity",  d.max_vorticity);
    std::cout << "divergence_rms : got=" << std::scientific << d.divergence_rms << "\n";
    if (d.divergence_rms >= 1e-10) {
        std::cerr << "FAIL divergence_rms: " << d.divergence_rms << " >= 1e-10\n";
        pass = false;
    }

    // ── Field shape check (nx*ny) ──────────────────────────────────────────────
    int expected_n = req.grid_nx * req.grid_ny;
    if ((int)result.u_velocity.size() != expected_n ||
        (int)result.vorticity.size()  != expected_n) {
        std::cerr << "FAIL field shape mismatch\n";
        pass = false;
    }

    if (pass) {
        std::cout << "\nAll tests PASSED\n";
        return 0;
    } else {
        std::cout << "\nSome tests FAILED\n";
        return 1;
    }
}
