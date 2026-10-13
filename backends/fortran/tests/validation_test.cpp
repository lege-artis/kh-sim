// validation_test.cpp — KH-006 Fortran backend acceptance test
//
// Validates the Fortran kernel (through the C++ shim) against the oracle golden
// vector V1 (shared/golden-vectors/v1-canonical-reference.json, class A2), which
// replaced the retired shared/physics/kh_reference_output.json on 2026-10-08.
// The request is read from V1's "parameters"; the tolerance is V1's A2 band:
// |got - ref| <= max(scalar_abs * f, scalar_rel * f * |ref|), f = the
// "kh-fortran" scalar_factor override (1 when absent).
// Incompressibility: divergence_rms < 1e-10.
//
// Usage: kh-sim-fortran-validate <path/to/v1-canonical-reference.json>

#include "kh_shim.hpp"
#include <nlohmann/json.hpp>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>

using json = nlohmann::json;

static const char* kBackend = "kh-fortran";

int main(int argc, char* argv[]) {
    std::string path = "../../shared/golden-vectors/v1-canonical-reference.json";
    if (argc >= 2) path = argv[1];

    std::ifstream f(path);
    if (!f) {
        std::cerr << "Cannot open golden vector: " << path << "\n";
        return 1;
    }
    auto v = json::parse(f);
    const json& p    = v["parameters"];
    const json& rd   = v["diagnostics"];
    const json& tols = v["tolerances"];

    double factor = 1.0;
    if (tols.contains("backend_overrides") && tols["backend_overrides"].contains(kBackend))
        factor = tols["backend_overrides"][kBackend].value("scalar_factor", 1.0);
    const double s_abs = tols["scalar_abs"].get<double>() * factor;
    const double s_rel = tols["scalar_rel"].get<double>() * factor;

    khsim::SimRequest req;
    req.nx    = p["grid_nx"];
    req.ny    = p["grid_ny"];
    req.steps = p["nsteps"];
    req.lx    = p["domain_lx"];
    req.ly    = p["domain_ly"];
    req.dt    = p["dt"];
    req.re    = p["reynolds_number"];
    req.u0    = p["velocity_shear"];
    req.amp   = p["perturbation_amplitude"];
    req.mode  = p["perturbation_mode"];

    auto r = khsim::simulate(req);
    const auto& d = r.diag;

    bool pass = true;
    std::cout << "--- Fortran vs golden vector V1 (A2 band, factor " << factor << ") ---\n";
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
    std::cout << "compute_ms     : " << std::defaultfloat << r.compute_ms << "\n";
    if (d.divergence_rms >= 1e-10) {
        std::cerr << "FAIL divergence_rms: " << d.divergence_rms << "\n"; pass = false;
    }
    if ((int)r.u.size() != req.nx * req.ny) {
        std::cerr << "FAIL field shape\n"; pass = false;
    }

    std::cout << (pass ? "\nAll tests PASSED\n" : "\nSome tests FAILED\n");
    return pass ? 0 : 1;
}
