#pragma once
// models.hpp -- REST API data structures (mirrors openapi.yaml)
// Phase 7f: added forcing fields (PHASE-7F-FORCING-API-DESIGN-v0.1.md)
#include <optional>
#include <string>
#include <vector>
#include <nlohmann/json.hpp>

namespace khsim {

// ---- Request ----------------------------------------------------------------

struct SimulationRequest {
    int    grid_nx                = 128;
    int    grid_ny                = 64;
    double domain_lx              = 1.0;
    double domain_ly              = 0.5;
    double dt                     = 0.001;
    int    steps                  = 100;
    double reynolds_number        = 1000.0;
    double velocity_shear         = 1.0;
    double perturbation_amplitude = 0.01;
    int    perturbation_mode      = 2;
    std::optional<std::vector<double>> initial_omega;

    // Phase 7f forcing extension (PHASE-7F-FORCING-API-DESIGN-v0.1.md)
    // Defaults reproduce pre-7f behaviour exactly (no forcing).
    std::string forcing_type      = "none";  // "none" | "mode" | "band"
    double forcing_amplitude      = 0.0;     // A
    int    forcing_kx             = 1;       // "mode" only: integer mode along x
    int    forcing_ky             = 1;       // "mode" only: integer mode along y
    double forcing_kmin           = 0.0;     // "band" only: inner radius
    double forcing_kmax           = 0.0;     // "band" only: outer radius
    double forcing_omega0         = 0.0;     // v0.2: phase-rotation rate scale
    double linear_drag            = 0.0;     // v0.2: Ekman linear drag coefficient
    double t_start                = 0.0;     // time of the initial field (2026-10-10)
};

inline void from_json(const nlohmann::json& j, SimulationRequest& r) {
    if (j.contains("grid_nx"))                r.grid_nx                = j["grid_nx"];
    if (j.contains("grid_ny"))                r.grid_ny                = j["grid_ny"];
    if (j.contains("domain_lx"))              r.domain_lx              = j["domain_lx"];
    if (j.contains("domain_ly"))              r.domain_ly              = j["domain_ly"];
    if (j.contains("dt"))                     r.dt                     = j["dt"];
    if (j.contains("steps"))                  r.steps                  = j["steps"];
    if (j.contains("reynolds_number"))        r.reynolds_number        = j["reynolds_number"];
    if (j.contains("velocity_shear"))         r.velocity_shear         = j["velocity_shear"];
    if (j.contains("perturbation_amplitude")) r.perturbation_amplitude = j["perturbation_amplitude"];
    if (j.contains("perturbation_mode"))      r.perturbation_mode      = j["perturbation_mode"];
    if (j.contains("initial_omega") && !j["initial_omega"].is_null())
        r.initial_omega = j["initial_omega"].get<std::vector<double>>();
    // forcing fields
    if (j.contains("forcing_type"))           r.forcing_type           = j["forcing_type"];
    if (j.contains("forcing_amplitude"))      r.forcing_amplitude      = j["forcing_amplitude"];
    if (j.contains("forcing_kx"))             r.forcing_kx             = j["forcing_kx"];
    if (j.contains("forcing_ky"))             r.forcing_ky             = j["forcing_ky"];
    if (j.contains("forcing_kmin"))           r.forcing_kmin           = j["forcing_kmin"];
    if (j.contains("forcing_kmax"))           r.forcing_kmax           = j["forcing_kmax"];
    if (j.contains("forcing_omega0"))         r.forcing_omega0         = j["forcing_omega0"];
    if (j.contains("linear_drag"))            r.linear_drag            = j["linear_drag"];
    if (j.contains("t_start"))                r.t_start                = j["t_start"];
}

// ---- Response ---------------------------------------------------------------

struct Diagnostics {
    double kinetic_energy  = 0.0;
    double enstrophy       = 0.0;
    double max_vorticity   = 0.0;
    double divergence_rms  = 0.0;
    double cfl_peak        = 0.0;  // max over all steps of dt*(max|u|/dx + max|v|/dy) -- KH-NM-5
};

NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(Diagnostics,
    kinetic_energy, enstrophy, max_vorticity, divergence_rms, cfl_peak)

struct SimulationResult {
    std::string          backend;
    std::string          language;
    int                  steps_completed;
    double               t_final;
    int                  grid_nx;
    int                  grid_ny;
    std::vector<double>  u_velocity;
    std::vector<double>  v_velocity;
    std::vector<double>  vorticity;
    std::vector<double>  pressure;
    Diagnostics          diagnostics;
    double               compute_time_ms;
    std::vector<std::string> warnings;   // KH-VAL-2 K5 (filled by the HTTP layer)
};

NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(SimulationResult,
    backend, language, steps_completed, t_final,
    grid_nx, grid_ny,
    u_velocity, v_velocity, vorticity, pressure,
    diagnostics, compute_time_ms, warnings)

} // namespace khsim
