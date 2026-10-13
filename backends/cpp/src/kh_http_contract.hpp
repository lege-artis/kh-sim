#pragma once
// kh_http_contract.hpp -- POST /simulate request contract for the C++-hosted servers
//
// Used by backends/cpp and backends/fortran. Each Docker build context copies only
// its own src/, so the file exists twice and the two copies are byte-identical;
// tests/integration/test_conformance.py locks the behaviour of both servers.
//
// Contract KH-VAL-1 (kh-sim/_audit/FINDING-PASCAL-HTTP-PARAMETER-DROP-2026-09-16.md
// v0.5; first implemented in backends/pascal/src/kh_server.pas):
//   - body must be a JSON object                          -> 400 field "body"
//   - a key present with the wrong JSON type               -> 400 field <key>
//   - grid_nx, grid_ny in [1, 65536]; power of two when the
//     FFT is radix-2; grid_nx*grid_ny <= 4194304          -> 400 field grid_nx|grid_ny
//   - domain_lx, domain_ly, reynolds_number non-zero      -> 400 field <key>
//   - initial_omega (legacy alias init_omega, not both)
//     has exactly grid_nx*grid_ny entries                 -> 400 field initial_omega
//   - forcing_type in {none, mode, band}; band needs
//     domain_lx == domain_ly                              -> 400 field forcing_type
//   - a computation that fails or produces a non-finite
//     number                                             -> 500 {error, backend}
//
// Contract KH-VAL-2 (kh-sim/_audit/KH-SIM-MODELLING-RULES-v0.1.md, ruling 2026-10-07)
// adds, after the KH-VAL-1 checks:
//   K3 outside the model (computable, but no physical state of it):
//     domain_lx, domain_ly, reynolds_number < 0; dt < 0; steps < 0;
//     linear_drag < 0; perturbation_mode < 0; band forcing with
//     forcing_kmin > forcing_kmax                         -> 400 field <key>
//   K4 not representable (analytic initial field only, amplitude != 0):
//     perturbation_mode 0, or beyond the 2/3-rule band
//     (3*mode > 2*(grid_nx/2))                            -> 400 field perturbation_mode
//   K4b (ruling 2026-10-09) a non-zero perturbation_amplitude sent together with
//     initial_omega, which replaces the analytic field: the amplitude would be
//     silently ignored                                    -> 400 field perturbation_amplitude
//   K4 forcing (amplitude != 0; ruled 2026-10-07, implemented with the CP-7F-B/D
//     harness 2026-10-10): mode forcing_kx / forcing_ky beyond the 2/3-rule band
//     (3*|k| > 2*(grid/2))                                -> 400 field forcing_kx|forcing_ky
//     band reaching beyond it (3*floor(forcing_kmax) > 2*(min(grid_nx, grid_ny)/2))
//                                                         -> 400 field forcing_kmax
//   K3 forcing: a band with no integer mode kmin <= |k| <= kmax (nothing would be
//     forced, R1)                                         -> 400 field forcing_kmax
//   K5 computed but possibly not resolved: "warnings" in the response (warnings())
//   K7 unknown keys: logged to stderr (refused in a later version, after a census)
//
// t_start (2026-10-10): the time of the initial field. Time-dependent forcing is
// evaluated at t_start + t, and t_final = t_start + steps*dt, so a run continued from
// a returned field (initial_omega = vorticity, t_start = t_final) keeps its phase.

#include <nlohmann/json.hpp>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdint>
#include <cstdlib>
#include <optional>
#include <stdexcept>
#include <string>
#include <vector>

namespace khcontract {

using json = nlohmann::json;

constexpr std::int64_t MAX_AXIS  = 65536;
constexpr std::int64_t MAX_CELLS = 4194304;

// A request the backend refuses: HTTP 400 {"error": what(), "field": field}.
struct Invalid : std::runtime_error {
    std::string field;
    Invalid(std::string f, const std::string& msg) : std::runtime_error(msg), field(std::move(f)) {}
};

struct Request {
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
    std::string forcing_type      = "none";
    double forcing_amplitude      = 0.0;
    int    forcing_kx             = 1;
    int    forcing_ky             = 1;
    double forcing_kmin           = 0.0;
    double forcing_kmax           = 0.0;
    double forcing_omega0         = 0.0;
    double linear_drag            = 0.0;
    double t_start                = 0.0;
};

inline void read_int(const json& j, const char* key, int& dst) {
    auto it = j.find(key);
    if (it == j.end() || it->is_null()) return;
    std::int64_t v = 0;
    if (it->is_number_integer()) {
        v = it->get<std::int64_t>();
    } else if (it->is_number_float()) {
        double d = it->get<double>();
        if (!std::isfinite(d) || d != std::floor(d) || std::abs(d) > 2147483647.0)
            throw Invalid(key, std::string(key) + " must be an integer");
        v = static_cast<std::int64_t>(d);
    } else {
        throw Invalid(key, std::string(key) + " must be an integer, got " + it->type_name());
    }
    if (v < -2147483647LL || v > 2147483647LL)
        throw Invalid(key, std::string(key) + " is outside the 32-bit integer range");
    dst = static_cast<int>(v);
}

inline void read_double(const json& j, const char* key, double& dst) {
    auto it = j.find(key);
    if (it == j.end() || it->is_null()) return;
    if (!it->is_number())
        throw Invalid(key, std::string(key) + " must be a number, got " + it->type_name());
    dst = it->get<double>();
}

inline void read_string(const json& j, const char* key, std::string& dst) {
    auto it = j.find(key);
    if (it == j.end() || it->is_null()) return;
    if (!it->is_string())
        throw Invalid(key, std::string(key) + " must be a string, got " + it->type_name());
    dst = it->get<std::string>();
}

inline std::optional<std::vector<double>> read_field(const json& j, const char* key) {
    auto it = j.find(key);
    if (it == j.end() || it->is_null()) return std::nullopt;
    if (!it->is_array())
        throw Invalid("initial_omega", std::string(key) + " must be an array of numbers");
    if (static_cast<std::int64_t>(it->size()) > MAX_CELLS)
        throw Invalid("initial_omega", std::string(key) + " has more entries than the "
                      + std::to_string(MAX_CELLS) + "-cell limit");
    std::vector<double> v;
    v.reserve(it->size());
    for (const auto& x : *it) {
        if (!x.is_number())
            throw Invalid("initial_omega", std::string(key) + " must contain numbers only");
        v.push_back(x.get<double>());
    }
    return v;
}

inline bool is_power_of_two(int n) { return n > 0 && (n & (n - 1)) == 0; }

inline void check_axis(const char* key, int n, bool radix2) {
    if (n < 1 || n > MAX_AXIS)
        throw Invalid(key, std::string(key) + " = " + std::to_string(n)
                      + " is out of range [1, " + std::to_string(MAX_AXIS) + "]");
    if (radix2 && !is_power_of_two(n))
        throw Invalid(key, std::string(key) + " = " + std::to_string(n)
                      + " is not a power of two; this backend uses a radix-2 FFT");
}

// Parse and validate a POST /simulate body. Throws Invalid on any request error.
// radix2: the backend's FFT handles power-of-two lengths only.
// forcing_supported: false refuses any forcing request instead of ignoring it.
inline Request parse(const std::string& body, bool radix2, bool forcing_supported) {
    json j;
    try {
        j = json::parse(body);
    } catch (const json::exception&) {
        throw Invalid("body", "request body is not valid JSON");
    }
    if (!j.is_object())
        throw Invalid("body", "request body must be a JSON object");

    Request r;
    read_int   (j, "grid_nx",                r.grid_nx);
    read_int   (j, "grid_ny",                r.grid_ny);
    read_double(j, "domain_lx",              r.domain_lx);
    read_double(j, "domain_ly",              r.domain_ly);
    read_double(j, "dt",                     r.dt);
    read_int   (j, "steps",                  r.steps);
    read_double(j, "reynolds_number",        r.reynolds_number);
    read_double(j, "velocity_shear",         r.velocity_shear);
    read_double(j, "perturbation_amplitude", r.perturbation_amplitude);
    read_int   (j, "perturbation_mode",      r.perturbation_mode);
    read_string(j, "forcing_type",           r.forcing_type);
    read_double(j, "forcing_amplitude",      r.forcing_amplitude);
    read_int   (j, "forcing_kx",             r.forcing_kx);
    read_int   (j, "forcing_ky",             r.forcing_ky);
    read_double(j, "forcing_kmin",           r.forcing_kmin);
    read_double(j, "forcing_kmax",           r.forcing_kmax);
    read_double(j, "forcing_omega0",         r.forcing_omega0);
    read_double(j, "linear_drag",            r.linear_drag);
    read_double(j, "t_start",                r.t_start);

    auto canonical = read_field(j, "initial_omega");
    auto legacy    = read_field(j, "init_omega");
    if (canonical && legacy)
        throw Invalid("initial_omega", "send initial_omega or the legacy init_omega, not both");
    r.initial_omega = canonical ? std::move(canonical) : std::move(legacy);

    check_axis("grid_nx", r.grid_nx, radix2);
    check_axis("grid_ny", r.grid_ny, radix2);
    const std::int64_t cells = static_cast<std::int64_t>(r.grid_nx) * r.grid_ny;
    if (cells > MAX_CELLS)
        throw Invalid("grid_nx", "grid " + std::to_string(r.grid_nx) + "x" + std::to_string(r.grid_ny)
                      + " has " + std::to_string(cells) + " cells, more than the "
                      + std::to_string(MAX_CELLS) + "-cell limit");
    if (r.domain_lx == 0.0)       throw Invalid("domain_lx", "domain_lx must be non-zero");
    if (r.domain_ly == 0.0)       throw Invalid("domain_ly", "domain_ly must be non-zero");
    if (r.reynolds_number == 0.0) throw Invalid("reynolds_number", "reynolds_number must be non-zero");
    if (r.initial_omega && static_cast<std::int64_t>(r.initial_omega->size()) != cells)
        throw Invalid("initial_omega", "initial_omega has " + std::to_string(r.initial_omega->size())
                      + " entries; grid " + std::to_string(r.grid_nx) + "x" + std::to_string(r.grid_ny)
                      + " needs " + std::to_string(cells));

    if (r.forcing_type != "none" && r.forcing_type != "mode" && r.forcing_type != "band")
        throw Invalid("forcing_type", "unknown forcing_type '" + r.forcing_type
                      + "'; expected none, mode or band");
    if (!forcing_supported) {
        if (r.forcing_type != "none" && r.forcing_amplitude != 0.0)
            throw Invalid("forcing_type", "this backend's HTTP server does not implement forcing");
        if (r.linear_drag != 0.0)
            throw Invalid("linear_drag", "this backend's HTTP server does not implement linear_drag");
    } else if (r.forcing_type == "band" && r.forcing_amplitude != 0.0
               && std::abs(r.domain_lx - r.domain_ly) > 1e-10 * std::abs(r.domain_lx)) {
        throw Invalid("forcing_type", "band forcing requires domain_lx == domain_ly");
    }

    // KH-VAL-2, K3: outside the model.
    if (r.domain_lx < 0.0) throw Invalid("domain_lx", "domain_lx must be positive");
    if (r.domain_ly < 0.0) throw Invalid("domain_ly", "domain_ly must be positive");
    if (r.reynolds_number < 0.0)
        throw Invalid("reynolds_number", "reynolds_number must be positive (a negative value is negative viscosity)");
    if (r.dt < 0.0) throw Invalid("dt", "dt must be >= 0 (backward integration of a viscous flow is ill-posed)");
    if (r.steps < 0) throw Invalid("steps", "steps must be >= 0");
    if (r.linear_drag < 0.0) throw Invalid("linear_drag", "linear_drag must be >= 0 (negative drag injects energy)");
    if (r.perturbation_mode < 0) throw Invalid("perturbation_mode", "perturbation_mode must be >= 0");
    if (r.forcing_type == "band" && r.forcing_amplitude != 0.0 && r.forcing_kmin > r.forcing_kmax)
        throw Invalid("forcing_kmin", "forcing_kmin > forcing_kmax: the forcing band is empty");
    // KH-VAL-2, K4 (forcing): the forced modes must lie inside the 2/3-rule band.
    if (r.forcing_type == "mode" && r.forcing_amplitude != 0.0) {
        if (3LL * std::llabs(r.forcing_kx) > 2LL * (r.grid_nx / 2))
            throw Invalid("forcing_kx", "forcing_kx " + std::to_string(r.forcing_kx)
                          + " lies beyond the 2/3-rule band of grid_nx " + std::to_string(r.grid_nx));
        if (3LL * std::llabs(r.forcing_ky) > 2LL * (r.grid_ny / 2))
            throw Invalid("forcing_ky", "forcing_ky " + std::to_string(r.forcing_ky)
                          + " lies beyond the 2/3-rule band of grid_ny " + std::to_string(r.grid_ny));
    }
    if (r.forcing_type == "band" && r.forcing_amplitude != 0.0) {
        const double kfl = std::floor(r.forcing_kmax);
        if (3.0 * kfl > 2.0 * (std::min(r.grid_nx, r.grid_ny) / 2))
            throw Invalid("forcing_kmax", "forcing_kmax reaches beyond the 2/3-rule band of the "
                          + std::to_string(r.grid_nx) + "x" + std::to_string(r.grid_ny) + " grid");
        // K3: the band must contain an integer mode (half-plane, kmin <= |k| <= kmax).
        bool any = false;
        const int km = static_cast<int>(kfl);
        for (int ix = 0; ix <= km && !any; ++ix)
            for (int iy = -km; iy <= km; ++iy) {
                if (ix == 0 && iy <= 0) continue;
                const double kr = std::sqrt(static_cast<double>(ix * ix + iy * iy));
                if (kr >= r.forcing_kmin && kr <= r.forcing_kmax) { any = true; break; }
            }
        if (!any)
            throw Invalid("forcing_kmax", "the forcing band contains no integer mode: nothing would be forced");
    }
    // KH-VAL-2, K4: the analytic perturbation must exist on the grid, inside the 2/3-rule band.
    if (!r.initial_omega && r.perturbation_amplitude != 0.0) {
        if (r.perturbation_mode == 0)
            throw Invalid("perturbation_mode", "perturbation_mode 0 with a non-zero amplitude perturbs nothing");
        if (3LL * r.perturbation_mode > 2LL * (r.grid_nx / 2))
            throw Invalid("perturbation_mode", "perturbation_mode " + std::to_string(r.perturbation_mode)
                          + " lies beyond the 2/3-rule band of grid_nx " + std::to_string(r.grid_nx)
                          + " (largest resolved mode " + std::to_string(2 * (r.grid_nx / 2) / 3) + ")");
    }
    // KH-VAL-2, K4b: initial_omega replaces the analytic field; an amplitude sent with it would be ignored.
    if (r.initial_omega && j.contains("perturbation_amplitude") && r.perturbation_amplitude != 0.0)
        throw Invalid("perturbation_amplitude",
                      "perturbation_amplitude has no effect with initial_omega (send 0 or omit it)");

    // KH-VAL-2, K7: unknown keys are logged now and refused once no client sends one.
    static const char* known[] = {
        "grid_nx", "grid_ny", "domain_lx", "domain_ly", "dt", "steps", "reynolds_number",
        "velocity_shear", "perturbation_amplitude", "perturbation_mode", "initial_omega",
        "init_omega", "forcing_type", "forcing_amplitude", "forcing_kx", "forcing_ky",
        "forcing_kmin", "forcing_kmax", "forcing_omega0", "linear_drag", "t_start"};
    for (auto it = j.begin(); it != j.end(); ++it) {
        bool found = false;
        for (const char* k : known) if (it.key() == k) { found = true; break; }
        if (!found)
            std::fprintf(stderr, "KH-VAL-2 K7: unknown key '%s' ignored\n", it.key().c_str());
    }
    return r;
}

// KH-VAL-2, K5: computed, but possibly not resolved. Same conditions, order and
// prefixes in every backend and the oracle (kh_physics.request_warnings).
inline std::vector<std::string> warnings(const Request& r, double cfl_peak) {
    std::vector<std::string> out;
    char buf[160];
    const double nu = std::abs(r.velocity_shear) / r.reynolds_number;
    if (cfl_peak > 1.0) {
        std::snprintf(buf, sizeof buf, "cfl_advective: peak Courant number %.3g exceeds 1", cfl_peak);
        out.emplace_back(buf);
    }
    const double dmin = std::min(r.domain_lx / r.grid_nx, r.domain_ly / r.grid_ny);
    const double diff = 2.0 * nu * r.dt / (dmin * dmin);
    if (diff > 1.0) {
        std::snprintf(buf, sizeof buf, "cfl_diffusive: 2*nu*dt/min(dx,dy)^2 = %.3g exceeds 1", diff);
        out.emplace_back(buf);
    }
    if (nu == 0.0)
        out.emplace_back("inviscid: nu = 0, no physical dissipation; enstrophy accumulates at the 2/3 cutoff");
    return out;
}

// True when every value is finite; used to turn a blown-up run into a 500.
inline bool all_finite(const std::vector<double>& v) {
    for (double x : v) if (!std::isfinite(x)) return false;
    return true;
}

// __DATE__ ("Oct  3 2026") as ISO "2026-10-03": the image's build date in /info,
// so a deploy can tell a fresh build from a stale image.
inline std::string build_date_iso() {
    static const char* months = "JanFebMarAprMayJunJulAugSepOctNovDec";
    const std::string d = __DATE__;
    int m = 0;
    for (int i = 0; i < 12; ++i)
        if (d.compare(0, 3, months + 3 * i, 3) == 0) m = i + 1;
    std::string day = d.substr(4, 2);
    if (day[0] == ' ') day[0] = '0';
    char buf[11];
    std::snprintf(buf, sizeof buf, "%s-%02d-%s", d.substr(7, 4).c_str(), m, day.c_str());
    return buf;
}

} // namespace khcontract
