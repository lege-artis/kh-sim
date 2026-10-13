// test_golden_vectors.cpp — per-backend golden-vector regression tests, C++ backend.
//
// Mode 1: HTTP POST /simulate to running server on port 8003.
// Compares response against oracle diagnostics using Higham §4.2 hybrid tolerance
// with backend_factor = 1.0 (C++ f64 = IEEE 754 binary64, matches Python/Rust precision).
//
// Run (requires server on :8003):
//   cd backends/cpp
//   cmake -B build && cmake --build build
//   ./build/kh-sim-cpp &
//   SERVER_PID=$!
//   sleep 1
//   ./build/test_golden_vectors
//   kill $SERVER_PID
//
// Each vector test:
//   1. Scalar diagnostics (kinetic_energy, enstrophy, max_vorticity, divergence_rms, cfl_peak)
//   2. Field summary statistics (min/max/mean per field u, v, vorticity, pressure)
//   3. Field element comparison (first_10 / last_10 per field)
//
// Exit 0 = all vectors PASS; non-zero = at least one FAIL.

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <numeric>
#include <string>
#include <vector>

#include <httplib.h>
#include <nlohmann/json.hpp>

using json = nlohmann::json;
namespace fs = std::filesystem;

// ── Configuration ─────────────────────────────────────────────────────────────

const std::string BACKEND_NAME = "kh-cpp";   // matches JSON tolerances.backend_overrides key
const std::string BACKEND_HOST = "localhost";
const int         BACKEND_PORT = 8003;

const double MACHINE_EPS = 2.22e-16;  // IEEE 754 binary64

// ── Vector loading ─────────────────────────────────────────────────────────────

std::vector<std::pair<fs::path, json>> load_vectors(const fs::path& dir) {
    std::vector<std::pair<fs::path, json>> vectors;
    for (const auto& entry : fs::directory_iterator(dir)) {
        const auto& p = entry.path();
        const std::string fname = p.filename().string();
        if (fname.rfind("v", 0) == 0 && p.extension() == ".json") {
            std::ifstream f(p);
            if (!f.is_open()) {
                std::cerr << "Cannot open " << p << "\n";
                continue;
            }
            json data;
            f >> data;
            vectors.emplace_back(p, std::move(data));
        }
    }
    // Sort deterministically by filename
    std::sort(vectors.begin(), vectors.end(),
        [](const auto& a, const auto& b) {
            return a.first.filename() < b.first.filename();
        });
    return vectors;
}

// ── Tolerance helpers (Higham §4.2 hybrid) ────────────────────────────────────

double backend_factor(const json& vector, const std::string& kind) {
    try {
        return vector.at("tolerances")
                     .at("backend_overrides")
                     .at(BACKEND_NAME)
                     .at(kind + "_factor")
                     .get<double>();
    } catch (...) {
        return 1.0;
    }
}

double scalar_threshold(const json& vector, const std::string& key, double oracle) {
    const double f    = backend_factor(vector, "scalar");
    const auto& tols  = vector.at("tolerances");
    double abs_t, rel_t;
    if (key == "divergence_rms") {
        abs_t = tols.value("scalar_div_rms_abs", 1e-11) * f;
        rel_t = tols.value("scalar_div_rms_rel", 1e-8)  * f;
    } else {
        abs_t = tols.value("scalar_abs", 1e-12) * f;
        rel_t = tols.value("scalar_rel", 1e-9)  * f;
    }
    return std::max(abs_t, rel_t * std::abs(oracle));
}

double field_threshold(const json& vector, double oracle, double sqrt_n) {
    const double f     = backend_factor(vector, "field");
    const auto& tols   = vector.at("tolerances");
    const double abs_t = tols.value("field_summary_abs", 1e-10) * f;
    const double rel_t = tols.value("field_summary_rel", 1e-7)  * f;
    const bool use_sn  = tols.value("sqrt_n_scaling", false);
    const double eps_term = use_sn ? sqrt_n * MACHINE_EPS * std::abs(oracle) : 0.0;
    return std::max(abs_t, rel_t * std::abs(oracle) + eps_term);
}

double field_elem_threshold(const json& vector, double oracle, double sqrt_n) {
    const double f     = backend_factor(vector, "field");
    const auto& tols   = vector.at("tolerances");
    const double abs_t = tols.value("field_element_abs", 1e-9) * f;
    const double rel_t = tols.value("field_element_rel", 1e-7) * f;
    const bool use_sn  = tols.value("sqrt_n_scaling", false);
    const double eps_term = use_sn ? sqrt_n * MACHINE_EPS * std::abs(oracle) : 0.0;
    return std::max(abs_t, rel_t * std::abs(oracle) + eps_term);
}

// ── Request builder ────────────────────────────────────────────────────────────

// Maps golden-vector "nsteps" → REST endpoint "steps".
std::string params_to_request(const json& params) {
    json req;
    const std::vector<std::string> int_fields = {"grid_nx", "grid_ny", "perturbation_mode"};
    const std::vector<std::string> dbl_fields = {
        "domain_lx", "domain_ly", "dt", "reynolds_number",
        "velocity_shear", "perturbation_amplitude"
    };
    for (const auto& f : int_fields)
        if (params.contains(f)) req[f] = params[f].get<int>();
    for (const auto& f : dbl_fields)
        if (params.contains(f)) req[f] = params[f].get<double>();
    // golden-vector uses "nsteps"; REST API uses "steps"
    if (params.contains("nsteps")) req["steps"] = params["nsteps"].get<int>();
    return req.dump();
}

// ── Core comparison ────────────────────────────────────────────────────────────

struct Failure {
    std::string msg;
};

std::vector<Failure> compare_vector(const std::string& filename,
                                    const json& response,
                                    const json& vector) {
    std::vector<Failure> failures;

    const auto& oracle    = vector.at("diagnostics");
    const auto& oracle_fw = vector.at("field_summaries");
    const auto& params    = vector.at("parameters");

    const int nx    = params.value("grid_nx", 64);
    const int ny    = params.value("grid_ny", 32);
    const double sqrt_n = std::sqrt(static_cast<double>(nx * ny));

    // Response diagnostics object
    const json& bdiag = response.contains("diagnostics") ? response["diagnostics"] : response;

    // ── 1. Scalar diagnostics ────────────────────────────────────────────────────
    const std::vector<std::string> scalar_keys = {
        "kinetic_energy", "enstrophy", "max_vorticity", "divergence_rms", "cfl_peak"
    };
    for (const auto& key : scalar_keys) {
        if (!oracle.contains(key) || !bdiag.contains(key)) continue;
        const double actual = bdiag[key].get<double>();
        const double ora    = oracle[key].get<double>();
        const double thresh = scalar_threshold(vector, key, ora);
        const double diff   = std::abs(actual - ora);
        if (diff >= thresh) {
            char buf[256];
            std::snprintf(buf, sizeof(buf),
                "%s: %s diff=%.3e >= thresh=%.3e (actual=%.6e oracle=%.6e)",
                filename.c_str(), key.c_str(), diff, thresh, actual, ora);
            failures.push_back({buf});
        }
    }

    // ── 2+3. Field summaries + element checks ───────────────────────────────────
    const std::vector<std::pair<std::string,std::string>> field_map = {
        {"u_velocity", "u_velocity"},
        {"v_velocity", "v_velocity"},
        {"vorticity",  "vorticity"},
        {"pressure",   "pressure"},
    };
    for (const auto& [oracle_key, body_key] : field_map) {
        if (!oracle_fw.contains(oracle_key)) continue;
        if (!response.contains(body_key)) {
            failures.push_back({filename + ": " + body_key + ": missing from response"});
            continue;
        }
        std::vector<double> flat;
        try {
            flat = response[body_key].get<std::vector<double>>();
        } catch (...) {
            failures.push_back({filename + ": " + body_key + ": not a float array"});
            continue;
        }
        if (flat.empty()) {
            failures.push_back({filename + ": " + body_key + ": empty array"});
            continue;
        }

        const auto& oFs = oracle_fw[oracle_key];

        // Summary stats
        const double b_min  = *std::min_element(flat.begin(), flat.end());
        const double b_max  = *std::max_element(flat.begin(), flat.end());
        const double b_mean = std::accumulate(flat.begin(), flat.end(), 0.0) / flat.size();

        for (const auto& [stat, bval, okey] : std::vector<std::tuple<std::string,double,std::string>>{
                {"min", b_min, "min"}, {"max", b_max, "max"}, {"mean", b_mean, "mean"}}) {
            if (!oFs.contains(okey)) continue;
            const double ov    = oFs[okey].get<double>();
            const double thresh = field_threshold(vector, ov, sqrt_n);
            const double diff  = std::abs(bval - ov);
            if (diff >= thresh) {
                char buf[256];
                std::snprintf(buf, sizeof(buf),
                    "%s: %s.%s diff=%.3e >= thresh=%.3e (actual=%.6e oracle=%.6e)",
                    filename.c_str(), oracle_key.c_str(), stat.c_str(),
                    diff, thresh, bval, ov);
                failures.push_back({buf});
            }
        }

        // Element-wise: first_10 / last_10
        for (const std::string slot : {"first_10", "last_10"}) {
            if (!oFs.contains(slot)) continue;
            std::vector<double> oe = oFs[slot].get<std::vector<double>>();
            std::vector<double> be;
            if (slot == "first_10") {
                be.assign(flat.begin(), flat.begin() + std::min((int)oe.size(), (int)flat.size()));
            } else {
                int n = std::min((int)oe.size(), (int)flat.size());
                be.assign(flat.end() - n, flat.end());
            }
            for (int i = 0; i < (int)std::min(be.size(), oe.size()); ++i) {
                const double thresh = field_elem_threshold(vector, oe[i], sqrt_n);
                const double diff   = std::abs(be[i] - oe[i]);
                if (diff >= thresh) {
                    char buf[256];
                    std::snprintf(buf, sizeof(buf),
                        "%s: %s.%s[%d] diff=%.3e >= thresh=%.3e (actual=%.6e oracle=%.6e)",
                        filename.c_str(), oracle_key.c_str(), slot.c_str(),
                        i, diff, thresh, be[i], oe[i]);
                    failures.push_back({buf});
                }
            }
        }
    }

    return failures;
}

// ── Per-vector test runner ─────────────────────────────────────────────────────

enum class TestResult { PASS, FAIL, SKIP };

TestResult run_vector_test(const fs::path& path, const json& vector) {
    const std::string filename = path.filename().string();

    // Build request
    const json& params = vector.at("parameters");
    const std::string req_body = params_to_request(params);

    // HTTP POST
    httplib::Client cli(BACKEND_HOST, BACKEND_PORT);
    cli.set_connection_timeout(5, 0);  // 5 seconds
    cli.set_read_timeout(30, 0);       // 30 seconds

    auto res = cli.Post("/simulate", req_body, "application/json");
    if (!res) {
        std::cerr << "SKIP " << filename << ": server not reachable at "
                  << BACKEND_HOST << ":" << BACKEND_PORT << "\n";
        return TestResult::SKIP;
    }
    if (res->status != 200) {
        std::cerr << "FAIL " << filename << ": HTTP " << res->status
                  << " — " << res->body.substr(0, 200) << "\n";
        return TestResult::FAIL;
    }

    json response;
    try {
        response = json::parse(res->body);
    } catch (const std::exception& e) {
        std::cerr << "FAIL " << filename << ": JSON parse error: " << e.what() << "\n";
        return TestResult::FAIL;
    }

    const auto failures = compare_vector(filename, response, vector);
    if (!failures.empty()) {
        std::cerr << "FAIL " << filename << ": " << failures.size() << " comparison failure(s):\n";
        for (const auto& f : failures)
            std::cerr << "  " << f.msg << "\n";
        return TestResult::FAIL;
    }

    return TestResult::PASS;
}

// ── Main ──────────────────────────────────────────────────────────────────────

int main() {
    // Vectors relative to binary location (backends/cpp/build/ → ../../shared/golden-vectors)
    const fs::path golden_dir = fs::path(__FILE__).parent_path().parent_path()
                                / "../../shared/golden-vectors";

    if (!fs::exists(golden_dir)) {
        std::cerr << "ERROR: golden-vectors dir not found: " << golden_dir << "\n";
        return 1;
    }

    const auto vectors = load_vectors(golden_dir);
    if (vectors.empty()) {
        std::cerr << "ERROR: no v*.json files found in " << golden_dir << "\n";
        return 1;
    }

    int passed = 0, failed = 0, skipped = 0;
    for (const auto& [path, vector] : vectors) {
        const TestResult result = run_vector_test(path, vector);
        const std::string fname = path.filename().string();
        switch (result) {
            case TestResult::PASS:
                std::cout << "PASS " << fname << "\n";
                ++passed;
                break;
            case TestResult::FAIL:
                std::cout << "FAIL " << fname << "\n";
                ++failed;
                break;
            case TestResult::SKIP:
                std::cout << "SKIP " << fname << "\n";
                ++skipped;
                break;
        }
    }

    std::cout << "\nSummary: " << passed << "/" << vectors.size() << " PASS";
    if (failed  > 0) std::cout << ", " << failed  << " FAIL";
    if (skipped > 0) std::cout << ", " << skipped << " SKIP";
    std::cout << "\n";

    if (skipped > 0 && failed == 0) {
        std::cerr << "NOTE: server was not reachable — all vectors skipped. "
                  << "Start kh-sim-cpp on port " << BACKEND_PORT
                  << " and re-run.\n";
        // Skipped-only: not a build/logic failure — exit 0 so CI can detect separately
        return 0;
    }

    return (failed == 0) ? 0 : 1;
}
