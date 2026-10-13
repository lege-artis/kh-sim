// canonical_physics.cpp -- Phase 7f canonical-physics testbeds, C++ backend.
//
// Every CP-7F testbed now runs once, for the oracle and all five backends, from a
// shared harness in tests/integration (ruling 2, 2026-10-07):
//   CP-7F-A  test_canonical_physics_cp7fa.py   (MASTER-BRIEF s3.1, E9; moved 2026-10-12)
//   CP-7F-B  test_canonical_physics_cp7fbd.py  (s3.2, E7; moved 2026-10-10)
//   CP-7F-C  test_canonical_physics_cp7fc.py   (s3.3, E6; moved 2026-10-07)
//   CP-7F-D  test_canonical_physics_cp7fbd.py  (s3.4, E7; moved 2026-10-10)
//   CP-7F-E  the Pascal-binding anchor (s3.5, E8): asserted by the oracle and by the
//            Pascal and Fortran native suites; C++ takes no part.
// The C++ CP-7F-A that stood here ran a single 256^2 sweep of 12 modes on the double
// layer at delta = Ly/20, where the two layers couple (assessment 2026-09-15, s7), with
// the pre-correction gate [0.180, 0.200]. Its text is in the git history of this file.
//
// Usage: ./canonical_physics [A|B|C|D|E|ALL]   -- prints where each testbed runs; exit 0.

#include <iostream>
#include <map>
#include <string>

int main(int argc, char** argv) {
    static const std::map<std::string, std::string> moved = {
        {"A", "CP-7F-A: SKIP -- moved to tests/integration/test_canonical_physics_cp7fa.py (MASTER s3.1 E9)"},
        {"B", "CP-7F-B: SKIP -- moved to tests/integration/test_canonical_physics_cp7fbd.py (MASTER s3.2 E7)"},
        {"C", "CP-7F-C: SKIP -- moved to tests/integration/test_canonical_physics_cp7fc.py (MASTER s3.3 E6)"},
        {"D", "CP-7F-D: SKIP -- moved to tests/integration/test_canonical_physics_cp7fbd.py (MASTER s3.4 E7)"},
        {"E", "CP-7F-E: SKIP -- Pascal-binding anchor (MASTER s3.5 E8): the oracle, Pascal and Fortran assert it"},
    };
    std::string mode = argc >= 2 ? argv[1] : "ALL";
    if (mode == "B-GATE") mode = "B";
    std::cout << "=== kh-cpp phase 7f canonical physics testbeds ===\n";
    if (mode == "ALL") {
        for (const auto& [k, line] : moved) std::cout << line << "\n";
        return 0;
    }
    auto it = moved.find(mode);
    if (it == moved.end()) { std::cerr << "Unknown mode: " << mode << "\n"; return 2; }
    std::cout << it->second << "\n";
    return 0;
}
