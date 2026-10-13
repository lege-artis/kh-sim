#pragma once
// physics.hpp — KH instability solver declarations
// Reference: kh-sim/shared/physics/KH-PHYSICS.md
#include "models.hpp"

namespace khsim {
    SimulationResult simulate(const SimulationRequest& req);
} // namespace khsim
