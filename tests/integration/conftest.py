# conftest.py — shared fixtures for KH-SIM integration test suite
# Task: KH-018 | Device: ThinkPad (LDE stack must be running)
#
# Environment overrides (all optional):
#   KH_RUST_URL     default http://localhost:8001
#   KH_SCALA_URL    default http://localhost:8002
#   KH_CPP_URL      default http://localhost:8003
#   KH_FORTRAN_URL  default http://localhost:8004
#   KH_PASCAL_URL   default http://localhost:8005
#   KH_LOG_URL      default http://localhost:8006
#
# Run: pytest kh-sim/tests/integration/ -v
#      pytest kh-sim/tests/integration/ -v -k "health"     # subset

import json
import os
import uuid
from pathlib import Path

import pytest
import requests

# ── Backend manifest ──────────────────────────────────────────────────────────

BACKEND_META = {
    "kh-rust":    {"port": 8001, "language": "Rust",    "framework": "axum"},
    "kh-scala":   {"port": 8002, "language": "Scala",   "framework": "http4s"},
    "kh-cpp":     {"port": 8003, "language": "C++",     "framework": "cpp-httplib"},
    "kh-fortran": {"port": 8004, "language": "Fortran", "framework": "C-interop"},
    "kh-pascal":  {"port": 8005, "language": "Pascal",  "framework": "fphttpapp"},
}

BACKEND_NAMES = list(BACKEND_META.keys())


def _backend_url(name: str) -> str:
    env_key = f"KH_{name.replace('kh-', '').upper()}_URL"
    return os.environ.get(env_key, f"http://localhost:{BACKEND_META[name]['port']}")


BACKEND_URLS: dict[str, str] = {name: _backend_url(name) for name in BACKEND_NAMES}
LOG_SERVICE_URL: str = os.environ.get("KH_LOG_URL", "http://localhost:8006")

# ── Reference data ────────────────────────────────────────────────────────────
# The oracle golden vector V1 (class A2, _audit/KH-SIM-MODELLING-RULES-v0.1.md
# s4). It replaced shared/physics/kh_reference_output.json, retired 2026-10-08:
# a pre-5c copy (A3, situations S2 + S9) that every backend was graded against
# at +-5 %, which hid the 8.5e-10 dealiasing difference. One source of truth:
# the request, the expected diagnostics and the tolerance all come from V1.

_REF_PATH = (
    Path(__file__).parents[2] / "shared" / "golden-vectors" / "v1-canonical-reference.json"
)

with open(_REF_PATH) as _f:
    REFERENCE_VECTOR: dict = json.load(_f)

REFERENCE_DIAGNOSTICS: dict[str, float] = REFERENCE_VECTOR["diagnostics"]

_P = REFERENCE_VECTOR["parameters"]
REFERENCE_PARAMS: dict = {
    "grid_nx":                 _P["grid_nx"],
    "grid_ny":                 _P["grid_ny"],
    "domain_lx":               _P["domain_lx"],
    "domain_ly":               _P["domain_ly"],
    "dt":                      _P["dt"],
    "steps":                   _P["nsteps"],
    "reynolds_number":         _P["reynolds_number"],
    "velocity_shear":          _P["velocity_shear"],
    "perturbation_amplitude":  _P["perturbation_amplitude"],
    "perturbation_mode":       _P["perturbation_mode"],
}

EXPECTED_FIELD_LENGTH = REFERENCE_PARAMS["grid_nx"] * REFERENCE_PARAMS["grid_ny"]  # 2048

DIV_RMS_LIMIT: float = 1e-10      # incompressibility acceptance criterion


def reference_tolerance(key: str, backend: str) -> float:
    """A2 band of V1 for one diagnostic: max(abs * f, rel * f * |reference|),
    f = the backend's scalar_factor override (Higham s4.2 hybrid rule, the
    same rule as test_golden_vectors_cross_backend.py)."""
    tols = REFERENCE_VECTOR["tolerances"]
    f = tols.get("backend_overrides", {}).get(backend, {}).get("scalar_factor", 1.0)
    if key == "divergence_rms":
        a, r = tols["scalar_div_rms_abs"], tols["scalar_div_rms_rel"]
    else:
        a, r = tols["scalar_abs"], tols["scalar_rel"]
    return max(a * f, r * f * abs(REFERENCE_DIAGNOSTICS[key]))

# ── Fixtures ──────────────────────────────────────────────────────────────────

@pytest.fixture(scope="session")
def backend_urls() -> dict[str, str]:
    """Map of backend name -> base URL for all 5 backends."""
    return BACKEND_URLS


@pytest.fixture(scope="session")
def log_url() -> str:
    """Base URL for the kh-log-service."""
    return LOG_SERVICE_URL


@pytest.fixture(scope="session")
def test_session_id() -> str:
    """Unique session ID shared across all integration tests in this run."""
    return f"pytest-{uuid.uuid4().hex[:12]}"


@pytest.fixture(scope="session")
def simulate_responses(backend_urls: dict[str, str]) -> dict[str, requests.Response]:
    """
    Run POST /simulate on every backend once per pytest session.
    Cached so physics validation tests share the same HTTP calls.
    Tests that depend on this fixture will be skipped per-backend if the
    backend is unreachable.
    """
    results: dict[str, requests.Response] = {}
    for name, base_url in backend_urls.items():
        try:
            resp = requests.post(
                f"{base_url}/simulate",
                json=REFERENCE_PARAMS,
                timeout=30,
            )
            results[name] = resp
        except requests.exceptions.ConnectionError:
            results[name] = None  # type: ignore[assignment]
    return results
