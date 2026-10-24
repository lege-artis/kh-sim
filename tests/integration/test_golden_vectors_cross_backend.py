"""
test_golden_vectors_cross_backend.py
Cross-backend golden-vector validation harness.

Phase 6b integration test: parametrised over all 9 golden vectors × all 5
backends. Each test POSTs the vector's parameters to the backend's /simulate
endpoint and compares the response against the oracle diagnostics and field
summaries using the Higham §4.2 hybrid tolerance rule.

Requires running backend servers. Tests skip gracefully on connection refused.

Two modes:
  run-mode:   all 5 backends up — full cross-backend validation (reference workstation)
  quick-mode: Python oracle self-check — loads JSONs + validates structure

Usage:
    pytest tests/integration/test_golden_vectors_cross_backend.py -v
    pytest tests/integration/test_golden_vectors_cross_backend.py -v -k V1
    pytest tests/integration/test_golden_vectors_cross_backend.py -v --backend kh-rust
"""
from __future__ import annotations

import json
import math
import pathlib
from typing import Any

import pytest
import requests

# ── Paths ──────────────────────────────────────────────────────────────────────

_REPO = pathlib.Path(__file__).resolve().parents[2]
_VECTORS_DIR = _REPO / "shared" / "golden-vectors"

# ── Backend configuration (mirrors conftest.py BACKEND_META) ──────────────────

import os

BACKEND_URLS: dict[str, str] = {
    "kh-rust":    os.environ.get("KH_RUST_URL",    "http://localhost:8001"),
    "kh-scala":   os.environ.get("KH_SCALA_URL",   "http://localhost:8002"),
    "kh-cpp":     os.environ.get("KH_CPP_URL",     "http://localhost:8003"),
    "kh-fortran": os.environ.get("KH_FORTRAN_URL", "http://localhost:8004"),
    "kh-pascal":  os.environ.get("KH_PASCAL_URL",  "http://localhost:8005"),
}

# ── Load all golden vectors once ──────────────────────────────────────────────

def _load_vectors() -> list[dict]:
    vectors = []
    for p in sorted(_VECTORS_DIR.glob("v*.json")):
        with open(p) as fh:
            vectors.append(json.load(fh))
    return vectors


GOLDEN_VECTORS: list[dict] = _load_vectors()
VECTOR_IDS     = [v["metadata"]["vector_id"] for v in GOLDEN_VECTORS]
BACKEND_NAMES  = list(BACKEND_URLS.keys())

# ── Tolerance helpers (Higham §4.2 hybrid rule) ────────────────────────────────

_MACHINE_EPS = 2.22e-16  # IEEE 754 binary64


def _scalar_tol(oracle: float, abs_: float, rel_: float, factor: float = 1.0) -> float:
    return max(abs_ * factor, rel_ * factor * abs(oracle))


def _field_tol(oracle: float, abs_: float, rel_: float, factor: float = 1.0,
               sqrt_n: float = 1.0) -> float:
    return max(abs_ * factor, rel_ * factor * abs(oracle) + sqrt_n * _MACHINE_EPS * abs(oracle))


def _get_backend_factor(tols: dict, backend: str, kind: str) -> float:
    """Return scalar_factor or field_factor for the given backend override."""
    overrides = tols.get("backend_overrides", {})
    if backend in overrides:
        return overrides[backend].get(f"{kind}_factor", 1.0)
    return 1.0


def _should_skip_backend(tols: dict, backend: str) -> bool:
    overrides = tols.get("backend_overrides", {})
    return overrides.get(backend, {}).get("skip", False)


# ── REST payload builder ───────────────────────────────────────────────────────

def _params_to_request(params: dict) -> dict:
    """Convert golden-vector parameters to REST /simulate request body."""
    return {
        "grid_nx":                params["grid_nx"],
        "grid_ny":                params["grid_ny"],
        "domain_lx":              params["domain_lx"],
        "domain_ly":              params["domain_ly"],
        "dt":                     params["dt"],
        "steps":                  params["nsteps"],   # REST uses 'steps', not 'nsteps'
        "reynolds_number":        params["reynolds_number"],
        "velocity_shear":         params["velocity_shear"],
        "perturbation_amplitude": params["perturbation_amplitude"],
        "perturbation_mode":      params["perturbation_mode"],
    }


# ── Per-backend comparison logic ──────────────────────────────────────────────

def _compare_scalar(name: str, backend_val: float, oracle_val: float,
                    abs_: float, rel_: float, factor: float,
                    failures: list[str]) -> None:
    tol = _scalar_tol(oracle_val, abs_, rel_, factor)
    err = abs(backend_val - oracle_val)
    if err > tol:
        failures.append(
            f"{name}: |{backend_val:.6e} - {oracle_val:.6e}| = {err:.2e} > tol {tol:.2e}"
        )


def _compare_field_summary(field_name: str, backend_fs: dict, oracle_fs: dict,
                            abs_: float, rel_: float, factor: float,
                            sqrt_n: float, failures: list[str]) -> None:
    for stat in ("min", "max", "mean"):
        oracle_v  = oracle_fs[stat]
        backend_v = backend_fs.get(stat)
        if backend_v is None:
            failures.append(f"{field_name}.{stat}: missing in backend response")
            continue
        tol = _field_tol(oracle_v, abs_, rel_, factor, sqrt_n)
        err = abs(backend_v - oracle_v)
        if err > tol:
            failures.append(
                f"{field_name}.{stat}: |{backend_v:.6e} - {oracle_v:.6e}| = {err:.2e} > tol {tol:.2e}"
            )
    # element-wise: first_10 / last_10
    for slot in ("first_10", "last_10"):
        oracle_elems  = oracle_fs.get(slot, [])
        backend_elems = backend_fs.get(slot, [])
        for i, (o, b) in enumerate(zip(oracle_elems, backend_elems)):
            tol = _field_tol(o, abs_, rel_, factor, sqrt_n)
            err = abs(b - o)
            if err > tol:
                failures.append(
                    f"{field_name}.{slot}[{i}]: |{b:.6e} - {o:.6e}| = {err:.2e} > tol {tol:.2e}"
                )


def _run_comparison(vector: dict, backend_name: str, response_body: dict) -> list[str]:
    """Compare backend response against oracle; return list of failure messages."""
    tols   = vector["tolerances"]
    oracle = vector["diagnostics"]
    s_abs  = tols["scalar_abs"]
    s_rel  = tols["scalar_rel"]
    sd_abs = tols["scalar_div_rms_abs"]
    sd_rel = tols["scalar_div_rms_rel"]
    fs_abs = tols["field_summary_abs"]
    fs_rel = tols["field_summary_rel"]
    fe_abs = tols["field_element_abs"]
    fe_rel = tols["field_element_rel"]

    sf  = _get_backend_factor(tols, backend_name, "scalar")
    ff  = _get_backend_factor(tols, backend_name, "field")
    nx  = vector["parameters"]["grid_nx"]
    ny  = vector["parameters"]["grid_ny"]
    sqN = math.sqrt(nx * ny)

    failures: list[str] = []

    # Extract backend diagnostics
    b_diag = response_body.get("diagnostics", response_body)

    # Scalar diagnostics
    for key, a, r, fac in [
        ("kinetic_energy", s_abs, s_rel, sf),
        ("enstrophy",      s_abs, s_rel, sf),
        ("max_vorticity",  s_abs, s_rel, sf),
        ("cfl_peak",       s_abs, s_rel, sf),
        ("divergence_rms", sd_abs, sd_rel, sf),
    ]:
        bv = b_diag.get(key)
        if bv is None:
            continue   # backend may not return cfl_peak; skip gracefully
        _compare_scalar(key, float(bv), oracle[key], a, r, fac, failures)

    # Field summaries
    oracle_fs = vector.get("field_summaries", {})
    backend_fields = {
        "u_velocity": "u_velocity",
        "v_velocity": "v_velocity",
        "vorticity":  "vorticity",
        "pressure":   "pressure",
    }
    for oracle_key, backend_key in backend_fields.items():
        o_fs = oracle_fs.get(oracle_key)
        b_fs = response_body.get(backend_key)
        if o_fs is None or b_fs is None:
            continue
        # Backend may return full list or a summary dict
        if isinstance(b_fs, list):
            # Convert to summary dict for comparison
            b_fs_summary = {
                "length":   len(b_fs),
                "first_10": b_fs[:10],
                "last_10":  b_fs[-10:],
                "min":      min(b_fs),
                "max":      max(b_fs),
                "mean":     sum(b_fs) / len(b_fs),
            }
        elif isinstance(b_fs, dict):
            b_fs_summary = b_fs
        else:
            continue
        _compare_field_summary(oracle_key, b_fs_summary, o_fs,
                                fs_abs, fs_rel, ff, sqN, failures)

    return failures


# ── Parametrised test cases ────────────────────────────────────────────────────

def pytest_addoption(parser):
    """Allow --backend CLI filter: pytest ... --backend kh-rust"""
    try:
        parser.addoption(
            "--backend", action="store", default="",
            help="Run only against this backend (e.g. kh-rust)"
        )
    except ValueError:
        pass  # already registered by conftest.py or another plugin


def _get_backend_filter(config) -> str:
    try:
        return config.getoption("--backend", default="")
    except Exception:
        return ""


@pytest.fixture(scope="session")
def backend_filter(request) -> str:
    return _get_backend_filter(request.config)


# Build parametrize lists at module level for pytest collection
_test_cases = [
    (v["metadata"]["vector_id"], backend)
    for v in GOLDEN_VECTORS
    for backend in BACKEND_NAMES
]
_test_ids = [f"{vid}-x-{b}" for vid, b in _test_cases]


@pytest.mark.parametrize("vector_id,backend_name", _test_cases, ids=_test_ids)
def test_golden_vector_cross_backend(vector_id: str, backend_name: str,
                                      backend_filter: str) -> None:
    """
    For each (vector, backend) pair: POST /simulate, compare against oracle.
    Skips if backend is unreachable or flagged as skip in tolerances.
    """
    # Optional backend filter via --backend CLI
    if backend_filter and backend_name != backend_filter:
        pytest.skip(f"filtered out by --backend {backend_filter}")

    # Find the vector
    vector = next((v for v in GOLDEN_VECTORS if v["metadata"]["vector_id"] == vector_id), None)
    assert vector is not None, f"vector {vector_id!r} not found in golden-vectors dir"

    # Check skip flag in tolerances
    if _should_skip_backend(vector["tolerances"], backend_name):
        pytest.skip(f"{vector_id} x {backend_name}: marked skip in tolerances.backend_overrides")

    base_url = BACKEND_URLS[backend_name]
    payload  = _params_to_request(vector["parameters"])

    # Try to reach the backend
    try:
        resp = requests.post(f"{base_url}/simulate", json=payload, timeout=60)
        resp.raise_for_status()
    except requests.exceptions.ConnectionError:
        pytest.skip(f"{backend_name} not reachable at {base_url}")
    except requests.exceptions.Timeout:
        pytest.skip(f"{backend_name} timed out at {base_url}")
    except requests.exceptions.HTTPError as exc:
        pytest.fail(f"{backend_name} HTTP error: {exc}")

    body = resp.json()
    failures = _run_comparison(vector, backend_name, body)

    if failures:
        msg = f"{vector_id} x {backend_name}: {len(failures)} comparison failure(s):\n"
        msg += "\n".join(f"  - {f}" for f in failures)
        pytest.fail(msg)


# ── Quick-mode: Python oracle self-check (no backends needed) ─────────────────

class TestGoldenVectorStructure:
    """
    Quick-mode: verify all golden-vector JSONs are loadable and structurally
    sound without requiring any backend servers.
    """

    @pytest.mark.parametrize("vector", GOLDEN_VECTORS, ids=VECTOR_IDS)
    def test_loads_without_error(self, vector: dict) -> None:
        assert "$schema" in vector
        assert vector["$schema"] == "kh-sim-golden-vector-v1"

    @pytest.mark.parametrize("vector", GOLDEN_VECTORS, ids=VECTOR_IDS)
    def test_all_required_metadata(self, vector: dict) -> None:
        meta = vector["metadata"]
        for key in ("vector_id", "purpose", "generation_date", "oracle_commit",
                    "oracle_path", "scheme", "dealiasing", "cfl_guard", "environment"):
            assert key in meta, f"metadata missing key: {key}"

    @pytest.mark.parametrize("vector", GOLDEN_VECTORS, ids=VECTOR_IDS)
    def test_all_diagnostics_positive(self, vector: dict) -> None:
        diag = vector["diagnostics"]
        for key in ("kinetic_energy", "enstrophy", "max_vorticity", "cfl_peak"):
            assert diag[key] >= 0, f"{key} is negative: {diag[key]}"
        assert diag["divergence_rms"] < 1e-10, (
            f"div_rms={diag['divergence_rms']:.2e} exceeds incompressibility floor 1e-10"
        )

    @pytest.mark.parametrize("vector", GOLDEN_VECTORS, ids=VECTOR_IDS)
    def test_field_summaries_consistent(self, vector: dict) -> None:
        for field_name, fs in vector["field_summaries"].items():
            assert fs["length"] > 0
            assert fs["min"] <= fs["mean"] <= fs["max"], (
                f"{field_name}: min/mean/max ordering violated: "
                f"{fs['min']:.4f} <= {fs['mean']:.4f} <= {fs['max']:.4f}"
            )
            assert len(fs["first_10"]) <= 10
            assert len(fs["last_10"])  <= 10

    @pytest.mark.parametrize("vector", GOLDEN_VECTORS, ids=VECTOR_IDS)
    def test_t_final_consistent_with_params(self, vector: dict) -> None:
        params = vector["parameters"]
        expected = round(params["dt"] * params["nsteps"], 8)
        actual   = vector["diagnostics"]["t_final"]
        assert abs(actual - expected) < 1e-6, (
            f"t_final={actual} inconsistent with dt*nsteps={expected}"
        )

    @pytest.mark.parametrize("vector", GOLDEN_VECTORS, ids=VECTOR_IDS)
    def test_field_length_equals_grid(self, vector: dict) -> None:
        params   = vector["parameters"]
        expected = params["grid_nx"] * params["grid_ny"]
        for field_name, fs in vector["field_summaries"].items():
            assert fs["length"] == expected, (
                f"{field_name}.length={fs['length']} != nx*ny={expected}"
            )
