"""
test_golden_vectors_python.py — golden-vector regression tests, Python reference backend.

Loads each shared/golden-vectors/v*.json, calls simulate() in-process via Mode 2
(direct kernel invocation, no HTTP), compares the result against oracle diagnostics
and field_summaries using the Higham §4.2 hybrid tolerance rule.

Pass criterion per vector (all 27 cases = 9 vectors × 3 functions):
  - test_scalar_diagnostics : KE, enstrophy, max_vort, div_rms, cfl_peak all within
                              (abs|rel) hybrid threshold × backend_factor
  - test_field_summaries    : min/max/mean per field (u, v, ω, ψ) within field_summary
                              threshold
  - test_field_elements     : first_10 / last_10 per field within field_element
                              threshold (includes sqrt(N)·ε term when sqrt_n_scaling=true)

python-numpy has no backend_overrides entry in the JSON → all factors default to 1.0.
Since the oracle compares to itself all diffs should be ≤ 2e-15.

Usage (from kh-sim root with .venv active):
    pytest tests/integration/test_golden_vectors_python.py -v
"""
from __future__ import annotations

import functools
import json
import math
import sys
from pathlib import Path
from typing import Any

import numpy as np
import pytest

# ── Repo path bootstrap (allows running without editable install) ──────────────
_REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(_REPO))

from shared.physics.kh_physics import simulate, SimulationRequest, SimulationResult

# ── Constants ─────────────────────────────────────────────────────────────────
_VECTORS_DIR = _REPO / "shared" / "golden-vectors"
_BACKEND     = "python-numpy"
_MACHINE_EPS = 2.22e-16  # IEEE 754 binary64


# ── Vector loading ─────────────────────────────────────────────────────────────

def _load_vectors() -> list[tuple[Path, dict]]:
    return [
        (p, json.loads(p.read_text()))
        for p in sorted(_VECTORS_DIR.glob("v*.json"))
    ]


_VECTORS    = _load_vectors()
_VECTOR_IDS = [v["metadata"]["vector_id"] for _, v in _VECTORS]


# ── Simulation (cached to avoid running 3× per vector) ────────────────────────

@functools.lru_cache(maxsize=None)
def _simulate_cached(vector_id: str) -> SimulationResult:
    _, vector = next(
        (p, v) for p, v in _VECTORS
        if v["metadata"]["vector_id"] == vector_id
    )
    params = vector["parameters"]
    req = SimulationRequest(
        grid_nx               = params["grid_nx"],
        grid_ny               = params["grid_ny"],
        domain_lx             = params["domain_lx"],
        domain_ly             = params["domain_ly"],
        dt                    = params["dt"],
        steps                 = params["nsteps"],
        reynolds_number       = params["reynolds_number"],
        velocity_shear        = params["velocity_shear"],
        perturbation_amplitude= params["perturbation_amplitude"],
        perturbation_mode     = params["perturbation_mode"],
    )
    return simulate(req)


# ── Tolerance helpers (Higham §4.2 hybrid) ────────────────────────────────────

def _backend_factor(tols: dict, kind: str) -> float:
    """Return scalar_factor or field_factor for python-numpy (defaults to 1.0)."""
    return tols.get("backend_overrides", {}).get(_BACKEND, {}).get(f"{kind}_factor", 1.0)


def _scalar_tol(oracle: float, abs_: float, rel_: float, factor: float) -> float:
    return max(abs_ * factor, rel_ * factor * abs(oracle))


def _field_tol(oracle: float, abs_: float, rel_: float,
               factor: float, sqrt_n: float, use_sqrt_n: bool) -> float:
    eps_term = sqrt_n * _MACHINE_EPS * abs(oracle) if use_sqrt_n else 0.0
    return max(abs_ * factor, rel_ * factor * abs(oracle) + eps_term)


# ── Test cases (9 vectors × 3 functions = 27 total) ───────────────────────────

@pytest.mark.parametrize("_path,vector", _VECTORS, ids=_VECTOR_IDS)
def test_scalar_diagnostics(_path: Path, vector: dict) -> None:
    """KE, enstrophy, max_vorticity, cfl_peak, divergence_rms all within threshold."""
    vid    = vector["metadata"]["vector_id"]
    result = _simulate_cached(vid)
    diag   = result.diagnostics
    tols   = vector["tolerances"]
    oracle = vector["diagnostics"]
    sf     = _backend_factor(tols, "scalar")

    checks = [
        ("kinetic_energy",  tols["scalar_abs"],         tols["scalar_rel"]),
        ("enstrophy",       tols["scalar_abs"],         tols["scalar_rel"]),
        ("max_vorticity",   tols["scalar_abs"],         tols["scalar_rel"]),
        ("cfl_peak",        tols["scalar_abs"],         tols["scalar_rel"]),
        ("divergence_rms",  tols["scalar_div_rms_abs"], tols["scalar_div_rms_rel"]),
    ]
    failures: list[str] = []
    for key, abs_, rel_ in checks:
        actual = getattr(diag, key)
        ov     = oracle[key]
        tol    = _scalar_tol(ov, abs_, rel_, sf)
        diff   = abs(actual - ov)
        if diff >= tol:
            failures.append(
                f"{key}: |diff|={diff:.3e} >= tol={tol:.3e} "
                f"(actual={actual:.6e}, oracle={ov:.6e})"
            )
    if failures:
        pytest.fail(
            f"{_path.name}: {len(failures)} scalar failure(s):\n" +
            "\n".join(f"  {f}" for f in failures)
        )


@pytest.mark.parametrize("_path,vector", _VECTORS, ids=_VECTOR_IDS)
def test_field_summaries(_path: Path, vector: dict) -> None:
    """min/max/mean of u, v, vorticity, pressure within field_summary threshold."""
    vid      = vector["metadata"]["vector_id"]
    result   = _simulate_cached(vid)
    tols     = vector["tolerances"]
    oracle   = vector["field_summaries"]
    ff       = _backend_factor(tols, "field")
    nx, ny   = vector["parameters"]["grid_nx"], vector["parameters"]["grid_ny"]
    sqrt_n   = math.sqrt(nx * ny)
    use_sqn  = tols.get("sqrt_n_scaling", False)

    field_map = {
        "u_velocity": result.u_velocity,
        "v_velocity": result.v_velocity,
        "vorticity":  result.vorticity,
        "pressure":   result.pressure,
    }
    failures: list[str] = []
    for fname, flat in field_map.items():
        arr  = np.array(flat)
        o_fs = oracle[fname]
        for stat, bv, ov in [
            ("min",  float(np.min(arr)),  o_fs["min"]),
            ("max",  float(np.max(arr)),  o_fs["max"]),
            ("mean", float(np.mean(arr)), o_fs["mean"]),
        ]:
            tol  = _field_tol(ov, tols["field_summary_abs"], tols["field_summary_rel"],
                              ff, sqrt_n, use_sqn)
            diff = abs(bv - ov)
            if diff >= tol:
                failures.append(
                    f"{fname}.{stat}: |diff|={diff:.3e} >= tol={tol:.3e} "
                    f"(actual={bv:.6e}, oracle={ov:.6e})"
                )
    if failures:
        pytest.fail(
            f"{_path.name}: {len(failures)} field summary failure(s):\n" +
            "\n".join(f"  {f}" for f in failures)
        )


@pytest.mark.parametrize("_path,vector", _VECTORS, ids=_VECTOR_IDS)
def test_field_elements(_path: Path, vector: dict) -> None:
    """first_10 and last_10 elements of u, v, vorticity, pressure within threshold."""
    vid      = vector["metadata"]["vector_id"]
    result   = _simulate_cached(vid)
    tols     = vector["tolerances"]
    oracle   = vector["field_summaries"]
    ff       = _backend_factor(tols, "field")
    nx, ny   = vector["parameters"]["grid_nx"], vector["parameters"]["grid_ny"]
    sqrt_n   = math.sqrt(nx * ny)
    use_sqn  = tols.get("sqrt_n_scaling", False)

    field_map = {
        "u_velocity": result.u_velocity,
        "v_velocity": result.v_velocity,
        "vorticity":  result.vorticity,
        "pressure":   result.pressure,
    }
    failures: list[str] = []
    for fname, flat in field_map.items():
        o_fs = oracle[fname]
        for slot, backend_elems, oracle_elems in [
            ("first_10", flat[:10],  o_fs["first_10"]),
            ("last_10",  flat[-10:], o_fs["last_10"]),
        ]:
            for i, (bv, ov) in enumerate(zip(backend_elems, oracle_elems)):
                tol  = _field_tol(ov, tols["field_element_abs"], tols["field_element_rel"],
                                  ff, sqrt_n, use_sqn)
                diff = abs(bv - ov)
                if diff >= tol:
                    failures.append(
                        f"{fname}.{slot}[{i}]: |diff|={diff:.3e} >= tol={tol:.3e} "
                        f"(actual={bv:.6e}, oracle={ov:.6e})"
                    )
    if failures:
        pytest.fail(
            f"{_path.name}: {len(failures)} field element failure(s):\n" +
            "\n".join(f"  {f}" for f in failures)
        )
