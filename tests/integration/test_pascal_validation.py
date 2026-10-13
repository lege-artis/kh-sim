# test_pascal_validation.py -- Pascal backend: exact request-validation contract
# Task: KH-018 extension (2026-09-26)
# Finding: kh-sim/_audit/FINDING-PASCAL-HTTP-PARAMETER-DROP-2026-09-16.md (v0.4, B-3/B-4)
#
# Pascal-specific, unlike test_conformance.py: these pin exactly what THIS backend
# refuses and how, because its radix-2 FFT makes a non-power-of-two grid a request
# error here and a valid request elsewhere.
#
# Each case is a lock for one fix. Reverting the fix turns its case red; that was
# checked by mutation on 2026-09-26 before this file shipped. In particular the
# HTTP power-of-two check needs its own lock: without it the FFT1D guard still
# stops the grid, but as a 500 - so every backend-agnostic test stays green.
#
# A SKIP here does not mean "fine". If an earlier module in the same session
# damaged the backend (the pre-fix build wedged on grid_nx=48), this module finds
# it unreachable and skips; test_conformance.py reports that damage as a failure,
# and the APPLY.ps1 gate requires 0 skipped for exactly this reason.
#
# Inputs that are computable but physically questionable (dt < 0, steps < 0,
# Re < 0, negative domain) are deliberately absent: accepting them is current
# behaviour, and whether to refuse them is an open contract decision.

import json

import pytest
import requests

from conftest import BACKEND_URLS

URL = BACKEND_URLS["kh-pascal"]
BASE = {"grid_nx": 32, "grid_ny": 16, "steps": 5}

# (case id, request body, expected status, expected "field" for a 400)
CASES = [
    ("grid_nx_48_non_pow2",   dict(BASE, grid_nx=48),            400, "grid_nx"),
    ("grid_ny_24_non_pow2",   dict(BASE, grid_ny=24),            400, "grid_ny"),
    ("grid_nx_zero",          dict(BASE, grid_nx=0),             400, "grid_nx"),
    ("grid_nx_negative",      dict(BASE, grid_nx=-8),            400, "grid_nx"),
    ("grid_nx_above_axis_cap", dict(BASE, grid_nx=131072),       400, "grid_nx"),
    ("domain_lx_zero",        dict(BASE, domain_lx=0.0),         400, "domain_lx"),
    ("domain_ly_zero",        dict(BASE, domain_ly=0.0),         400, "domain_ly"),
    ("reynolds_number_zero",  dict(BASE, reynolds_number=0.0),   400, "reynolds_number"),
    ("initial_omega_short",   dict(BASE, initial_omega=[0.0] * 10), 400, "initial_omega"),
    ("init_omega_short",      dict(BASE, init_omega=[0.0] * 10),    400, "initial_omega"),
    ("both_ic_keys",          dict(BASE, init_omega=[0.0] * 512, initial_omega=[0.0] * 512), 400, "initial_omega"),
    ("band_forcing_lx_ne_ly", dict(BASE, forcing_type="band", forcing_amplitude=1.0,
                                   forcing_kmin=1, forcing_kmax=3), 400, "forcing_type"),
    ("unknown_forcing_type",  dict(BASE, forcing_type="bogus", forcing_amplitude=1.0), 400, "forcing_type"),
    ("body_not_json",         "this is not json",                400, "body"),
    ("body_json_array",       [1, 2, 3],                         400, "body"),
    # A computation that fails is the backend's error, not the caller's.
    ("dt_beyond_cfl",         dict(BASE, dt=1.0),                500, None),
]


@pytest.fixture(scope="module", autouse=True)
def _pascal_reachable():
    try:
        requests.get(f"{URL}/health", timeout=5)
    except requests.exceptions.RequestException:
        pytest.skip(f"kh-pascal not reachable at {URL}")


@pytest.mark.parametrize("case_id,body,status,field", CASES, ids=[c[0] for c in CASES])
def test_pascal_validation(case_id, body, status, field) -> None:
    data = body if isinstance(body, str) else json.dumps(body)
    resp = requests.post(f"{URL}/simulate", data=data,
                         headers={"Content-Type": "application/json"}, timeout=30)
    assert resp.status_code == status, (
        f"{case_id}: expected {status}, got {resp.status_code}: {resp.text[:200]}")
    payload = resp.json()
    assert isinstance(payload.get("error"), str) and payload["error"], (
        f"{case_id}: no 'error' message: {payload}")
    if status == 400:
        assert payload.get("field") == field, (
            f"{case_id}: field {payload.get('field')!r}, expected {field!r}: {payload}")
    else:
        assert payload.get("backend") == "kh-pascal", (
            f"{case_id}: a 500 must name the backend: {payload}")


def test_pascal_info_declares_grid_constraint() -> None:
    """A caller can discover the radix-2 limit before sending a grid."""
    info = requests.get(f"{URL}/info", timeout=5).json()
    assert info.get("grid_constraint") == "power-of-two", info
    assert isinstance(info.get("build_date"), str) and len(info["build_date"]) == 10, info
