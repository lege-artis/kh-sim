# test_simulate.py — POST /simulate physics validation for all 5 backends
# Task: KH-018  |  Spec: kh-sim/shared/physics/KH-PHYSICS.md section 7
#
# Canonical parameters: 64x32 grid, 100 steps, dt=0.001, Re=1000 (read from V1)
# Reference:            kh-sim/shared/golden-vectors/v1-canonical-reference.json
#                       (replaced the retired kh_reference_output.json, 2026-10-08)
# Acceptance criteria (modelling rules s4, class A2):
#   kinetic_energy  : within the A2 band of V1 (reference_tolerance)
#   enstrophy       : within the A2 band of V1
#   max_vorticity   : within the A2 band of V1
#   divergence_rms  : < DIV_RMS_LIMIT (incompressibility gate)

import pytest
import requests

from conftest import (
    BACKEND_NAMES,
    REFERENCE_PARAMS,
    REFERENCE_DIAGNOSTICS,
    EXPECTED_FIELD_LENGTH,
    DIV_RMS_LIMIT,
    reference_tolerance,
)

REQUIRED_SIM_FIELDS = {
    "backend", "language", "steps_completed", "t_final",
    "grid_nx", "grid_ny", "u_velocity", "v_velocity",
    "vorticity", "pressure", "diagnostics", "compute_time_ms",
}
REQUIRED_DIAG_FIELDS = {
    "kinetic_energy", "enstrophy", "max_vorticity", "divergence_rms",
}
VECTOR_FIELDS = ("u_velocity", "v_velocity", "vorticity", "pressure")


# ── Shared session fixture: run /simulate once per backend per session ────────

@pytest.fixture(scope="session")
def sim(simulate_responses: dict) -> dict:
    """Return simulate_responses — skip per-backend if response is None (unreachable)."""
    return simulate_responses


def _get_body(sim: dict, backend_name: str) -> dict:
    """Return parsed response body, skipping test if backend unreachable."""
    resp = sim.get(backend_name)
    if resp is None:
        pytest.skip(f"{backend_name} was not reachable during session fixture")
    assert resp.status_code == 200, (
        f"{backend_name}: POST /simulate returned {resp.status_code}: {resp.text[:200]}"
    )
    return resp.json()


# ── HTTP layer ────────────────────────────────────────────────────────────────

@pytest.mark.parametrize("backend_name", BACKEND_NAMES)
def test_simulate_status_200(backend_name: str, sim: dict) -> None:
    """POST /simulate returns HTTP 200."""
    resp = sim.get(backend_name)
    if resp is None:
        pytest.skip(f"{backend_name} not reachable")
    assert resp.status_code == 200, (
        f"{backend_name}: /simulate returned {resp.status_code}: {resp.text[:200]}"
    )


# ── Response schema ───────────────────────────────────────────────────────────

@pytest.mark.parametrize("backend_name", BACKEND_NAMES)
def test_simulate_top_level_fields(backend_name: str, sim: dict) -> None:
    """POST /simulate response contains all required top-level fields."""
    body = _get_body(sim, backend_name)
    missing = REQUIRED_SIM_FIELDS - set(body.keys())
    assert not missing, f"{backend_name}: /simulate response missing: {missing}"


@pytest.mark.parametrize("backend_name", BACKEND_NAMES)
def test_simulate_diagnostics_fields(backend_name: str, sim: dict) -> None:
    """POST /simulate diagnostics block contains all required fields."""
    body = _get_body(sim, backend_name)
    diag = body.get("diagnostics", {})
    missing = REQUIRED_DIAG_FIELDS - set(diag.keys())
    assert not missing, f"{backend_name}: diagnostics missing: {missing}"


@pytest.mark.parametrize("backend_name", BACKEND_NAMES)
def test_simulate_grid_dimensions_echo(backend_name: str, sim: dict) -> None:
    """Response echoes back the requested grid dimensions."""
    body = _get_body(sim, backend_name)
    assert body["grid_nx"] == REFERENCE_PARAMS["grid_nx"], (
        f"{backend_name}: grid_nx {body['grid_nx']} != requested {REFERENCE_PARAMS['grid_nx']}"
    )
    assert body["grid_ny"] == REFERENCE_PARAMS["grid_ny"], (
        f"{backend_name}: grid_ny {body['grid_ny']} != requested {REFERENCE_PARAMS['grid_ny']}"
    )


@pytest.mark.parametrize("backend_name", BACKEND_NAMES)
def test_simulate_steps_completed(backend_name: str, sim: dict) -> None:
    """steps_completed matches requested step count."""
    body = _get_body(sim, backend_name)
    expected = REFERENCE_PARAMS["steps"]
    actual = body.get("steps_completed")
    assert actual == expected, (
        f"{backend_name}: steps_completed {actual} != requested {expected}"
    )


@pytest.mark.parametrize("backend_name", BACKEND_NAMES)
def test_simulate_t_final(backend_name: str, sim: dict) -> None:
    """t_final is approximately steps * dt = 0.1."""
    body = _get_body(sim, backend_name)
    expected_t = REFERENCE_PARAMS["steps"] * REFERENCE_PARAMS["dt"]  # 0.1
    actual_t = body.get("t_final", -1.0)
    assert abs(actual_t - expected_t) < 1e-9, (
        f"{backend_name}: t_final {actual_t} != expected {expected_t}"
    )


@pytest.mark.parametrize("backend_name", BACKEND_NAMES)
@pytest.mark.parametrize("field_name", VECTOR_FIELDS)
def test_simulate_field_length(backend_name: str, field_name: str, sim: dict) -> None:
    """Each output field array has length grid_nx * grid_ny = 2048."""
    body = _get_body(sim, backend_name)
    actual_len = len(body.get(field_name, []))
    assert actual_len == EXPECTED_FIELD_LENGTH, (
        f"{backend_name}: {field_name} length {actual_len} != {EXPECTED_FIELD_LENGTH}"
    )


@pytest.mark.parametrize("backend_name", BACKEND_NAMES)
def test_simulate_compute_time_positive(backend_name: str, sim: dict) -> None:
    """compute_time_ms is a positive number."""
    body = _get_body(sim, backend_name)
    ct = body.get("compute_time_ms", -1)
    assert isinstance(ct, (int, float)) and ct > 0, (
        f"{backend_name}: compute_time_ms={ct} is not positive"
    )


# ── Physics diagnostics validation ───────────────────────────────────────────

def _assert_in_band(backend_name: str, body: dict, key: str) -> None:
    ref = REFERENCE_DIAGNOSTICS[key]
    actual = body["diagnostics"][key]
    tol = reference_tolerance(key, backend_name)
    assert abs(actual - ref) <= tol, (
        f"{backend_name}: {key} {actual:.17g} vs V1 {ref:.17g} "
        f"(|diff| {abs(actual - ref):.3e} > A2 band {tol:.3e})"
    )


@pytest.mark.parametrize("backend_name", BACKEND_NAMES)
def test_simulate_kinetic_energy(backend_name: str, sim: dict) -> None:
    """kinetic_energy within the A2 band of golden vector V1."""
    _assert_in_band(backend_name, _get_body(sim, backend_name), "kinetic_energy")


@pytest.mark.parametrize("backend_name", BACKEND_NAMES)
def test_simulate_enstrophy(backend_name: str, sim: dict) -> None:
    """enstrophy within the A2 band of golden vector V1."""
    _assert_in_band(backend_name, _get_body(sim, backend_name), "enstrophy")


@pytest.mark.parametrize("backend_name", BACKEND_NAMES)
def test_simulate_max_vorticity(backend_name: str, sim: dict) -> None:
    """max_vorticity within the A2 band of golden vector V1."""
    _assert_in_band(backend_name, _get_body(sim, backend_name), "max_vorticity")


@pytest.mark.parametrize("backend_name", BACKEND_NAMES)
def test_simulate_divergence_rms(backend_name: str, sim: dict) -> None:
    """divergence_rms < 1e-10 — incompressibility acceptance criterion."""
    body = _get_body(sim, backend_name)
    div = body["diagnostics"]["divergence_rms"]
    assert div < DIV_RMS_LIMIT, (
        f"{backend_name}: divergence_rms {div:.2e} >= incompressibility limit {DIV_RMS_LIMIT:.2e}"
    )
