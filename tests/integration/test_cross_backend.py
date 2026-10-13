# test_cross_backend.py — schema parity and consistency across all 5 backends
# Task: KH-018
#
# These tests verify that all backends present a consistent API contract —
# they do NOT compare physics output values across backends (that is KH-VAL).
# Checks:
#   - /simulate top-level JSON keys are identical for all backends
#   - /simulate diagnostics keys are identical for all backends
#   - Vector field lengths are identical for same grid params
#   - t_final is identical (same timestep and step count)
#   - compute_time_ms is present and positive across all backends

import pytest
import requests

from conftest import (
    BACKEND_NAMES,
    REFERENCE_PARAMS,
    EXPECTED_FIELD_LENGTH,
)

VECTOR_FIELDS = ("u_velocity", "v_velocity", "vorticity", "pressure")


def _bodies(sim: dict) -> dict[str, dict]:
    """
    Return {backend_name: response_body} for all reachable backends.
    Skips the test suite if fewer than 2 backends are reachable.
    """
    result = {}
    for name in BACKEND_NAMES:
        resp = sim.get(name)
        if resp is not None and resp.status_code == 200:
            result[name] = resp.json()
    if len(result) < 2:
        pytest.skip(
            f"Cross-backend tests require >= 2 reachable backends; "
            f"got {len(result)}: {list(result.keys())}"
        )
    return result


def test_all_backends_same_top_level_keys(simulate_responses: dict) -> None:
    """All reachable backends return the same set of top-level JSON keys from /simulate."""
    bodies = _bodies(simulate_responses)
    key_sets = {name: frozenset(body.keys()) for name, body in bodies.items()}
    reference_name = next(iter(key_sets))
    reference_keys = key_sets[reference_name]
    for name, keys in key_sets.items():
        assert keys == reference_keys, (
            f"Top-level key mismatch between {reference_name} and {name}: "
            f"extra={keys - reference_keys}, missing={reference_keys - keys}"
        )


def test_all_backends_same_diagnostics_keys(simulate_responses: dict) -> None:
    """All reachable backends return the same diagnostics keys."""
    bodies = _bodies(simulate_responses)
    diag_key_sets = {name: frozenset(body["diagnostics"].keys()) for name, body in bodies.items()}
    reference_name = next(iter(diag_key_sets))
    reference_keys = diag_key_sets[reference_name]
    for name, keys in diag_key_sets.items():
        assert keys == reference_keys, (
            f"Diagnostics key mismatch between {reference_name} and {name}: "
            f"extra={keys - reference_keys}, missing={reference_keys - keys}"
        )


@pytest.mark.parametrize("field_name", VECTOR_FIELDS)
def test_all_backends_same_field_length(field_name: str, simulate_responses: dict) -> None:
    """All reachable backends return the same vector field length for the same grid."""
    bodies = _bodies(simulate_responses)
    for name, body in bodies.items():
        actual = len(body.get(field_name, []))
        assert actual == EXPECTED_FIELD_LENGTH, (
            f"{name}: {field_name} length {actual} != expected {EXPECTED_FIELD_LENGTH} "
            f"(grid {REFERENCE_PARAMS['grid_nx']}x{REFERENCE_PARAMS['grid_ny']})"
        )


def test_all_backends_same_t_final(simulate_responses: dict) -> None:
    """All backends report the same t_final for the same steps and dt."""
    bodies = _bodies(simulate_responses)
    expected_t = REFERENCE_PARAMS["steps"] * REFERENCE_PARAMS["dt"]
    for name, body in bodies.items():
        actual = body.get("t_final", -1.0)
        assert abs(actual - expected_t) < 1e-9, (
            f"{name}: t_final {actual} != expected {expected_t:.4f}"
        )


def test_all_backends_same_steps_completed(simulate_responses: dict) -> None:
    """All backends report the same steps_completed for the same request."""
    bodies = _bodies(simulate_responses)
    expected = REFERENCE_PARAMS["steps"]
    for name, body in bodies.items():
        actual = body.get("steps_completed")
        assert actual == expected, (
            f"{name}: steps_completed {actual} != expected {expected}"
        )


def test_all_backends_compute_time_positive(simulate_responses: dict) -> None:
    """compute_time_ms is positive for all reachable backends."""
    bodies = _bodies(simulate_responses)
    for name, body in bodies.items():
        ct = body.get("compute_time_ms", -1)
        assert isinstance(ct, (int, float)) and ct > 0, (
            f"{name}: compute_time_ms={ct} is not positive"
        )


def test_all_backends_grid_dimensions_consistent(simulate_responses: dict) -> None:
    """All backends echo back consistent grid_nx and grid_ny."""
    bodies = _bodies(simulate_responses)
    for name, body in bodies.items():
        assert body.get("grid_nx") == REFERENCE_PARAMS["grid_nx"] and \
               body.get("grid_ny") == REFERENCE_PARAMS["grid_ny"], (
            f"{name}: grid dimensions {body.get('grid_nx')}x{body.get('grid_ny')} "
            f"!= requested {REFERENCE_PARAMS['grid_nx']}x{REFERENCE_PARAMS['grid_ny']}"
        )


def test_all_backends_diagnostics_are_numeric(simulate_responses: dict) -> None:
    """All diagnostic values are finite floats for all reachable backends."""
    import math
    bodies = _bodies(simulate_responses)
    diag_fields = ("kinetic_energy", "enstrophy", "max_vorticity", "divergence_rms")
    for name, body in bodies.items():
        diag = body.get("diagnostics", {})
        for field in diag_fields:
            val = diag.get(field)
            assert isinstance(val, (int, float)), (
                f"{name}: diagnostics.{field} is not numeric: {val!r}"
            )
            assert math.isfinite(val), (
                f"{name}: diagnostics.{field} is non-finite: {val}"
            )
