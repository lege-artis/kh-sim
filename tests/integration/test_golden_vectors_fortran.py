"""test_golden_vectors_fortran.py — per-backend golden-vector regression tests, Fortran backend.

Mode 1: HTTP POST /simulate to running server on port 8004.
Compares backend response against oracle diagnostics + field_summaries using
Higham §4.2 hybrid tolerance with backend_factor = 1.0 (kh-fortran has no
override entry in JSON; defaults to 1.0×).

Run (requires server running on :8004):
    cd kh-sim/backends/fortran
    cmake -B build -G Ninja && cmake --build build
    ./build/kh-sim-fortran.exe &
    sleep 2
    cd ../..
    pytest tests/integration/test_golden_vectors_fortran.py -v
    kill %1

Each test covers all three check types for that vector:
  1. Scalar diagnostics (KE, enstrophy, max_vort, div_rms, cfl_peak)
  2. Field summary statistics (min/max/mean per field u, v, ω, ψ)
  3. Field element comparison (first_10 / last_10 per field)

If the server is unreachable the test is skipped (not failed) unless
KH_GOLDEN_REQUIRE_SERVER=1 is set in the environment.
"""

from __future__ import annotations

import json
import math
import os
from pathlib import Path

import pytest
import requests

# ── Constants ──────────────────────────────────────────────────────────────────

# Honour KH_FORTRAN_URL as conftest.py does (2026-10-03, finding v0.5 F-5b): the
# same hardcoded URL that made the Pascal gate grade the live container on
# 2026-09-26 (F-5) was here too.
BACKEND_URL  = os.environ.get("KH_FORTRAN_URL", "http://localhost:8004").rstrip("/") + "/simulate"
BACKEND_NAME = "kh-fortran"
MACHINE_EPS  = 2.22e-16

GOLDEN_DIR = (
    Path(__file__).parents[2] / "shared" / "golden-vectors"
)

VECTOR_FILES = sorted(GOLDEN_DIR.glob("v*.json"))

# ── Helpers ────────────────────────────────────────────────────────────────────

def load_vector(filename: str) -> dict:
    path = GOLDEN_DIR / filename
    with open(path) as f:
        return json.load(f)


def backend_factor(tols: dict, kind: str) -> float:
    overrides = tols.get("backend_overrides", {}).get(BACKEND_NAME, {})
    return overrides.get(f"{kind}_factor", 1.0)


def scalar_tol(oracle: float, abs_: float, rel_: float, factor: float) -> float:
    return max(abs_ * factor, rel_ * factor * abs(oracle))


def field_tol(oracle: float, abs_: float, rel_: float, factor: float,
              sqrt_n: float, use_sqn: bool) -> float:
    eps_term = sqrt_n * MACHINE_EPS * abs(oracle) if use_sqn else 0.0
    return max(abs_ * factor, rel_ * factor * abs(oracle) + eps_term)


def params_to_request(params: dict) -> dict:
    return {
        "grid_nx":                params["grid_nx"],
        "grid_ny":                params["grid_ny"],
        "domain_lx":              params["domain_lx"],
        "domain_ly":              params["domain_ly"],
        "dt":                     params["dt"],
        "steps":                  params["nsteps"],
        "reynolds_number":        params["reynolds_number"],
        "velocity_shear":         params["velocity_shear"],
        "perturbation_amplitude": params["perturbation_amplitude"],
        "perturbation_mode":      params["perturbation_mode"],
    }


def require_server() -> bool:
    return os.environ.get("KH_GOLDEN_REQUIRE_SERVER", "") == "1"


# ── Core comparison ────────────────────────────────────────────────────────────

def check_scalars(vector: dict, body: dict) -> list[str]:
    tols   = vector["tolerances"]
    oracle = vector["diagnostics"]
    sf     = backend_factor(tols, "scalar")
    s_abs  = tols["scalar_abs"]
    s_rel  = tols["scalar_rel"]
    sd_abs = tols["scalar_div_rms_abs"]
    sd_rel = tols["scalar_div_rms_rel"]

    b_diag   = body["diagnostics"]
    failures = []

    checks = [
        ("kinetic_energy", s_abs,  s_rel),
        ("enstrophy",      s_abs,  s_rel),
        ("max_vorticity",  s_abs,  s_rel),
        ("divergence_rms", sd_abs, sd_rel),
    ]
    for key, abs_, rel_ in checks:
        if key not in b_diag:
            failures.append(f"{key}: missing from response")
            continue
        bv  = float(b_diag[key])
        ov  = float(oracle[key])
        tol = scalar_tol(ov, abs_, rel_, sf)
        err = abs(bv - ov)
        if err >= tol:
            failures.append(
                f"{key}: |diff|={err:.3e} >= tol={tol:.3e} "
                f"(actual={bv:.6e}, oracle={ov:.6e})"
            )

    # cfl_peak — skip gracefully if absent
    if "cfl_peak" in b_diag:
        bv  = float(b_diag["cfl_peak"])
        ov  = float(oracle["cfl_peak"])
        tol = scalar_tol(ov, s_abs, s_rel, sf)
        err = abs(bv - ov)
        if err >= tol:
            failures.append(
                f"cfl_peak: |diff|={err:.3e} >= tol={tol:.3e} "
                f"(actual={bv:.6e}, oracle={ov:.6e})"
            )

    return failures


def check_field_summaries_and_elements(vector: dict, body: dict) -> list[str]:
    tols    = vector["tolerances"]
    o_fsum  = vector["field_summaries"]
    ff      = backend_factor(tols, "field")
    fs_abs  = tols["field_summary_abs"]
    fs_rel  = tols["field_summary_rel"]
    fe_abs  = tols["field_element_abs"]
    fe_rel  = tols["field_element_rel"]
    nx      = float(vector["parameters"]["grid_nx"])
    ny      = float(vector["parameters"]["grid_ny"])
    sqrt_n  = math.sqrt(nx * ny)
    use_sqn = tols.get("sqrt_n_scaling", False)

    fields = [
        ("u_velocity", "u_velocity"),
        ("v_velocity", "v_velocity"),
        ("vorticity",  "vorticity"),
        ("pressure",   "pressure"),
    ]
    failures = []

    for oracle_key, body_key in fields:
        o_fs    = o_fsum[oracle_key]
        b_field = body.get(body_key, [])

        if not isinstance(b_field, list) or not b_field:
            failures.append(f"{body_key}: not a non-empty float array")
            continue

        b_min  = min(b_field)
        b_max  = max(b_field)
        b_mean = sum(b_field) / len(b_field)

        for stat, bv, ov_key in [
            ("min",  b_min,  "min"),
            ("max",  b_max,  "max"),
            ("mean", b_mean, "mean"),
        ]:
            ov  = float(o_fs[ov_key])
            tol = field_tol(ov, fs_abs, fs_rel, ff, sqrt_n, use_sqn)
            err = abs(bv - ov)
            if err >= tol:
                failures.append(
                    f"{oracle_key}.{stat}: |diff|={err:.3e} >= tol={tol:.3e} "
                    f"(actual={bv:.6e}, oracle={ov:.6e})"
                )

        # Element checks: first_10 / last_10
        o_first = o_fs["first_10"]
        o_last  = o_fs["last_10"]
        b_first = b_field[:len(o_first)]
        b_last  = b_field[len(b_field) - len(o_last):]

        for slot, b_elems, o_elems in [
            ("first_10", b_first, o_first),
            ("last_10",  b_last,  o_last),
        ]:
            for i, (bv, ov) in enumerate(zip(b_elems, o_elems)):
                tol = field_tol(ov, fe_abs, fe_rel, ff, sqrt_n, use_sqn)
                err = abs(bv - ov)
                if err >= tol:
                    failures.append(
                        f"{oracle_key}.{slot}[{i}]: |diff|={err:.3e} >= tol={tol:.3e} "
                        f"(actual={bv:.6e}, oracle={ov:.6e})"
                    )

    return failures


def run_vector(filename: str) -> tuple[str, list[str]]:
    """POST vector to backend, return (status, failures). status is 'skip', 'pass', or 'fail'."""
    vector  = load_vector(filename)
    payload = params_to_request(vector["parameters"])

    try:
        resp = requests.post(BACKEND_URL, json=payload, timeout=120)
    except requests.exceptions.ConnectionError as e:
        return "skip", [f"SKIP: server not reachable at {BACKEND_URL}: {e}"]

    if not resp.ok:
        return "fail", [f"HTTP {resp.status_code}: {resp.text[:200]}"]

    body = resp.json()

    failures = check_scalars(vector, body)
    failures += check_field_summaries_and_elements(vector, body)
    return ("pass" if not failures else "fail"), failures


# ── Parametrised test cases ────────────────────────────────────────────────────

VECTOR_IDS = [f.name for f in VECTOR_FILES]


@pytest.mark.parametrize("filename", VECTOR_IDS)
def test_golden_vector(filename: str) -> None:
    status, failures = run_vector(filename)
    if status == "skip":
        if require_server():
            pytest.fail(failures[0])
        else:
            pytest.skip(failures[0])
    elif status == "fail":
        pytest.fail(f"{filename}: {len(failures)} failure(s):\n" +
                    "\n".join(f"  {f}" for f in failures))
