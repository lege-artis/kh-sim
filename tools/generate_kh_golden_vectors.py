#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
"""
generate_kh_golden_vectors.py
Generate kh-sim/shared/golden-vectors/*.json from the Python reference oracle.

Per phase 6a design: _audit/PHASE-6A-GOLDEN-VECTOR-MATRIX-DESIGN-v0.1.md §2-§7

Usage:
    python tools/generate_kh_golden_vectors.py [--out DIR] [--only V1,V2,...] [--force]

Exit codes:
    0 — all vectors generated + schema-valid
    1 — generation error (oracle ran but failed)
    2 — schema validation failed
    3 — output directory not writable
"""
from __future__ import annotations

import argparse
import datetime
import json
import pathlib
import subprocess
import sys
import time
import warnings

import numpy as np

# ── Locate oracle ──────────────────────────────────────────────────────────────

HERE  = pathlib.Path(__file__).resolve().parent
REPO  = HERE.parent
sys.path.insert(0, str(REPO))

from shared.physics.kh_physics import (  # noqa: E402
    simulate,
    SimulationRequest,
)

DEFAULT_OUT  = REPO / "shared" / "golden-vectors"
SCHEMA_PATH  = DEFAULT_OUT / "schema.json"
ORACLE_PATH  = "shared/physics/kh_physics.py"


# ── 9-vector matrix (phase 6a §2) ─────────────────────────────────────────────

VECTOR_MATRIX = [
    {
        "id":       "V1",
        "filename": "v1-canonical-reference.json",
        "purpose":  "Canonical reference. Phase-5 baseline.",
        "params": {
            "grid_nx": 64, "grid_ny": 32, "domain_lx": 1.0, "domain_ly": 0.5,
            "dt": 0.001, "nsteps": 100,
            "reynolds_number": 1000.0, "velocity_shear": 1.0,
            "perturbation_amplitude": 0.01, "perturbation_mode": 2,
        },
    },
    {
        "id":       "V2",
        "filename": "v2-low-re-diffusive.json",
        "purpose":  "Low-Re diffusive. nu=0.01; strong viscous damping; KE decreases monotonically.",
        "params": {
            "grid_nx": 64, "grid_ny": 32, "domain_lx": 1.0, "domain_ly": 0.5,
            "dt": 0.001, "nsteps": 100,
            "reynolds_number": 100.0, "velocity_shear": 1.0,
            "perturbation_amplitude": 0.01, "perturbation_mode": 2,
        },
    },
    {
        "id":       "V3",
        "filename": "v3-high-re-nonlinear.json",
        "purpose":  "High-Re nearly-inviscid. nu=0.0001; nonlinear cascade; dealiasing stress test.",
        "params": {
            "grid_nx": 64, "grid_ny": 32, "domain_lx": 1.0, "domain_ly": 0.5,
            "dt": 0.001, "nsteps": 100,
            "reynolds_number": 10000.0, "velocity_shear": 1.0,
            "perturbation_amplitude": 0.01, "perturbation_mode": 2,
        },
    },
    {
        "id":       "V4",
        "filename": "v4-double-resolution.json",
        "purpose":  "Double-resolution refinement. Same T_end as V1; 4x spatial points; convergence check.",
        "params": {
            "grid_nx": 128, "grid_ny": 64, "domain_lx": 1.0, "domain_ly": 0.5,
            "dt": 0.0005, "nsteps": 200,
            "reynolds_number": 1000.0, "velocity_shear": 1.0,
            "perturbation_amplitude": 0.01, "perturbation_mode": 2,
        },
    },
    {
        "id":       "V5",
        "filename": "v5-half-resolution-fast.json",
        "purpose":  "Half-resolution fast. Same T_end as V1; 1/4 spatial points; fast CI unit-test.",
        "params": {
            "grid_nx": 32, "grid_ny": 16, "domain_lx": 1.0, "domain_ly": 0.5,
            "dt": 0.0005, "nsteps": 200,
            "reynolds_number": 1000.0, "velocity_shear": 1.0,
            "perturbation_amplitude": 0.01, "perturbation_mode": 2,
        },
    },
    {
        "id":       "V6",
        "filename": "v6-long-time-evolution.json",
        "purpose":  "Long-time KH evolution. T_end=0.5; nonlinear saturation phase; vortex rollup.",
        "params": {
            "grid_nx": 64, "grid_ny": 32, "domain_lx": 1.0, "domain_ly": 0.5,
            "dt": 0.001, "nsteps": 500,
            "reynolds_number": 1000.0, "velocity_shear": 1.0,
            "perturbation_amplitude": 0.01, "perturbation_mode": 2,
        },
    },
    {
        "id":       "V7",
        "filename": "v7-alt-mode-k1.json",
        "purpose":  "Alternate mode k=1. Different unstable wavenumber; checks scheme symmetry.",
        "params": {
            "grid_nx": 64, "grid_ny": 32, "domain_lx": 1.0, "domain_ly": 0.5,
            "dt": 0.001, "nsteps": 100,
            "reynolds_number": 1000.0, "velocity_shear": 1.0,
            "perturbation_amplitude": 0.01, "perturbation_mode": 1,
        },
    },
    {
        "id":       "V8",
        "filename": "v8-dealiasing-edge-k4.json",
        "purpose":  "Dealiasing-edge mode k=4. Energy near 2/3 cutoff; dealiasing precision test.",
        "params": {
            "grid_nx": 64, "grid_ny": 32, "domain_lx": 1.0, "domain_ly": 0.5,
            "dt": 0.001, "nsteps": 100,
            "reynolds_number": 1000.0, "velocity_shear": 1.0,
            "perturbation_amplitude": 0.01, "perturbation_mode": 4,
        },
    },
    {
        "id":       "V9",
        "filename": "v9-zero-perturbation.json",
        "purpose":  "Zero-perturbation base-flow. amp=0; pure tanh shear IC; no KH growth; integrator stability test.",
        "params": {
            "grid_nx": 64, "grid_ny": 32, "domain_lx": 1.0, "domain_ly": 0.5,
            "dt": 0.001, "nsteps": 100,
            "reynolds_number": 1000.0, "velocity_shear": 1.0,
            "perturbation_amplitude": 0.0, "perturbation_mode": 2,
        },
    },
]


# ── Default tolerances (phase 6a §4) ──────────────────────────────────────────

DEFAULT_TOLERANCES = {
    "scalar_abs":            1e-12,
    "scalar_rel":            1e-9,
    "scalar_div_rms_abs":    1e-11,
    "scalar_div_rms_rel":    1e-8,
    "field_summary_abs":     1e-10,
    "field_summary_rel":     1e-7,
    "field_element_abs":     1e-9,
    "field_element_rel":     1e-7,
    "sqrt_n_scaling":        True,
    "backend_overrides": {
        "kh-scala":  {"scalar_factor": 2.0,  "field_factor": 2.0},
        "kh-pascal": {"scalar_factor": 10.0, "field_factor": 10.0},
    },
}


# ── Oracle runner ─────────────────────────────────────────────────────────────

def _run_vector(params: dict) -> dict:
    """
    Run one golden vector via simulate() and return all data for JSON serialisation.
    cfl_peak is now provided by simulate() directly (Diagnostics.cfl_peak field).
    """
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

    t0     = time.perf_counter()
    result = simulate(req)
    elapsed_ms = (time.perf_counter() - t0) * 1000.0

    nx, ny = params["grid_nx"], params["grid_ny"]
    u   = np.array(result.u_velocity).reshape(nx, ny)
    v   = np.array(result.v_velocity).reshape(nx, ny)
    omega = np.array(result.vorticity).reshape(nx, ny)
    psi = np.array(result.pressure).reshape(nx, ny)

    return {
        "u": u, "v": v, "omega": omega, "psi": psi,
        "ke":       result.diagnostics.kinetic_energy,
        "enstrophy": result.diagnostics.enstrophy,
        "max_vort": result.diagnostics.max_vorticity,
        "div":      result.diagnostics.divergence_rms,
        "cfl_peak": result.diagnostics.cfl_peak,
        "t_final":  result.t_final,
        "elapsed_ms": elapsed_ms,
        "nx": nx, "ny": ny,
    }


def _field_summary(arr: np.ndarray) -> dict:
    flat = arr.flatten()
    return {
        "length":   int(len(flat)),
        "first_10": [float(x) for x in flat[:10]],
        "last_10":  [float(x) for x in flat[-10:]],
        "min":      float(np.min(flat)),
        "max":      float(np.max(flat)),
        "mean":     float(np.mean(flat)),
    }


def _oracle_commit() -> str:
    try:
        return subprocess.check_output(
            ["git", "rev-parse", "HEAD"],
            cwd=str(REPO), stderr=subprocess.DEVNULL
        ).decode().strip()
    except Exception:
        return "unknown"


def _env_block() -> dict:
    try:
        import scipy
        scipy_ver = scipy.__version__
    except ImportError:
        scipy_ver = "unavailable"
    return {
        "python": f"{sys.version_info.major}.{sys.version_info.minor}.{sys.version_info.micro}",
        "numpy":  np.__version__,
        "scipy":  scipy_ver,
    }


# ── JSON emitter ──────────────────────────────────────────────────────────────

def _emit_vector(spec: dict, result: dict, oracle_commit: str, env: dict) -> dict:
    params = spec["params"]
    return {
        "$schema": "kh-sim-golden-vector-v1",
        "metadata": {
            "vector_id":       f"{spec['id']}-{spec['filename'].removeprefix(spec['id'].lower() + '-').removesuffix('.json')}",
            "purpose":         spec["purpose"],
            "generation_date": datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
            "oracle_commit":   oracle_commit,
            "oracle_path":     ORACLE_PATH,
            "scheme":          "RK4 (canonical)",
            "dealiasing":      "2/3-rule (Orszag 1971; KH-NM-4)",
            "cfl_guard":       "advective+diffusive bound (KH-NM-5)",
            "environment":     env,
        },
        "parameters": {
            "grid_nx":                params["grid_nx"],
            "grid_ny":                params["grid_ny"],
            "domain_lx":              params["domain_lx"],
            "domain_ly":              params["domain_ly"],
            "dt":                     params["dt"],
            "nsteps":                 params["nsteps"],
            "reynolds_number":        params["reynolds_number"],
            "velocity_shear":         params["velocity_shear"],
            "perturbation_amplitude": params["perturbation_amplitude"],
            "perturbation_mode":      params["perturbation_mode"],
        },
        "diagnostics": {
            "kinetic_energy":  result["ke"],
            "enstrophy":       result["enstrophy"],
            "max_vorticity":   result["max_vort"],
            "divergence_rms":  result["div"],
            "cfl_peak":        result["cfl_peak"],
            "t_final":         result["t_final"],
        },
        "field_summaries": {
            "u_velocity": _field_summary(result["u"]),
            "v_velocity": _field_summary(result["v"]),
            "vorticity":  _field_summary(result["omega"]),
            "pressure":   _field_summary(result["psi"]),
        },
        "tolerances": DEFAULT_TOLERANCES,
    }


# ── Schema validation ──────────────────────────────────────────────────────────

def _validate_json(doc: dict, schema_path: pathlib.Path) -> None:
    """Validate doc against schema.json; raise if invalid."""
    if not schema_path.exists():
        return   # schema not yet present (generated before schema.json exists)
    try:
        import jsonschema
        with open(schema_path) as fh:
            schema = json.load(fh)
        jsonschema.validate(instance=doc, schema=schema)
    except ImportError:
        # Fallback: manual key-presence check (per SONNET-BRIEF §4)
        _manual_key_check(doc)


def _manual_key_check(doc: dict) -> None:
    required_top = {"$schema", "metadata", "parameters", "diagnostics", "field_summaries", "tolerances"}
    missing = required_top - set(doc.keys())
    if missing:
        raise ValueError(f"manual key-check: missing top-level keys: {missing}")
    diag_required = {"kinetic_energy", "enstrophy", "max_vorticity", "divergence_rms", "cfl_peak", "t_final"}
    missing_diag = diag_required - set(doc["diagnostics"].keys())
    if missing_diag:
        raise ValueError(f"manual key-check: missing diagnostics keys: {missing_diag}")


# ── CLI ────────────────────────────────────────────────────────────────────────

def _parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Generate kh-sim golden vectors from the Python reference oracle."
    )
    p.add_argument("--out",   type=pathlib.Path, default=DEFAULT_OUT,
                   help=f"Output directory (default: {DEFAULT_OUT})")
    p.add_argument("--only",  type=str, default="",
                   help="Comma-separated vector IDs to generate (e.g. V1,V3)")
    p.add_argument("--force", action="store_true",
                   help="Overwrite existing JSON files")
    return p.parse_args()


def main() -> int:
    args = _parse_args()
    out_dir: pathlib.Path = args.out

    # Exit-code 3: output not writable
    try:
        out_dir.mkdir(parents=True, exist_ok=True)
        (out_dir / ".write_test").touch()
        (out_dir / ".write_test").unlink()
    except OSError as exc:
        print(f"ERROR: output directory not writable: {exc}", file=sys.stderr)
        return 3

    only_ids = set(v.strip().upper() for v in args.only.split(",") if v.strip())
    vectors  = [v for v in VECTOR_MATRIX if (not only_ids or v["id"] in only_ids)]

    oracle_commit = _oracle_commit()
    env           = _env_block()

    print(f"kh-sim golden-vector generator")
    print(f"  oracle commit : {oracle_commit[:12]}")
    print(f"  python        : {env['python']}")
    print(f"  numpy         : {env['numpy']}")
    print(f"  scipy         : {env['scipy']}")
    print(f"  out dir       : {out_dir}")
    print(f"  vectors       : {[v['id'] for v in vectors]}")
    print()

    schema_path = out_dir / "schema.json"
    any_failed  = False
    validation_failed = False

    for spec in vectors:
        out_path = out_dir / spec["filename"]
        if out_path.exists() and not args.force:
            print(f"  [{spec['id']}] skip (exists; use --force to overwrite)")
            continue

        print(f"  [{spec['id']}] generating {spec['filename']} ...", end=" ", flush=True)
        t0 = time.perf_counter()
        try:
            result = _run_vector(spec["params"])
        except Exception as exc:
            import traceback
            print(f"FAILED")
            print(f"    ERROR: {exc}")
            traceback.print_exc()
            any_failed = True
            continue

        doc = _emit_vector(spec, result, oracle_commit, env)

        # Schema validation before writing
        try:
            _validate_json(doc, schema_path)
        except Exception as exc:
            print(f"SCHEMA-FAIL")
            print(f"    SCHEMA ERROR: {exc}")
            validation_failed = True
            any_failed = True
            continue

        with open(out_path, "w") as fh:
            json.dump(doc, fh, indent=2)

        wall = time.perf_counter() - t0
        print(
            f"done ({wall:.1f}s)  KE={result['ke']:.6f}  "
            f"enstrophy={result['enstrophy']:.4f}  "
            f"div_rms={result['div']:.2e}  "
            f"cfl_peak={result['cfl_peak']:.4f}"
        )

    print()
    if validation_failed:
        print("FAIL: schema validation failed on one or more vectors.")
        return 2
    if any_failed:
        print("FAIL: generation error on one or more vectors.")
        return 1
    print("All vectors generated successfully.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
