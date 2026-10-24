"""
test_canonical_physics_cp7fc.py -- CP-7F-C co-rotating vortex-pair merger, ONE
specification run against the Python oracle and all five backends.

Specification: MASTER-BRIEF.md s3.3, amendment E6 (2026-10-07), which replaces the
E3 windows and the six per-suite implementations. Why: finding
development record, finding 2026-10-06 -- the five backend detectors
could not see a merged vortex, the oracle's Gaussian variant started from a
single vorticity maximum, and the six suites ran six different experiments.

WHAT IS FIXED HERE, ONCE
  Configuration  two near-top-hat vortices omega ~ exp(-(r/a)^8), each calibrated
                 to circulation Gamma = 1 on the grid; L = 1, a = 0.05 L, 256^2,
                 dt = 0.0005, nu = 1/5000, so Re_Gamma = Gamma/nu = 5000 (the
                 classification is Re_Gamma-dependent: at 1837 a 3.5a pair merges).
                 Separations d/a = 2.5, 3.0, 3.3, 3.5, 4.0 along x, centred.
  Time unit      tau = pi a^2 / Gamma; one sample per tau, horizon 80 tau.
  Detector       the field is smoothed by a periodic Gaussian of width 0.5a, so
                 the flat top and filament ripples of a patch are not peaks; the
                 local maxima above half the smoothed maximum are counted. A
                 sample reads MERGED when there is one maximum, or two closer
                 than 0.5a. The pair has merged when two consecutive samples read
                 MERGED; t_merge is the first of them.
  Assertions     2.5a and 3.0a merge; 3.5a and 4.0a do not merge by 80 tau;
                 t_merge(2.5a) <= t_merge(3.0a); 3.3a is reported, not asserted
                 (the boundary). No absolute-tau window (task 04-C, 2026-07-04).
  Anchor         Melander, Zabusky & McWilliams (1988) JFM 195:303-340 and
                 Saffman (1992) s7.4: d_c ~ 3.3a for uniform patches. Measured on
                 2026-10-06 with this configuration: merger at ~7 / 13 / 50 tau for
                 2.5a / 3.0a / 3.3a, none at 3.5a and 4.0a.

The detector itself is tested first (test_cp7fc_detector_*), on synthetic fields
with known answers, because a detector that cannot fire is what this replaces.

Targets: "oracle" runs shared/physics/kh_physics.py in-process; the five backends
are reached over HTTP at conftest.BACKEND_URLS (KH_<B>_URL overrides). With
KH_GOLDEN_REQUIRE_SERVER=1 an unreachable backend FAILS instead of skipping.
KH_CP7F_OUT=<dir> writes one JSON record per target (the per-sample readings).

Cost: about 230 tau of simulation per target (merged cases stop early), i.e.
~3 700 steps at 256^2 -- a few minutes per backend.
"""
import json
import math
import os
import pathlib
import sys

import numpy as np
import pytest
import requests

from conftest import BACKEND_NAMES, BACKEND_URLS

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[2] / "shared" / "physics"))
import kh_physics  # noqa: E402

# --------------------------------------------------------------------------- spec (E6)
N = 256
L = 1.0
A_CORE = 0.05 * L
GAMMA = 1.0
DT = 0.0005
RE = 5000.0
U = 1.0                                   # nu = U / Re = 1/5000
TAU = math.pi * A_CORE ** 2 / GAMMA
HORIZON_TAU = 80
SEPARATIONS = (2.5, 3.0, 3.3, 3.5, 4.0)
MUST_MERGE = (2.5, 3.0)
MUST_NOT_MERGE = (3.5, 4.0)
SMOOTH_SIGMA = 0.5 * A_CORE
MERGE_SEP = 0.5 * A_CORE
CONSECUTIVE = 2

TIMEOUT = 600


# --------------------------------------------------------------------------- IC
def pair_ic(d_over_a: float, n: int = N) -> np.ndarray:
    """Two exp(-(r/a)^8) vortices at (L/2 -/+ d/2, L/2), each with circulation
    GAMMA on the grid. Flat, x-major: index i*n + j, x = i*dx."""
    dx = L / n
    x = np.arange(n) * dx
    X, Y = np.meshgrid(x, x, indexing="ij")
    d = d_over_a * A_CORE
    f1 = np.exp(-(np.hypot(X - (L / 2 - d / 2), Y - L / 2) / A_CORE) ** 8)
    f2 = np.exp(-(np.hypot(X - (L / 2 + d / 2), Y - L / 2) / A_CORE) ** 8)
    return (GAMMA / (f1.sum() * dx * dx) * (f1 + f2)).ravel()


# --------------------------------------------------------------------------- detector
def smoothed_peaks(omega_flat, n: int = N, lx: float = L, sigma: float = SMOOTH_SIGMA):
    """(number of local maxima above half the smoothed maximum, separation of the
    two highest). Periodic in both directions."""
    w = np.asarray(omega_flat, dtype=float).reshape(n, n)
    k = 2.0 * np.pi * np.fft.fftfreq(n, d=lx / n)
    KX, KY = np.meshgrid(k, k, indexing="ij")
    s = np.real(np.fft.ifft2(np.fft.fft2(w) * np.exp(-0.5 * sigma ** 2 * (KX ** 2 + KY ** 2))))
    nb = [np.roll(np.roll(s, i, 0), j, 1) for i in (-1, 0, 1) for j in (-1, 0, 1) if i or j]
    is_max = np.all([s > q for q in nb], axis=0) & (s >= 0.5 * s.max())
    ii, jj = np.nonzero(is_max)
    order = np.argsort(s[ii, jj])[::-1]
    ii, jj = ii[order], jj[order]
    if len(ii) < 2:
        return int(len(ii)), 0.0
    dx = lx / n
    ddx = (ii[0] - ii[1]) * dx
    ddy = (jj[0] - jj[1]) * dx
    ddx -= round(ddx / lx) * lx
    ddy -= round(ddy / lx) * lx
    return int(len(ii)), float(math.hypot(ddx, ddy))


def reads_merged(omega_flat, n: int = N, lx: float = L, a: float = A_CORE) -> bool:
    npk, sep = smoothed_peaks(omega_flat, n, lx, 0.5 * a)
    return npk <= 1 or sep < 0.5 * a


# --------------------------------------------------------------------------- engines
def _request(omega, steps: int) -> dict:
    return {
        "grid_nx": N, "grid_ny": N, "domain_lx": L, "domain_ly": L,
        "dt": DT, "steps": int(steps), "reynolds_number": RE, "velocity_shear": U,
        "perturbation_amplitude": 0.0, "perturbation_mode": 0,
        "initial_omega": [float(v) for v in omega],
    }


class _Oracle:
    name = "oracle"

    def advance(self, omega, steps):
        req = kh_physics.SimulationRequest(**_request(omega, steps))
        return np.asarray(kh_physics.simulate(req).vorticity, dtype=float).ravel()


class _Http:
    def __init__(self, name, url):
        self.name, self.url, self.s = name, url, requests.Session()

    def advance(self, omega, steps):
        r = self.s.post(f"{self.url}/simulate", json=_request(omega, steps), timeout=TIMEOUT)
        assert r.status_code == 200, f"{self.name}: HTTP {r.status_code}: {r.text[:300]}"
        body = r.json()
        assert body.get("steps_completed") == steps, \
            f"{self.name}: steps_completed {body.get('steps_completed')} != {steps}"
        w = np.asarray(body["vorticity"], dtype=float).ravel()
        assert w.size == N * N and np.all(np.isfinite(w)), f"{self.name}: bad vorticity field"
        return w


def _engine(target):
    if target == "oracle":
        return _Oracle()
    url = BACKEND_URLS[target]
    try:
        requests.get(f"{url}/health", timeout=5).raise_for_status()
    except requests.RequestException as exc:
        if os.environ.get("KH_GOLDEN_REQUIRE_SERVER") == "1":
            pytest.fail(f"{target} not reachable at {url}: {exc}")
        pytest.skip(f"{target} not reachable at {url}")
    return _Http(target, url)


def run_case(engine, d_over_a: float) -> dict:
    """Advance one tau per call to HORIZON_TAU; stop once the merger is confirmed."""
    omega = pair_ic(d_over_a)
    done, run, t_merge, samples = 0, 0, None, []
    for k in range(1, HORIZON_TAU + 1):
        target_step = int(round(k * TAU / DT))
        omega = engine.advance(omega, target_step - done)
        done = target_step
        npk, sep = smoothed_peaks(omega)
        merged = npk <= 1 or sep < MERGE_SEP
        samples.append({"t_tau": k, "peaks": npk, "sep_a": round(sep / A_CORE, 4), "merged": merged})
        run = run + 1 if merged else 0
        if run == CONSECUTIVE:
            t_merge = k - CONSECUTIVE + 1
            break
    return {"d_over_a": d_over_a, "t_merge_tau": t_merge, "horizon_tau": HORIZON_TAU,
            "steps": done, "samples": samples}


# --------------------------------------------------------------------------- detector tests
def test_cp7fc_detector_two_patches_are_two():
    for d in (2.5, 3.0, 3.5, 4.0):
        npk, sep = smoothed_peaks(pair_ic(d))
        assert npk == 2 and abs(sep / A_CORE - d) < 0.15, f"d={d}a read as {npk} peaks, sep {sep / A_CORE:.2f}a"
        assert not reads_merged(pair_ic(d))


def test_cp7fc_detector_one_patch_is_merged_wherever_it_sits():
    # The half-domain detector this replaces could not fire for a vortex sitting
    # on x = L/2; this one must, at any position.
    dx = L / N
    x = np.arange(N) * dx
    X, Y = np.meshgrid(x, x, indexing="ij")
    for cx, cy in ((0.5, 0.5), (0.3, 0.7), (0.02, 0.5)):
        ddx = X - cx; ddx -= np.round(ddx / L) * L
        ddy = Y - cy; ddy -= np.round(ddy / L) * L
        w = np.exp(-(np.hypot(ddx, ddy) / (1.4 * A_CORE)) ** 8).ravel()
        assert reads_merged(w), f"single vortex at ({cx},{cy}) not read as merged"


def test_cp7fc_detector_is_rotation_invariant():
    # A vertical pair must read like a horizontal one (the old detector read a
    # vertical pair as merged).
    w = pair_ic(3.0).reshape(N, N).T.ravel()
    npk, sep = smoothed_peaks(w)
    assert npk == 2 and abs(sep / A_CORE - 3.0) < 0.15


# --------------------------------------------------------------------------- the testbed
@pytest.mark.slow
@pytest.mark.parametrize("target", ["oracle"] + BACKEND_NAMES)
def test_cp7fc_vortex_pair_merger(target):
    engine = _engine(target)
    results = {d: run_case(engine, d) for d in SEPARATIONS}

    out = os.environ.get("KH_CP7F_OUT")
    if out:
        pathlib.Path(out).mkdir(parents=True, exist_ok=True)
        with open(pathlib.Path(out) / f"cp7fc-{target}.json", "w") as fh:
            json.dump({"target": target, "spec": "MASTER s3.3 E6", "tau": TAU,
                       "cases": [results[d] for d in SEPARATIONS]}, fh, indent=1)
    line = ", ".join(f"{d}a: {results[d]['t_merge_tau'] if results[d]['t_merge_tau'] else 'no merge'}"
                     for d in SEPARATIONS)
    print(f"CP-7F-C {target}: t_merge/tau  {line}")

    for d in MUST_MERGE:
        assert results[d]["t_merge_tau"] is not None, \
            f"{target}: d={d}a did not merge by {HORIZON_TAU} tau ({line})"
    for d in MUST_NOT_MERGE:
        assert results[d]["t_merge_tau"] is None, \
            f"{target}: d={d}a merged at {results[d]['t_merge_tau']} tau, expected none by {HORIZON_TAU} ({line})"
    assert results[2.5]["t_merge_tau"] <= results[3.0]["t_merge_tau"], \
        f"{target}: closer pair merged later ({line})"
