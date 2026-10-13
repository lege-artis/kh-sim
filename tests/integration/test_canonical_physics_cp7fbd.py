"""
test_canonical_physics_cp7fbd.py -- CP-7F-B (forced 2-D turbulence, enstrophy-range
slope) and CP-7F-D (forced Taylor-Green steady state): ONE specification each, run
against the Python oracle and all five backends.

Specification: MASTER-BRIEF.md s3.2 and s3.4, amendment E7 (2026-10-10). Ruling 2
(2026-10-07) replaced the per-language CP-7F suites by one shared HTTP harness; ruling 8
made the Fortran HTTP server compute forcing when that harness arrives. Before it the
six suites ran six different B and D experiments (128^2 probes reported
NATIVE-REQUIRED, a 64^2 D at Re 500..5000 in C++, another at Re 50..1000 in the
oracle), and the Fortran HTTP server refused every forcing request.

CP-7F-D  forced Taylor-Green, a solver-consistency check (ruling D1, 2026-06-15; not C1)
  L = 2 pi, 64^2, dt = 0.005, velocity_shear 1, Re in {50, 100, 200, 500, 1000}.
  Forcing "mode" (1, 1), amplitude f0 = sqrt(2 nu eps), eps = 0.1. Initial field the
  analytic balance A_ss cos x cos y, A_ss = f0 / (2 nu): the single Taylor-Green mode
  has no nonlinear transfer, so A_ss is a fixed point of the forced RK4 map.
  500 steps, then 10 samples 50 steps apart; omega_rms of the returned field.
  Assert: every sample within 1e-6 of A_ss / 2 (the fixed point is kept), and the
  log-log slope of omega_rms over Re in [0.45, 0.55] (MASTER s3.4).

CP-7F-B  forced 2-D turbulence, regime + reproducibility gate (E5, 2026-06-13; not C1)
  CP7F-B-TUNED-PARAMS.json: L = 2 pi, Re 5000 (nu = 2e-4), dt = 0.001, band forcing
  A = 40 on 10 <= |k| <= 12, rotating phases omega0 = 40, linear drag 0.05.
  Initial field: random phases, E(k) ~ k^4 exp(-2 (k/8)^2), omega_rms 1, generated
  here once (seed 20261010) and sent to every target.
  Stage 1: 128^2, t 0 -> 80, eight requests of 10 000 steps.
  Stage 2: the field interpolated spectrally to 256^2, t 80 -> 110, thirty requests
  of 1 000 steps; every request continues the previous one with t_start = its t_final,
  so the rotating forcing keeps its phase (t_start, added 2026-10-10 for this harness).
  Spectrum E(k) = 2 pi k * (mean over the shell round(|k|) = k of 0.5 |omega_hat|^2 / |k|^2),
  the angle-averaged form: a plain shell sum carries the lattice's mode count per shell,
  which jumps by up to 15 % between neighbouring shells and alone breaks the monotone
  guard (measured on the oracle run of 2026-10-10: rises of 1.15 at k 20..63 in the sum,
  none in the average). Sampled at t = 91, 92, ..., 110; mean over (90, 110], windows
  (90, 100] and (100, 110].
  Assert on k in [16, 64]: least-squares log-log slope of the mean in [-4.6, -3.6];
  window slopes within 0.12 (stationarity; 0.05 until 2026-10-11, see B_STATIONARITY);
  E(k+1) <= 1.01 E(k) (monotone).
  Measured by the oracle at these parameters on 2026-06-12: -4.10, stationarity 0.044;
  with this specification on 2026-10-10 (workspace): -4.18, stationarity 0.026, no rise;
  on the station 2026-10-11, all six targets: slopes -4.12 .. -4.16, stationarity
  0.0004 .. 0.057, no rise (RESULTS-CP7FB-20261010-234912, 315 min).
  Cost: 80 000 steps at 128^2 + 30 000 at 256^2 per target - from about an hour
  (oracle, Rust) to several hours (Pascal). Run only with KH_CP7FB=1; resumable with
  KH_CP7F_CACHE=<dir> (one state file per target after every request).

Targets: "oracle" runs shared/physics/kh_physics.py in-process; the five backends are
reached over HTTP at conftest.BACKEND_URLS (KH_<B>_URL overrides). With
KH_GOLDEN_REQUIRE_SERVER=1 an unreachable backend FAILS instead of skipping.
KH_CP7F_OUT=<dir> writes one JSON record per target and test (samples, spectra,
per-request diagnostics), whether the assertions pass or not.
"""
import json
import math
import os
import pathlib
import sys
import time
import warnings

import numpy as np
import pytest
import requests

from conftest import BACKEND_NAMES, BACKEND_URLS

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[2] / "shared" / "physics"))
import kh_physics  # noqa: E402

TWO_PI = 2.0 * math.pi
TIMEOUT = 7200


# --------------------------------------------------------------------------- engines
class _Oracle:
    name = "oracle"

    def advance(self, req: dict) -> dict:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            res = kh_physics.simulate(kh_physics.SimulationRequest(**req))
        return {"vorticity": res.vorticity, "t_final": res.t_final,
                "steps_completed": res.steps_completed, "warnings": res.warnings,
                "diagnostics": dict(vars(res.diagnostics)), "compute_time_ms": res.compute_time_ms}


class _Http:
    # One connection per request, as test_conformance.py: the Pascal server
    # (fphttpserver) resets a kept-alive connection that a Session reuses
    # ("Connection reset by peer", seen 2026-10-10 on a 64^2 D run).
    def __init__(self, name, url):
        self.name, self.url = name, url

    def advance(self, req: dict) -> dict:
        r = requests.post(f"{self.url}/simulate", json=req, timeout=TIMEOUT)
        assert r.status_code == 200, f"{self.name}: HTTP {r.status_code}: {r.text[:300]}"
        return r.json()


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


def _field(body: dict, n: int, req: dict, who: str) -> np.ndarray:
    assert body.get("steps_completed") == req["steps"], \
        f"{who}: steps_completed {body.get('steps_completed')} != {req['steps']}"
    want_t = req.get("t_start", 0.0) + req["steps"] * req["dt"]
    assert abs(body.get("t_final", float("nan")) - want_t) < 1e-9, \
        f"{who}: t_final {body.get('t_final')} != t_start + steps*dt = {want_t}"
    w = np.asarray(body["vorticity"], dtype=float).reshape(n, n)
    assert np.all(np.isfinite(w)), f"{who}: non-finite vorticity"
    return w


def _write(name: str, record: dict) -> None:
    out = os.environ.get("KH_CP7F_OUT")
    if out:
        pathlib.Path(out).mkdir(parents=True, exist_ok=True)
        with open(pathlib.Path(out) / name, "w") as fh:
            json.dump(record, fh, indent=1)


def _slope(x, y) -> float:
    return float(np.polyfit(np.log(np.asarray(x, float)), np.log(np.asarray(y, float)), 1)[0])


# =========================================================================== CP-7F-D
D_N = 64
D_DT = 0.005
D_EPS = 0.1
D_RE = (50.0, 100.0, 200.0, 500.0, 1000.0)
D_SPIN, D_EVERY, D_SAMPLES = 500, 50, 10
D_FIXED_POINT_REL = 1e-6
D_SLOPE = (0.45, 0.55)


def d_case(engine, re: float) -> dict:
    nu = 1.0 / re
    f0 = math.sqrt(2.0 * nu * D_EPS)
    a_ss = f0 / (2.0 * nu)
    x = np.arange(D_N) * (TWO_PI / D_N)
    X, Y = np.meshgrid(x, x, indexing="ij")
    omega = a_ss * np.cos(X) * np.cos(Y)
    base = {"grid_nx": D_N, "grid_ny": D_N, "domain_lx": TWO_PI, "domain_ly": TWO_PI,
            "dt": D_DT, "reynolds_number": re, "velocity_shear": 1.0,
            "perturbation_amplitude": 0.0, "perturbation_mode": 1,
            "forcing_type": "mode", "forcing_amplitude": f0, "forcing_kx": 1, "forcing_ky": 1}
    t, samples = 0.0, []
    for k in range(D_SAMPLES + 1):
        req = dict(base, steps=D_SPIN if k == 0 else D_EVERY, t_start=t,
                   initial_omega=[float(v) for v in omega.ravel()])
        body = engine.advance(req)
        omega = _field(body, D_N, req, f"{engine.name} Re={re}")
        t = body["t_final"]
        if k:
            samples.append(float(np.sqrt(np.mean(omega ** 2))))
    return {"re": re, "f0": f0, "a_ss": a_ss, "expected_rms": a_ss / 2.0, "samples": samples,
            "omega_rms": float(np.mean(samples))}


@pytest.mark.parametrize("target", ["oracle"] + BACKEND_NAMES)
def test_cp7fd_forced_taylor_green(target):
    engine = _engine(target)
    cases = [d_case(engine, re) for re in D_RE]
    slope = _slope(D_RE, [c["omega_rms"] for c in cases])
    _write(f"cp7fd-{target}.json", {"target": target, "spec": "MASTER s3.4 E7", "slope": slope,
                                    "cases": cases})
    print(f"CP-7F-D {target}: slope {slope:.6f}; omega_rms "
          + ", ".join(f"Re {int(c['re'])}: {c['omega_rms']:.6f}" for c in cases))
    for c in cases:
        dev = max(abs(s - c["expected_rms"]) for s in c["samples"]) / c["expected_rms"]
        assert dev <= D_FIXED_POINT_REL, (
            f"{target}: Re={c['re']:g} omega_rms left the forced fixed point A_ss/2 = "
            f"{c['expected_rms']:.9f} by {dev:.2e} (samples {c['samples']})")
    assert D_SLOPE[0] <= slope <= D_SLOPE[1], f"{target}: slope {slope:.4f} outside {D_SLOPE}"


# =========================================================================== CP-7F-B
B_PARAMS = {"domain_lx": TWO_PI, "domain_ly": TWO_PI, "dt": 0.001, "reynolds_number": 5000.0,
            "velocity_shear": 1.0, "perturbation_amplitude": 0.0, "perturbation_mode": 1,
            "forcing_type": "band", "forcing_amplitude": 40.0, "forcing_kmin": 10.0,
            "forcing_kmax": 12.0, "forcing_omega0": 40.0, "linear_drag": 0.05}
B_SEED = 20261010
B_K0 = 8.0
B_STAGE1 = (128, 8, 10000)      # grid, requests, steps per request: t 0 -> 80
B_STAGE2 = (256, 30, 1000)      # t 80 -> 110
B_SAMPLE_FROM = 90.0            # spectra at t = 91 ... 110
B_FIT = (16, 64)
B_SLOPE = (-4.6, -3.6)
# Window-to-window slope agreement, |slope(w1) - slope(w2)|. Each window averages 10
# spectra one time unit apart; a single spectrum's slope scatters by 0.072 (sd, pooled
# over the six station targets of 2026-10-11) with lag-1 autocorrelation ~ +0.2, so a
# window slope carries ~0.027 and the difference of two ~0.038 (1 sigma) on a flow that
# IS stationary. The E5 bound 0.05 (set from the one oracle value 0.044) is 1.3 sigma:
# a stationary run failed it with ~20 % probability per target, and Fortran did (0.057;
# its slope -4.16 sits inside the other five). 0.12 = 3 sigma of the procedure, R5:
# the tolerance of a procedure observable is derived from the procedure, not chosen.
# (Ruling 2026-10-11.)
B_STATIONARITY = 0.12
B_MONOTONE = 1.01


def b_initial_field(n: int = B_STAGE1[0]) -> np.ndarray:
    """Random-phase vorticity, E(k) ~ k^4 exp(-2 (k/k0)^2), zero mean, omega_rms 1."""
    rng = np.random.Generator(np.random.PCG64(B_SEED))
    k1 = np.fft.fftfreq(n, d=1.0 / n)
    KX, KY = np.meshgrid(k1, k1, indexing="ij")
    k = np.hypot(KX, KY)
    amp = np.zeros_like(k)
    nz = k > 0
    # |omega_hat|^2 per mode ~ k^2 E(k) / (2 pi k)
    amp[nz] = np.sqrt(k[nz] ** 4 * np.exp(-2.0 * (k[nz] / B_K0) ** 2) * k[nz] / (2 * math.pi))
    w = np.real(np.fft.ifft2(amp * np.exp(1j * rng.uniform(0.0, 2 * math.pi, k.shape))))
    w -= w.mean()
    return w / np.sqrt(np.mean(w ** 2))


def spectral_upsample(w: np.ndarray, m: int) -> np.ndarray:
    """Exact interpolation of a periodic field from n^2 to m^2 points (zero padding)."""
    n = w.shape[0]
    wh = np.fft.fft2(w) / n ** 2
    out = np.zeros((m, m), dtype=complex)
    idx = np.r_[0:n // 2, m - n // 2:m]
    src = np.r_[0:n // 2, n - n // 2:n]
    out[np.ix_(idx, idx)] = wh[np.ix_(src, src)]
    return np.real(np.fft.ifft2(out)) * m ** 2


def _shells(n: int):
    k1 = np.fft.fftfreq(n, d=1.0 / n)
    KX, KY = np.meshgrid(k1, k1, indexing="ij")
    k2 = KX ** 2 + KY ** 2
    return k2, np.rint(np.sqrt(k2)).astype(int)


def shell_sums(w: np.ndarray) -> np.ndarray:
    """Sum of 0.5 |omega_hat|^2 / |k|^2 over the integer shells round(|k|), L = 2 pi.
    Sums to the kinetic energy per unit area (Parseval)."""
    n = w.shape[0]
    wh = np.fft.fft2(w) / n ** 2
    k2, shells = _shells(n)
    e = np.zeros_like(k2)
    e[k2 > 0] = 0.5 * np.abs(wh[k2 > 0]) ** 2 / k2[k2 > 0]
    return np.bincount(shells.ravel(), weights=e.ravel())[: n // 2 + 1]


def shell_spectrum(w: np.ndarray) -> np.ndarray:
    """Angle-averaged E(k) = 2 pi k * (shell sum / number of lattice modes in the shell):
    the shell sum without the lattice's mode-count fluctuation."""
    n = w.shape[0]
    _, shells = _shells(n)
    count = np.bincount(shells.ravel())[: n // 2 + 1].astype(float)
    k = np.arange(n // 2 + 1, dtype=float)
    E = np.zeros(n // 2 + 1)
    nz = (count > 0) & (k > 0)
    E[nz] = 2.0 * math.pi * k[nz] * shell_sums(w)[nz] / count[nz]
    return E


def b_fit(E: np.ndarray) -> float:
    k = np.arange(B_FIT[0], B_FIT[1] + 1)
    return _slope(k, E[k])


def _cache(target: str):
    d = os.environ.get("KH_CP7F_CACHE")
    return pathlib.Path(d) / f"cp7fb-{target}-state.npz" if d else None


def b_run(engine, log=print) -> dict:
    """Run the CP-7F-B experiment on one target; resumes from KH_CP7F_CACHE."""
    cache = _cache(engine.name)
    state = {"i": 0, "t": 0.0, "omega": b_initial_field(), "spectra": [], "times": [], "rows": []}
    if cache and cache.exists():
        z = np.load(cache, allow_pickle=False)
        state = {"i": int(z["i"]), "t": float(z["t"]), "omega": z["omega"],
                 "spectra": [list(r) for r in z["spectra"]], "times": list(z["times"]),
                 "rows": json.loads(str(z["rows"]))}
        log(f"CP-7F-B {engine.name}: resumed at request {state['i']}, t = {state['t']}")
    plan = [(B_STAGE1[0], B_STAGE1[2])] * B_STAGE1[1] + [(B_STAGE2[0], B_STAGE2[2])] * B_STAGE2[1]
    while state["i"] < len(plan):
        n, steps = plan[state["i"]]
        omega = state["omega"]
        if omega.shape[0] != n:
            omega = spectral_upsample(omega, n)
        req = dict(B_PARAMS, grid_nx=n, grid_ny=n, steps=steps, t_start=state["t"],
                   initial_omega=[float(v) for v in omega.ravel()])
        t0 = time.time()
        body = engine.advance(req)
        omega = _field(body, n, req, f"{engine.name} request {state['i']}")
        d = body.get("diagnostics", {})
        row = {"i": state["i"], "grid": n, "t_final": body["t_final"], "wall_s": round(time.time() - t0, 1),
               "kinetic_energy": d.get("kinetic_energy"), "enstrophy": d.get("enstrophy"),
               "cfl_peak": d.get("cfl_peak"), "omega_rms": float(np.sqrt(np.mean(omega ** 2))),
               "warnings": body.get("warnings", [])}
        state["rows"].append(row)
        state["t"], state["omega"], state["i"] = float(body["t_final"]), omega, state["i"] + 1
        if n == B_STAGE2[0] and state["t"] > B_SAMPLE_FROM + 1e-9:
            state["spectra"].append([float(v) for v in shell_spectrum(omega)])
            state["times"].append(state["t"])
        log(f"CP-7F-B {engine.name}: request {row['i'] + 1}/{len(plan)} t = {row['t_final']:g} "
            f"omega_rms {row['omega_rms']:.4f} cfl {row['cfl_peak']} ({row['wall_s']} s)")
        if cache:
            cache.parent.mkdir(parents=True, exist_ok=True)
            np.savez(cache, i=state["i"], t=state["t"], omega=state["omega"],
                     spectra=np.asarray(state["spectra"], float) if state["spectra"] else np.zeros((0, 0)),
                     times=np.asarray(state["times"], float), rows=json.dumps(state["rows"]))
    S = np.asarray(state["spectra"])
    T = np.asarray(state["times"])
    mean, w1, w2 = S.mean(axis=0), S[T <= 100.0 + 1e-9].mean(axis=0), S[T > 100.0 + 1e-9].mean(axis=0)
    k = np.arange(B_FIT[0], B_FIT[1])
    rises = [int(kk) for kk in k if mean[kk + 1] > B_MONOTONE * mean[kk]]
    return {"slope": b_fit(mean), "slope_w1": b_fit(w1), "slope_w2": b_fit(w2),
            "stationarity": abs(b_fit(w1) - b_fit(w2)), "monotone_rises_at_k": rises,
            "samples": len(T), "times": T.tolist(), "E_mean": mean.tolist(),
            "E_w1": w1.tolist(), "E_w2": w2.tolist(), "requests": state["rows"]}


def test_cp7fb_helpers():
    """The shell sums add up to the kinetic energy (Parseval); the interpolation is exact
    for a band-limited field; a field built to E(k) = k^-4 per mode reads as exactly
    k^-4, monotone, in the angle-averaged spectrum, while its shell sums rise between
    neighbouring shells (the lattice's mode count); the slope fit recovers a power law."""
    w = b_initial_field(64)
    req = kh_physics.SimulationRequest(grid_nx=64, grid_ny=64, domain_lx=TWO_PI, domain_ly=TWO_PI,
                                       steps=0, perturbation_amplitude=0.0,
                                       initial_omega=[float(v) for v in w.ravel()])
    ke = kh_physics.simulate(req).diagnostics.kinetic_energy
    S = shell_sums(w)
    assert abs(S.sum() - ke) <= 1e-12 * ke, (S.sum(), ke)
    up = spectral_upsample(w, 128)
    assert np.max(np.abs(up[::2, ::2] - w)) < 1e-12
    assert abs(shell_sums(up).sum() - S.sum()) <= 1e-12 * S.sum()
    n = 256
    k2, shells = _shells(n)
    k = np.sqrt(k2)
    amp = np.zeros_like(k)
    nz = k > 0
    amp[nz] = np.sqrt(2.0 * k2[nz] * k[nz] ** -4 / (2 * math.pi * k[nz])) * n ** 2
    # Hermitian phases (those of a real random field), so every mode has exactly amp.
    rng = np.random.Generator(np.random.PCG64(1))
    phase = np.exp(1j * np.angle(np.fft.fft2(rng.standard_normal(k.shape))))
    field = np.real(np.fft.ifft2(amp * phase))
    kk = np.arange(B_FIT[0], B_FIT[1] + 1)
    Ea, Es = shell_spectrum(field), shell_sums(field)
    assert abs(b_fit(Ea) + 4.0) < 0.02, b_fit(Ea)
    assert not [int(x) for x in kk[:-1] if Ea[x + 1] > B_MONOTONE * Ea[x]]
    assert [int(x) for x in kk[:-1] if Es[x + 1] > B_MONOTONE * Es[x]], "shell sums expected to rise somewhere"
    kk = np.arange(B_FIT[0], B_FIT[1] + 1)
    fake = np.zeros(129)
    fake[kk] = 3.0 * kk ** -4.1
    assert abs(b_fit(fake) + 4.1) < 1e-12


@pytest.mark.slow
@pytest.mark.skipif(os.environ.get("KH_CP7FB") != "1",
                    reason="CP-7F-B runs for hours per target; set KH_CP7FB=1 (see module docstring)")
@pytest.mark.parametrize("target", ["oracle"] + BACKEND_NAMES)
def test_cp7fb_enstrophy_range_slope(target):
    engine = _engine(target)
    r = b_run(engine)
    _write(f"cp7fb-{target}.json", dict({"target": target, "spec": "MASTER s3.2 E7"}, **r))
    line = (f"slope {r['slope']:.4f} (w1 {r['slope_w1']:.4f}, w2 {r['slope_w2']:.4f}, "
            f"stationarity {r['stationarity']:.4f}), monotone rises at k {r['monotone_rises_at_k']}")
    print(f"CP-7F-B {target}: {line}")
    assert r["samples"] == 20, f"{target}: {r['samples']} spectra sampled, expected 20"
    assert B_SLOPE[0] <= r["slope"] <= B_SLOPE[1], f"{target}: {line}; gate {B_SLOPE}"
    assert r["stationarity"] <= B_STATIONARITY, f"{target}: {line}; stationarity gate {B_STATIONARITY}"
    assert not r["monotone_rises_at_k"], f"{target}: {line}; monotone gate {B_MONOTONE}"
