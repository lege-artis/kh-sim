"""
test_canonical_physics_cp7fa.py -- CP-7F-A quantitative (Kelvin-Helmholtz growth rate
against Michalke 1964): ONE specification, run against the Python oracle and all five
backends through the request API.

Specification: MASTER-BRIEF.md s3.1, amendment E9 (2026-10-11), under ruling 2
(2026-10-07: the per-language CP-7F suites give way to one shared harness, as C did with
E6 and B, D with E7). The oracle's method of 2026-09-15
(development record, CP-7F-A assessment 2026-09-15), in the box and with the
seed chosen by the three-box study of 2026-10-11 (assessment s7, s8):

  Box Lx = pi^2, Ly = 2 pi, 256^2, dt = 0.001, velocity_shear 1, delta = Ly/40, so
  k.delta = 0.1 km. Base flow: the periodic-consistent DOUBLE tanh layer
  U = tanh((y - Ly/4)/delta) - tanh((y - 3Ly/4)/delta) - 1 (a single tanh layer in a
  periodic box carries a return shear -2/Ly that lowers the growth rate by ~9 %).
  Why delta = Ly/40 and not Ly/20 (2026-09-15): the two layers of the periodic box couple.
  At delta = Ly/20 (layers 10 delta apart) the linear spectrum at k.delta 0.4 has two
  growing modes, 0.1943 and 0.1815, +-3 % around Michalke's 0.1881 - not the -2e-4
  "image correction" assumed on 2026-09-15. At delta = Ly/40 (20 delta apart) they are
  0.1883 and 0.1881 (0.1875 / 0.1875 at k.delta 0.5): Michalke's single layer.
  Seed: the box's own leading inviscid mode (Fourier collocation of the Rayleigh operator
  on the grid's ny points), 1e-4 x Re(w(y) exp(i k x)). The single-layer Chebyshev seed
  of 2026-09-15 (its real part times cos k x) projects onto more than one mode and cost
  ~1e-3 in sigma_inf; with the box's own mode the measurement lands on linear theory.
  Twenty requests of 500 steps (0.5 s = 3.2 delta/U, as 1 000 steps were at delta = Ly/20),
  each continuing the returned field; after each, amp = max_ky |omega_hat(km, ky)|, the
  instantaneous rate ln(amp_n / amp_n-1) / 0.5 s, and the growth rate is the median of the
  positive rates within 30 % of their median (the adaptive window). sigma_nd = sigma*delta.
  Re in {4000, 10000, 20000, 40000, 80000} (Re_delta as {2000 .. 40000} at delta = Ly/20)
  at km 4 and 5 (k.delta 0.4, 0.5); per km a least-squares fit sigma_nd = a/Re + sigma_inf.
  The 1/Re term carries the viscous correction and the diffusion of the base flow during
  the window (the per-chunk rate falls by ~15 % over a run at Re 4000).
  Assert: max over km of sigma_inf in [0.1831, 0.1931] = Michalke's sigma(k.delta 0.4)
  = 0.1881 (Rayleigh solver on Michalke's problem; the box cannot sample 0.4446) +- 5e-3,
  Michalke's precision; and R^2 of each fit >= 0.99 (the 1/Re model holds).
  Study of 2026-10-11 (oracle, km 4, sigma_inf / R^2): this box and seed 0.18818 / 0.99962
  (own leading mode 0.18830); the same box with the Chebyshev seed 0.18719; Ly = 4 pi,
  256 x 512, delta = pi/10 (the same problem at twice the y resolution, Chebyshev seed)
  0.18746 / 0.99877; delta = Ly/20 with its own leading mode 0.19435 / 0.99997 against
  that mode's 0.19432.

Before E9 the five native suites ran their own A (a single tanh layer, sigma from two
amplitudes at t = 10 and 30, gate [0.180, 0.200]) and the oracle's A quantitative lived
in test_canonical_physics_python.py; the oracle's A qualitative (the stable band) stays
there.

Cost: 10 measurements x 10 000 steps at 256^2 per target. On the station 2026-10-11 one
step at 256^2 took 0.094 s (oracle) .. 0.33 s (Pascal), so one measurement is 16 min ..
55 min and a target 2.6 .. 9 h serially. Measurements are independent: a target given
several URLs (KH_CP7FA_URLS_<NAME>, e.g. KH_CP7FA_URLS_KH_PASCAL=http://localhost:8605,
http://localhost:8615) runs one measurement per URL at a time; the oracle runs
KH_CP7FA_ORACLE_PROCS processes. Run only with KH_CP7FA=1 (RUN-CP7FA.ps1); resumable
with KH_CP7F_CACHE=<dir> (one file per target, km and Re, written when the measurement
completes). KH_CP7F_OUT=<dir> receives one JSON record per target.
"""
import concurrent.futures as cf
import json
import math
import multiprocessing
import os
import pathlib
import sys
import threading
import time
import warnings

import numpy as np
import pytest
import requests

from conftest import BACKEND_NAMES, BACKEND_URLS

_ROOT = pathlib.Path(__file__).resolve().parents[2]
sys.path.insert(0, str(_ROOT / "shared" / "physics"))
import kh_physics  # noqa: E402
import functools  # noqa: E402

TIMEOUT = 7200

A_NX = A_NY = 256
A_LX = math.pi ** 2               # k.delta = 0.1 km with delta = Ly/40
A_LY = 2.0 * math.pi
A_DELTA = A_LY / 40.0             # the two layers 20 delta apart (see the module docstring)
A_DT = 0.001
A_CHUNK = 500                     # 0.5 s = 3.2 delta/U per request
A_CHUNKS = 20
A_AMP = 1e-4
A_RE = (4000.0, 10000.0, 20000.0, 40000.0, 80000.0)
A_KM = (4, 5)
A_ANCHOR = 0.1881                         # Michalke sigma(k.delta = 0.4), Rayleigh solver
A_GATE = (A_ANCHOR - 5e-3, A_ANCHOR + 5e-3)
A_REF_BY_KM = {4: 0.1881, 5: 0.1875}
A_R2_MIN = 0.99


# --------------------------------------------------------------------------- engines
class _Oracle:
    name = "oracle"

    def advance(self, req: dict) -> dict:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            res = kh_physics.simulate(kh_physics.SimulationRequest(**req))
        return {"vorticity": res.vorticity, "t_final": res.t_final,
                "steps_completed": res.steps_completed}


class _Http:
    # A plain requests.post per call, no Session (Pascal resets a kept-alive connection).
    def __init__(self, name, url):
        self.name, self.url = name, url

    def advance(self, req: dict) -> dict:
        r = requests.post(f"{self.url}/simulate", json=req, timeout=TIMEOUT)
        assert r.status_code == 200, f"{self.name}: HTTP {r.status_code}: {r.text[:300]}"
        return r.json()


def _urls(target: str) -> list:
    """The URLs that serve `target`: KH_CP7FA_URLS_<NAME> (comma-separated) or the one
    of conftest. Each must answer /health."""
    env = os.environ.get("KH_CP7FA_URLS_" + target.upper().replace("-", "_"), "")
    urls = [u.strip().rstrip("/") for u in env.split(",") if u.strip()] or [BACKEND_URLS[target]]
    for u in urls:
        try:
            requests.get(f"{u}/health", timeout=5).raise_for_status()
        except requests.RequestException as exc:
            if os.environ.get("KH_GOLDEN_REQUIRE_SERVER") == "1":
                pytest.fail(f"{target} not reachable at {u}: {exc}")
            pytest.skip(f"{target} not reachable at {u}")
    return urls


# --------------------------------------------------------------------------- method
def base_vorticity_double(ny: int = A_NY) -> np.ndarray:
    """omega = -dU/dy of U = tanh((y-Ly/4)/delta) - tanh((y-3Ly/4)/delta) - 1."""
    y = np.linspace(0, A_LY, ny, endpoint=False)
    s1 = 1.0 / np.cosh((y - A_LY / 4.0) / A_DELTA) ** 2
    s2 = 1.0 / np.cosh((y - 3.0 * A_LY / 4.0) / A_DELTA) ** 2
    return -1.0 / A_DELTA * (s1 - s2)            # the oracle's expression, bit for bit


@functools.lru_cache(maxsize=None)
def leading_mode(km: int, ny: int = A_NY, n_modes: int = 1):
    """Inviscid linear modes of the periodic double layer on ny points (Fourier collocation
    of the Rayleigh operator (U - c)(D^2 - k^2) psi = U'' psi, written as
    d/dt Lap psi = -ik U Lap psi + ik U'' psi): the n_modes fastest growth rates x delta,
    and the leading vorticity eigenfunction w = Lap psi, max |w| = 1 and real there."""
    k = km * 2.0 * math.pi / A_LX
    y = np.arange(ny) * A_LY / ny
    U = np.tanh((y - A_LY / 4.0) / A_DELTA) - np.tanh((y - 3.0 * A_LY / 4.0) / A_DELTA) - 1.0
    ky = np.fft.fftfreq(ny, 1.0 / ny) * 2.0 * math.pi / A_LY
    D2 = np.fft.ifft(np.diag(-(ky ** 2)) @ np.fft.fft(np.eye(ny), axis=0), axis=0)
    lap = D2 - k * k * np.eye(ny)
    A = -1j * k * np.diag(U) @ lap + 1j * k * np.diag(np.real(D2 @ U))
    lam, vec = np.linalg.eig(np.linalg.solve(lap, A))
    order = np.argsort(-lam.real)
    w = lap @ vec[:, order[0]]
    w = w / w[np.argmax(np.abs(w))]                 # max |w| = 1, real and positive there
    return tuple(float(lam[i].real) * A_DELTA for i in order[:n_modes]), w


def initial_field(km: int, nx: int = A_NX, ny: int = A_NY) -> np.ndarray:
    """Base vorticity + A_AMP x Re(w(y) exp(i k x)), w the box's leading inviscid mode."""
    k = km * 2.0 * math.pi / A_LX
    _, w = leading_mode(km, ny)
    x = np.linspace(0, A_LX, nx, endpoint=False)
    return np.tile(base_vorticity_double(ny), (nx, 1)) + A_AMP * np.real(w[np.newaxis, :] * np.exp(1j * k * x)[:, np.newaxis])


def adaptive_rate(rates: list) -> float:
    """Median of the positive instantaneous rates within 30 % of their median."""
    arr = np.asarray(rates, dtype=float)
    pos = arr[arr > 0]
    if len(pos) == 0:
        return float("nan")
    med = np.median(pos)
    stable = pos[np.abs(pos - med) / med < 0.30]
    return float(np.median(stable)) if len(stable) else float(np.median(pos))


def measure(engine, km: int, re: float, log=print, n: int = A_NX, chunks: int = A_CHUNKS) -> dict:
    """One sigma_nd measurement: `chunks` requests of A_CHUNK steps from initial_field
    (n x n; the specification is A_NX, A_CHUNKS)."""
    omega = initial_field(km, nx=n, ny=n)
    rates, amps, t, t0 = [], [], 0.0, time.time()
    amp_prev = None
    for c in range(chunks):
        req = {"grid_nx": n, "grid_ny": n, "domain_lx": A_LX, "domain_ly": A_LY,
               "dt": A_DT, "steps": A_CHUNK, "reynolds_number": re, "velocity_shear": 1.0,
               "perturbation_amplitude": 0.0, "perturbation_mode": 1, "t_start": t,
               "initial_omega": [float(v) for v in omega.ravel()]}
        body = engine.advance(req)
        assert body.get("steps_completed") == A_CHUNK, \
            f"{engine.name}: steps_completed {body.get('steps_completed')} != {A_CHUNK}"
        omega = np.asarray(body["vorticity"], dtype=float).reshape(n, n)
        assert np.all(np.isfinite(omega)), f"{engine.name}: non-finite vorticity (km {km}, Re {re:g})"
        t = float(body["t_final"])
        amp = float(np.abs(np.fft.rfft2(omega)[km, :]).max())
        amps.append(amp)
        if amp_prev is not None and amp > 1e-15 and amp_prev > 1e-15:
            rates.append(math.log(amp / amp_prev) / (A_CHUNK * A_DT))
        amp_prev = amp
    sigma = adaptive_rate(rates)
    out = {"km": km, "Re": re, "sigma_nd": sigma * A_DELTA, "rates": rates, "amps": amps,
           "wall_s": round(time.time() - t0, 1)}
    log(f"CP-7F-A {engine.name}: km {km} Re {re:g}: sigma_nd {out['sigma_nd']:.6f} ({out['wall_s']} s)")
    return out


def extrapolate(rows: list) -> dict:
    """Per km: least-squares sigma_nd = a / Re + sigma_inf, its R^2 and intercept s.e."""
    per_km = {}
    for km in A_KM:
        pts = [(1.0 / r["Re"], r["sigma_nd"]) for r in rows if r["km"] == km and np.isfinite(r["sigma_nd"])]
        x = np.array([p[0] for p in pts]); y = np.array([p[1] for p in pts])
        if len(x) < 3:
            per_km[km] = {"n": len(x), "sigma_inf": float("nan"), "r2": float("nan"), "se": float("nan")}
            continue
        a, b = np.polyfit(x, y, 1)
        res = y - (a * x + b)
        ss_res = float(np.sum(res ** 2)); ss_tot = float(np.sum((y - y.mean()) ** 2))
        sxx = float(np.sum((x - x.mean()) ** 2))
        se = math.sqrt(ss_res / (len(x) - 2) * (1.0 / len(x) + x.mean() ** 2 / sxx)) if len(x) > 2 and sxx > 0 else float("nan")
        per_km[km] = {"n": len(x), "slope": float(a), "sigma_inf": float(b),
                      "r2": 1.0 - ss_res / ss_tot if ss_tot > 0 else float("nan"), "se": se,
                      "reference": A_REF_BY_KM[km]}
    return per_km


def _cache_file(target: str, km: int, re: float):
    d = os.environ.get("KH_CP7F_CACHE")
    return pathlib.Path(d) / f"cp7fa-{target}-km{km}-Re{int(re)}.json" if d else None


def _oracle_job(km: int, re: float, n: int = A_NX, chunks: int = A_CHUNKS) -> dict:
    """One oracle measurement in a worker process (top level: picklable)."""
    return measure(_Oracle(), km, re, log=lambda *_: None, n=n, chunks=chunks)


def a_run(target: str, log=print) -> list:
    """All A_KM x A_RE measurements of one target, cached and resumed per measurement."""
    jobs = [(km, re) for km in A_KM for re in A_RE]
    rows, todo = {}, []
    for km, re in jobs:
        f = _cache_file(target, km, re)
        if f and f.exists():
            rows[(km, re)] = json.loads(f.read_text())
            log(f"CP-7F-A {target}: km {km} Re {re:g} from cache: sigma_nd {rows[(km, re)]['sigma_nd']:.6f}")
        else:
            todo.append((km, re))

    def keep(km, re, r):
        rows[(km, re)] = r
        f = _cache_file(target, km, re)
        if f:
            f.parent.mkdir(parents=True, exist_ok=True)
            f.write_text(json.dumps(r))

    if todo and target == "oracle":
        procs = max(1, int(os.environ.get("KH_CP7FA_ORACLE_PROCS", "1")))
        if procs == 1:
            for km, re in todo:
                keep(km, re, measure(_Oracle(), km, re, log))
        else:
            # spawn on every OS: what Windows does anyway, so the station path is the tested one
            with cf.ProcessPoolExecutor(procs, mp_context=multiprocessing.get_context("spawn")) as ex:
                futs = {ex.submit(_oracle_job, km, re): (km, re) for km, re in todo}
                for fu in cf.as_completed(futs):
                    km, re = futs[fu]; r = fu.result(); keep(km, re, r)
                    log(f"CP-7F-A oracle: km {km} Re {re:g}: sigma_nd {r['sigma_nd']:.6f} ({r['wall_s']} s)")
    elif todo:
        urls = _urls(target)
        lock = threading.Lock()
        queue = list(todo)

        def worker(url):
            engine = _Http(target, url)
            while True:
                with lock:
                    if not queue:
                        return
                    km, re = queue.pop(0)
                r = measure(engine, km, re, log)
                with lock:
                    keep(km, re, r)

        with cf.ThreadPoolExecutor(len(urls)) as ex:
            for fu in [ex.submit(worker, u) for u in urls]:
                fu.result()
    return [rows[j] for j in jobs]


def _write(name: str, record: dict) -> None:
    out = os.environ.get("KH_CP7F_OUT")
    if out:
        pathlib.Path(out).mkdir(parents=True, exist_ok=True)
        (pathlib.Path(out) / name).write_text(json.dumps(record, indent=1))


# --------------------------------------------------------------------------- tests
def test_cp7fa_helpers():
    """The double layer has (to its non-periodic tails) zero mean vorticity and the tanh
    profile's velocity jump (2 per layer, U0 = 1); k.delta = 0.1 km; the adaptive window
    drops a transient and a saturated chunk; the 1/Re fit recovers an exact law."""
    w = base_vorticity_double()
    assert abs(w.mean()) < 1e-5 * abs(w).max()          # the sech^2 tails are not periodic
    dy = A_LY / A_NY
    jump = -np.cumsum(w) * dy                            # U(y) - U(0), U = -int omega dy
    assert abs((jump.max() - jump.min()) - 2.0) < 1e-3, jump.max() - jump.min()
    assert abs(4 * 2.0 * math.pi / A_LX * A_DELTA - 0.4) < 1e-15        # k.delta = 0.1 km
    assert abs(adaptive_rate([0.05, 0.6, 0.61, 0.59, 0.6, 0.2]) - 0.6) < 0.006
    rows = [{"km": km, "Re": re, "sigma_nd": 0.18 - 3.0 / re} for km in A_KM for re in A_RE]
    ex = extrapolate(rows)
    assert all(abs(ex[km]["sigma_inf"] - 0.18) < 1e-12 and ex[km]["r2"] > 0.999999 for km in A_KM)


def test_cp7fa_box_is_michalke():
    """The box's two layers do not couple: its two leading inviscid modes are Michalke's
    single-layer rate (0.1881 at k.delta 0.4, 0.1875 at 0.5; the box cannot sample the
    maximum at 0.4446) to 3e-4, on 128 points as on 256 (2026-10-11: 0.1883 / 0.1881,
    0.1875 / 0.1875). At delta = Ly/20 the same code gives 0.1943 / 0.1815: the coupling
    that moved the box to Ly/40. The seed is a growing mode of this operator."""
    for km, ref in ((4, 0.1881), (5, 0.1875)):
        rates, w = leading_mode(km, 128, n_modes=2)
        assert all(abs(r - ref) < 3e-4 for r in rates), (km, rates)
        assert abs(np.abs(w).max() - 1.0) < 1e-12


def test_cp7fa_window_equals_the_oracle_code():
    """The growth-rate window equals the oracle code of 2026-09-15
    (test_canonical_physics_python._adaptive_sigma_from_chunks) bit for bit."""
    import test_canonical_physics_python as legacy
    rng = np.random.default_rng(7)
    for _ in range(20):
        rates = list(rng.normal(0.6, 0.2, 19))
        assert adaptive_rate(rates) == legacy._adaptive_sigma_from_chunks(rates)


@pytest.mark.parametrize("target", BACKEND_NAMES)
def test_cp7fa_request_path_equals_the_oracle(target):
    """Smoke for the long run: one measurement at 32^2, three chunks, through the
    backend's request API equals the oracle's to 1e-9 (measured 2026-10-11 on Rust, C++,
    Fortran, Pascal: 1e-14 .. 2e-12; the flow is laminar, so the backends agree with the
    oracle far below the gate, unlike CP-7F-B)."""
    url = _urls(target)[0]
    quiet = lambda *_: None  # noqa: E731
    mine = measure(_Http(target, url), 4, 2000.0, log=quiet, n=32, chunks=3)["sigma_nd"]
    ref = measure(_Oracle(), 4, 2000.0, log=quiet, n=32, chunks=3)["sigma_nd"]
    assert abs(mine - ref) <= 1e-9 * abs(ref), f"{target}: sigma_nd {mine!r} vs oracle {ref!r}"


@pytest.mark.slow
@pytest.mark.skipif(os.environ.get("KH_CP7FA") != "1",
                    reason="CP-7F-A runs for hours per target; set KH_CP7FA=1 (see module docstring)")
@pytest.mark.parametrize("target", ["oracle"] + BACKEND_NAMES)
def test_cp7fa_dispersion_michalke(target):
    rows = a_run(target)
    ex = extrapolate(rows)
    best = max(ex[km]["sigma_inf"] for km in A_KM)
    _write(f"cp7fa-{target}.json", {"target": target, "spec": "MASTER s3.1 E9", "gate": list(A_GATE),
                                    "sigma_inf_max": best, "per_km": {str(k): v for k, v in ex.items()},
                                    "measurements": rows})
    line = ", ".join(f"km {km}: sigma_inf {ex[km]['sigma_inf']:.5f} +- {ex[km]['se']:.1e} "
                     f"(R^2 {ex[km]['r2']:.5f}, Michalke {A_REF_BY_KM[km]})" for km in A_KM)
    print(f"CP-7F-A {target}: {line}; max {best:.5f}, gate [{A_GATE[0]:.4f}, {A_GATE[1]:.4f}]")
    for km in A_KM:
        assert ex[km]["n"] == len(A_RE), f"{target}: km {km}: {ex[km]['n']} finite rates of {len(A_RE)}"
        assert ex[km]["r2"] >= A_R2_MIN, f"{target}: km {km}: 1/Re fit R^2 {ex[km]['r2']:.4f} < {A_R2_MIN}; {line}"
    assert A_GATE[0] <= best <= A_GATE[1], f"{target}: {line}; max {best:.5f} outside {A_GATE}"
