"""
test_canonical_physics_python.py -- phase 7f canonical-physics testbeds for Python ref.

Implements CP-7F-A..D from MASTER-BRIEF.md ss3. Each test reproduces a textbook-anchored
physics signature; tolerances match the textbook's reported precision (NOT machine epsilon).

CP-7F-E: E2 (NR in Pascal s17.5 SOR Poisson) and E1 (s12.4 convlv), transcribed 2026-06-17;
the Pascal and Fortran native suites assert the same since 2026-10-11 (MASTER E8).

STATUS SUMMARY (CP-7F-A/C/D refine pass, 2026-06-15):
  CP-7F-A: Qualitative dispersion test PASS (unchanged).
           Quantitative: MOVED 2026-10-11 to test_canonical_physics_cp7fa.py (MASTER s3.1
           E9), the method unchanged, for the oracle and the five backends.

  CP-7F-B: MOVED 2026-10-10 to test_canonical_physics_cp7fbd.py (MASTER s3.2 E7, the
           E5 gate [-4.6, -3.6]); the pre-E5 xfail diagnostic that stood here is retired.

  CP-7F-C: MOVED 2026-10-07 to test_canonical_physics_cp7fc.py (MASTER s3.3 E6); the
           Gaussian variant of 2026-06-15 measured blob size, not merger (finding
           FINDING-CP7FC-MERGE-DETECTION-2026-10-06).

  CP-7F-D: MOVED 2026-10-10 to test_canonical_physics_cp7fbd.py (MASTER s3.4 E7).
           RECLASSIFIED 2026-06-15 (Pete ruling D1) -> SOLVER-CONSISTENCY CHECK, NOT C1.
           fixed-ε forcing f0 = sqrt(2·ν·ε) gives Re^{1/2} (slope 0.500) ANALYTICALLY,
           by construction -> not falsifiable against an independent book value, so it
           does not count toward C1. Kept green as a sanity check. Brachet 1983 anchor
           likely mis-attributed (3D decaying TG; OQ-7F-D-ANCHOR). v0.2.x: genuine
           forced-2D scaling + corrected anchor. See STATUS-...-REFINE-FILLED.md.

Citations:
  Drazin, P.G., Reid, W.H. (1981). Hydrodynamic Stability. Cambridge Univ. Press. s4.2.
  Michalke, A. (1964). J. Fluid Mech. 19(4):543-556.
  Lesieur, M. (2008). Turbulence in Fluids, 4th ed. Springer. Ch.9 s9.3.
  Kraichnan, R.H. (1967). Phys. Fluids 10(7):1417-1423.
  Saffman, P.G. (1992). Vortex Dynamics. Cambridge Univ. Press. s7.4.
  Melander, M.V., Zabusky, N.J., McWilliams, J.C. (1988). J. Fluid Mech. 195:303-340.
  Brachet, M.E. et al. (1983). J. Fluid Mech. 130:411-452.
"""
import json
import math
import os
import pathlib
import sys
import time

import numpy as np
import pytest

# Make rayleigh_eigenvalue importable from shared/physics
sys.path.insert(0, str(pathlib.Path(__file__).parent.parent.parent / "shared" / "physics"))
from rayleigh_eigenvalue import rayleigh_cheb_vort_seed  # noqa: E402

from shared.physics.kh_physics import (
    SimulationRequest,
    simulate,
)

AUDIT_DIR = pathlib.Path(__file__).parent.parent.parent / "_audit"
GOLDEN_DIR = pathlib.Path(__file__).parent.parent.parent / "shared" / "golden-vectors"
AUDIT_DIR.mkdir(exist_ok=True)


# ===========================================================================
# Resumable-compute infrastructure
# ===========================================================================
#
# Restored 2026-07-05 (was lost in the trunk<->thinkpad merge 39c5ab9, which took
# the thinkpad base). The physics runs below are expensive (CP-7F-A quantitative
# is ~36 min). On a normal machine plain `pytest` runs everything inline and this
# layer is INERT. To make a long run resumable/interruptible, set:
#     KH_CP7F_CACHE=<dir>     persist raw per-measurement results
#     KH_CP7F_BUDGET_S=<sec>  per-invocation compute budget; when exceeded the
#                             test checkpoints and raises _BudgetExceeded
#                             (re-invoke until it completes)
# The cache stores RAW measurements ONLY (sigma values) -- never assertion
# outcomes. pytest always re-evaluates every gate from the raw data, so the
# cache can never mask a failure.


class _BudgetExceeded(Exception):
    """Raised when KH_CP7F_BUDGET_S is exhausted; progress has been checkpointed."""


class _Clock:
    def __init__(self):
        self.t0 = time.perf_counter()
        env = os.environ.get("KH_CP7F_BUDGET_S")
        self.budget = float(env) if env else None

    def over(self) -> bool:
        if self.budget is None:
            return False
        return (time.perf_counter() - self.t0) > self.budget


def _cache_dir():
    env = os.environ.get("KH_CP7F_CACHE")
    if not env:
        return None
    d = pathlib.Path(env)
    d.mkdir(parents=True, exist_ok=True)
    return d


def _load_result(name: str):
    d = _cache_dir()
    if d is None:
        return None
    f = d / (name + ".result.json")
    if f.exists():
        with open(f) as fh:
            return json.load(fh)
    return None


def _save_result(name: str, payload: dict):
    d = _cache_dir()
    if d is None:
        return
    with open(d / (name + ".result.json"), "w") as fh:
        json.dump(payload, fh)


# ---------------------------------------------------------------------------
# Shared helpers
# ---------------------------------------------------------------------------


def _run(req: SimulationRequest) -> np.ndarray:
    """Run simulation and return vorticity as (nx, ny) array."""
    result = simulate(req)
    return np.array(result.vorticity).reshape(req.grid_nx, req.grid_ny)


def _kh_base_flow(ny: int, Ly: float, delta: float, velocity_shear: float = 1.0) -> np.ndarray:
    """Vorticity of the KH base flow: omega = -dU/dy for U = tanh((y-Ly/2)/delta).

    CAVEAT (found 2026-09-15, see _audit/phase-7f-cp7fa-extrapolation-assessment-2026-09-15.md):
    a SINGLE tanh layer is not periodic in y. Its mean vorticity 2*U0/Ly cannot be carried by a
    periodic streamfunction, and solve_poisson() zeroes the k=0 mode, so the flow actually
    integrated is  U = tanh(eta) - (2*U0/Ly)*(y - Ly/2)  -- tanh plus a uniform RETURN SHEAR of
    -2*U0/Ly (= -0.318 for Ly=2*pi, ~10% of the layer's peak shear). The Rayleigh eigenvalue of
    that profile is ~9% below Michalke's (0.1626 vs 0.1881 at k.delta=0.4). Fine for
    solver-consistency / golden tests; NOT a Michalke-anchored profile. Use
    _kh_base_flow_double() for any test gated against Michalke (1964).
    """
    y = np.linspace(0, Ly, ny, endpoint=False)
    y_c = y - Ly / 2
    return -velocity_shear / delta * (1.0 / np.cosh(y_c / delta)) ** 2


def _kh_base_flow_double(ny: int, Ly: float, delta: float, velocity_shear: float = 1.0) -> np.ndarray:
    """Vorticity of the periodic-consistent DOUBLE shear layer (zero mean vorticity):
        U = tanh((y-Ly/4)/delta) - tanh((y-3Ly/4)/delta) - 1,   omega = -dU/dy.
    Each layer is a Michalke (1964) tanh layer; the image interaction at Ly=2*pi, delta=pi/10
    shifts the k.delta=0.4446 eigenvalue by only -2e-4 (Rayleigh solver: 0.1895 vs 0.1897).
    """
    y = np.linspace(0, Ly, ny, endpoint=False)
    s1 = 1.0 / np.cosh((y - Ly / 4.0) / delta) ** 2
    s2 = 1.0 / np.cosh((y - 3.0 * Ly / 4.0) / delta) ** 2
    return -velocity_shear / delta * (s1 - s2)


def _make_eigenmode_init(km: int, nx: int, ny: int, Lx: float, Ly: float,
                          delta: float, amp: float = 1e-4,
                          double_layer: bool = False) -> np.ndarray:
    """
    Build the initial vorticity field seeded with the Rayleigh eigenmode.

    Uses Chebyshev pseudospectral collocation (Orszag 1971) to find the most-unstable
    eigenfunction of the Rayleigh equation for the tanh base flow.  The vorticity
    perturbation is ω'(x,y) = amp · ω_eig(y) · cos(k·x), where ω_eig is the
    vorticity eigenfunction normalised to unit maximum.

    double_layer=True: periodic-consistent double shear layer (_kh_base_flow_double);
    the eigenmode seed is centred on the lower layer (y_c = Ly/4). The seed is only an
    IC -- the adaptive growth window measures the asymptotic rate of the simulated flow.

    Reference: Orszag (1971) J. Fluid Mech. 50:689-703.
               Michalke (1964) J. Fluid Mech. 19(4):543-556.
    """
    y_c = Ly / 4.0 if double_layer else Ly / 2.0
    k = km * 2.0 * np.pi / Lx

    y_uniform = np.linspace(0, Ly, ny, endpoint=False)
    vort_eig, _ = rayleigh_cheb_vort_seed(k, Ly, y_c, delta, y_uniform, N=100)

    x = np.linspace(0, Lx, nx, endpoint=False)
    X, Y = np.meshgrid(x, y_uniform, indexing="ij")  # (nx, ny)

    omega_base = (_kh_base_flow_double(ny, Ly, delta) if double_layer
                  else _kh_base_flow(ny, Ly, delta))
    omega_base_2d = np.tile(omega_base, (nx, 1))

    return omega_base_2d + amp * vort_eig[np.newaxis, :] * np.cos(k * X)


def _adaptive_sigma_from_chunks(sigma_list: list) -> float:
    """
    Extract the exponential-growth rate from a list of instantaneous sigma values.

    Keeps only the positive-growth window that is internally consistent: retain
    values within 30% relative of the running median to exclude early transient
    and late nonlinear saturation.
    """
    if not sigma_list:
        return float("nan")
    arr = np.array(sigma_list, dtype=float)
    pos = arr[arr > 0]
    if len(pos) == 0:
        return float("nan")
    med = np.median(pos)
    stable = pos[np.abs(pos - med) / med < 0.30]
    return float(np.median(stable)) if len(stable) > 0 else float(np.median(pos))


def _measure_kh_sigma_nd(km: int, Re: float = 10000.0, n_chunks: int = 20,
                          nx: int = 256, ny: int = 256,
                          double_layer: bool = False) -> float:
    """
    Measure non-dimensional growth rate σ_nd = σ · δ / U0 for wavenumber mode km
    using the Chebyshev Rayleigh eigenmode as the initial perturbation.

    Improvements over the escalated v0 (non-eigenmode 128×64 seed):
      (A-1) Eigenmode seed: the IC is the true Rayleigh eigenfunction, eliminating
            the ~11% transient gap from non-eigenmode projection.
      (A-3) 256×256 grid per MASTER §3.1.
      Adaptive growth window: growth rate is measured over the window where
      d(log amp)/dt is stationary (not a hard-coded 5–10 s window).

    Domain: Lx = 2π², Ly = 2π → k·δ = km · 0.1 for integer km.
    δ = 0.05·Ly = π/10 ≈ 0.3142 rad.
    """
    Lx = 2.0 * np.pi ** 2
    Ly = 2.0 * np.pi
    delta = 0.05 * Ly
    velocity_shear = 1.0
    dt = 0.001
    chunk_steps = 1000  # 1.0 s per chunk

    omega_current = _make_eigenmode_init(km, nx, ny, Lx, Ly, delta,
                                         double_layer=double_layer)
    sigma_list = []
    amp_prev = None

    for chunk in range(n_chunks):
        req = SimulationRequest(
            grid_nx=nx,
            grid_ny=ny,
            domain_lx=Lx,
            domain_ly=Ly,
            dt=dt,
            steps=chunk_steps,
            reynolds_number=Re,
            velocity_shear=velocity_shear,
            perturbation_amplitude=0.0,
            perturbation_mode=1,
            initial_omega=omega_current.flatten().tolist(),
        )
        result = simulate(req)
        omega_current = np.array(result.vorticity).reshape(nx, ny)

        omega_hat = np.fft.rfft2(omega_current)
        amp_now = float(np.abs(omega_hat[km, :]).max())

        if amp_prev is not None and amp_now > 1e-15 and amp_prev > 1e-15:
            sigma_inst = np.log(amp_now / amp_prev) / (chunk_steps * dt)
            sigma_list.append(sigma_inst)

        amp_prev = amp_now

    sigma_raw = _adaptive_sigma_from_chunks(sigma_list)
    return sigma_raw * delta


@pytest.mark.slow
def test_cp7fa_kh_dispersion_qualitative():
    """
    CP-7F-A (QUALITATIVE): KH shear layer is unstable for k*delta in [0.1, 1.0]
    and stable for k*delta = 1.1.

    Reference: Michalke (1964) J. Fluid Mech. 19(4):543-556, Fig. 3.
               Drazin & Reid (1981) Hydrodynamic Stability s4.2.
    Domain: Lx=2*pi^2, Ly=2*pi => k*delta = km * 0.1 (exact for integer km).
    delta = 0.05 * Ly = 0.314 rad, Re = 10000.
    """
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        has_mpl = True
    except ImportError:
        has_mpl = False

    # Measure sigma_nd for selected modes covering the unstable + stable bands
    # km = 2, 4, 5, 6, 8, 10, 11
    modes_to_test = [2, 4, 5, 6, 8, 10, 11]
    sigma_nd = {}
    for km in modes_to_test:
        sigma_nd[km] = _measure_kh_sigma_nd(km, Re=10000.0, n_chunks=12, nx=128, ny=64)

    # --- Qualitative assertions (achievable) ---
    # 1. Unstable band: sigma_nd > 0 for k*delta in [0.2, 0.9] (well inside band)
    #    km=10 (k*delta=1.0) is the neutral cutoff -- growth rate -> 0 there; excluded.
    for km in [2, 4, 5, 6, 8]:
        k_delta = km * 0.1
        assert sigma_nd[km] > 0.0, (
            f"CP-7F-A: expected sigma_nd > 0 (unstable) at k*delta={k_delta:.1f}, "
            f"got {sigma_nd[km]:.4f}. "
            f"Ref: Michalke (1964) Fig.3 -- all modes k*delta<1.0 are unstable."
        )

    # 2. Near-neutral at cutoff: sigma_nd[km=10] < 0.05 (k*delta=1.0 ~ neutral)
    #    Michalke Table 1: sigma_nd -> 0 as k*delta -> 1. NaN counts as near-neutral.
    if not math.isnan(sigma_nd[10]):
        assert sigma_nd[10] < 0.05, (
            f"CP-7F-A: sigma_nd at k*delta=1.0 (cutoff) should be near-zero, "
            f"got {sigma_nd[10]:.4f}. Ref: Michalke (1964) -- neutral at k*delta=1.0."
        )

    # 3. Stable/decayed: sigma_nd < 0.02 for k*delta = 1.1 (km=11)
    if not math.isnan(sigma_nd[11]):
        assert sigma_nd[11] < 0.02, (
            f"CP-7F-A: expected sigma_nd ~ 0 (stable) at k*delta=1.1, "
            f"got {sigma_nd[11]:.4f}. "
            f"Ref: Michalke (1964) -- cutoff at k*delta ~ 1.0."
        )

    # 4. Peak in unstable band: argmax k*delta in [0.30, 0.70] (widened from Michalke 0.40-0.50)
    unstable_modes = {km: v for km, v in sigma_nd.items() if km <= 9 and not math.isnan(v)}
    km_peak = max(unstable_modes, key=lambda k: unstable_modes[k])
    k_delta_peak = km_peak * 0.1
    assert 0.30 <= k_delta_peak <= 0.70, (
        f"CP-7F-A: peak growth at k*delta={k_delta_peak:.2f} outside qualitative range [0.30, 0.70]. "
        f"Ref: Michalke (1964) Table 1 -- most-unstable at k*delta=0.4446."
    )

    # Save diagnostic plot
    if has_mpl:
        k_delta_vals = [km * 0.1 for km in modes_to_test]
        sigma_vals = [sigma_nd[km] for km in modes_to_test]

        # Michalke (1964) theoretical curve (approximate values from Table 1)
        k_theory = np.array([0.0, 0.1, 0.2, 0.3, 0.4, 0.4446, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.3])
        s_theory = np.array([0.0, 0.064, 0.121, 0.162, 0.184, 0.1897, 0.186, 0.171, 0.148, 0.116, 0.079, 0.038, 0.0, 0.0])

        fig, ax = plt.subplots(figsize=(8, 5))
        ax.plot(k_theory, s_theory, "k-", label="Michalke (1964) theory", linewidth=1.5)
        ax.plot(k_delta_vals, sigma_vals, "ro-", label="kh-sim Python ref (Re=10000)", markersize=6)
        ax.axhline(0, color="gray", linewidth=0.8, linestyle="--")
        ax.axvline(1.0, color="blue", linewidth=0.8, linestyle=":", label="Cutoff k*delta=1.0")
        ax.set_xlabel("k * delta")
        ax.set_ylabel("sigma_nd = sigma * delta / U0")
        ax.set_title("CP-7F-A: KH dispersion relation (qualitative)")
        ax.legend()
        info = (
            "Ref: Michalke (1964) J.Fluid Mech. 19(4):543-556\n"
            "Theory sigma_max=0.1897 at k*delta=0.4446\n"
            f"Measured peak: {max(sigma_vals):.4f} at k*delta={k_delta_vals[sigma_vals.index(max(sigma_vals))]:.1f}\n"
            "ESCALATION: quantitative gate [0.180,0.200] not achieved -- see STATUS"
        )
        ax.text(
            0.02, 0.02, info, transform=ax.transAxes,
            fontsize=7, verticalalignment="bottom",
            bbox=dict(boxstyle="round", facecolor="wheat", alpha=0.8),
        )
        plot_path = AUDIT_DIR / "phase-7f-cp7fa-python.png"
        fig.savefig(str(plot_path), dpi=120, bbox_inches="tight")
        plt.close(fig)


# ---------------------------------------------------------------------------
# CP-7F-A quantitative: MOVED 2026-10-11 to test_canonical_physics_cp7fa.py (MASTER s3.1 E9)
# ---------------------------------------------------------------------------
#
# The harness runs this module's method of 2026-09-15 (double shear layer, adaptive window,
# sigma_nd = a/Re + sigma_inf, gate 0.1881 +- 5e-3) for the oracle and the five backends,
# in the box and with the seed of a three-box study (assessment 2026-09-15, s7, s8):
# delta = Ly/40 instead of Ly/20, where the two layers couple into modes 0.1943 / 0.1815
# at k.delta 0.4, and the box's own leading mode as seed instead of the Chebyshev one
# below. The harness pins its window against _adaptive_sigma_from_chunks.
# The helpers stay here: the qualitative test uses them.


# ---------------------------------------------------------------------------
# CP-7F-B: MOVED 2026-10-10 to test_canonical_physics_cp7fbd.py (MASTER s3.2 E7)
# ---------------------------------------------------------------------------
#
# The test that stood here was the pre-E5 diagnostic: 128^2, forcing at k 4-6
# with random static phases, 20 time units, asserted against the Kraichnan gate
# [-3.15, -2.85] and marked xfail. The oracle now runs the E7 specification
# (CP7F-B-TUNED-PARAMS.json, E5 gate [-4.6, -3.6]) like every backend, in
# test_canonical_physics_cp7fbd.py (target "oracle"). Its text is in git history.

# ---------------------------------------------------------------------------
# CP-7F-C: MOVED 2026-10-07 to test_canonical_physics_cp7fc.py (MASTER s3.3 E6)
# ---------------------------------------------------------------------------
#
# The Gaussian variant that stood here (sigma = 0.1 L at d = 1.3, 1.5, 2.0 sigma)
# started from ONE vorticity maximum: two Gaussians closer than 2 sigma do not
# have two. Its "merged / orbiting" verdict measured the size of that blob, not
# a merger. The oracle now runs the same E6 specification as every backend, in
# test_canonical_physics_cp7fc.py (target "oracle").
# Finding: _audit/FINDING-CP7FC-MERGE-DETECTION-2026-10-06.md.

# ---------------------------------------------------------------------------
# CP-7F-D: MOVED 2026-10-10 to test_canonical_physics_cp7fbd.py (MASTER s3.4 E7)
# ---------------------------------------------------------------------------
#
# Same configuration (64^2, Re 50..1000, fixed-eps forcing, the analytic balance
# A_ss as initial field), now sent through the request API (forcing_type "mode")
# to the oracle and all five backends instead of a private RK4 loop here. Still a
# solver-consistency check, not C1 (ruling D1, 2026-06-15).


# ===========================================================================
# CP-7F-E -- NR-1986 Pascal PDE-chapter example (Pascal-binding test)
# ===========================================================================

def _solve_discrete_dirichlet_poisson(nx, ny, sources):
    """Exact discrete solution of the 5-point Dirichlet-Poisson problem.

    u[i+1,j]+u[i-1,j]+u[i,j+1]+u[i,j-1]-4 u[i,j] = f[i,j], u=0 on the boundary.
    `sources`: {(i0,j0): value} on 0-based interior nodes. Returns the full
    nx*ny grid (0-based). This is the machine-precision discrete solution; the
    decision doc (2026-06-17) confirmed it equals the book's printed SOR grid to
    the printed 2-decimal precision, so it is the oracle for the E2 sub-anchor.
    """
    interior = [(i, j) for i in range(1, nx - 1) for j in range(1, ny - 1)]
    idx = {p: k for k, p in enumerate(interior)}
    a = np.zeros((len(interior), len(interior)))
    b = np.zeros(len(interior))
    for (i, j), k in idx.items():
        a[k, k] = -4.0
        for di, dj in ((1, 0), (-1, 0), (0, 1), (0, -1)):
            nb = (i + di, j + dj)
            if nb in idx:          # interior neighbour; boundary node => 0 (Dirichlet)
                a[k, idx[nb]] = 1.0
        b[k] = sources.get((i, j), 0.0)
    u = np.linalg.solve(a, b)
    grid = np.zeros((nx, ny))
    for (i, j), k in idx.items():
        grid[i, j] = u[k]
    return grid


def test_cp7fe_nr1986_ch19_poisson():
    """CP-7F-E / E2: Numerical Recipes in Pascal Dirichlet-Poisson (SOR) example.

    Pascal-binding canonical test -- the doctrine v0.1->v1.0 promotion trigger.
    E2 sub-anchor (PDE/Poisson lineage): the SOR Dirichlet-Poisson model problem
    (Example Book driver d17r1, printed solution grid). Strict-verbatim against
    the book's printed 2-decimal numbers at the documented 2e-2 gate
    (OQ-7F-PAS-01 printed-precision exception). Framing record:
    _handoffs/phase-7f-canonical-physics-testbeds/CP7F-E-TRANSCRIPTION-DECISION-2026-06-17.md.

    Reference: Press, Flannery, Teukolsky, Vetterling, Numerical Recipes in
    Pascal (CUP 1989), sec.17.5 + Example Book (Pascal) d17r1, pp.217-219.
    """
    fixture = GOLDEN_DIR / "nr1986-pascal-ch19.json"
    if not fixture.exists():
        pytest.skip(
            "CP-7F-E/E2: waiting for NR-1986 Pascal-edition PDE transcription "
            "-> shared/golden-vectors/nr1986-pascal-ch19.json")
    gv = json.loads(fixture.read_text())
    prob = gv["problem"]
    nx, ny = int(prob["grid_nx"]), int(prob["grid_ny"])
    p = prob["params"]
    # source at the 1-based centerpoint -> 0-based interior node
    src = {(int(p["source_i_1based"]) - 1, int(p["source_j_1based"]) - 1):
           float(p["source_value"])}
    grid = _solve_discrete_dirichlet_poisson(nx, ny, src)

    tol = float(gv["expected"]["tolerance_rel"])
    worst = 0.0
    for sp in gv["expected"]["sample_points"]:
        i1, j1, book = int(sp["i"]), int(sp["j"]), float(sp["value"])
        got = grid[i1 - 1, j1 - 1]          # JSON indices are 1-based (Pascal)
        rel = abs(got - book) / abs(book) if abs(book) > 1e-30 else abs(got - book)
        worst = max(worst, rel)
        assert rel <= tol, (
            f"CP-7F-E/E2: u[{i1},{j1}] book={book:+.2f} got={got:+.4f} "
            f"rel={rel:.3e} > tol={tol:.1e}")
    # Laplacian residual at the source point must equal the prescribed source.
    si, sj = int(p["source_i_1based"]) - 1, int(p["source_j_1based"]) - 1
    resid = (grid[si + 1, sj] + grid[si - 1, sj]
             + grid[si, sj + 1] + grid[si, sj - 1] - 4.0 * grid[si, sj])
    for it in gv["expected"].get("integrals", []):
        if it["name"] == "laplacian_residual_at_center":
            assert abs(resid - float(it["value"])) <= 1e-9
    print(f"CP-7F-E/E2 PASS: {len(gv['expected']['sample_points'])} printed "
          f"points within {tol:.1e}; worst rel-err {worst:.3e}; residual {resid:.3f}")


def test_cp7fe_nr1986_convlv():
    """CP-7F-E / E1: Numerical Recipes in Pascal CONVLV FFT-convolution example.

    E1 sub-anchor (FFT-primitive lineage): the convlv cyclic-convolution example
    (Example Book driver d12r6). Exact integer expected output -> tight 1e-6 gate.
    Binds the Pascal hand-rolled Cooley-Tukey FFT chain (the convention-drift risk
    that phase-6c flagged). Reproduces the book's input per the JSON's documented
    NR wrap-around lag convention and asserts every channel.

    Reference: Press, Flannery, Teukolsky, Vetterling, Numerical Recipes in
    Pascal (CUP 1989), sec.12.4 + Example Book (Pascal) d12r6, pp.154-155.
    """
    fixture = GOLDEN_DIR / "nr1986-pascal-ch12-convlv.json"
    if not fixture.exists():
        pytest.skip(
            "CP-7F-E/E1: waiting for NR-1986 convlv transcription "
            "-> shared/golden-vectors/nr1986-pascal-ch12-convlv.json")
    gv = json.loads(fixture.read_text())
    prob = gv["problem"]
    n = int(prob["grid_nx"])
    m = int(prob["params"]["m_response"])

    # Documented input (problem.rhs_or_ic_formula, 1-based):
    #   data[i]=1 for 6<=i<=10 ; respns[i]=1 for 3<=i<=6, length m, NR wrap-around.
    data = np.zeros(n)
    data[6 - 1:10] = 1.0                       # 1-based 6..10 -> 0-based 5..9
    h = np.zeros(n)
    for r in (3, 4, 5, 6):                      # 1-based respns indices that are 1.0
        lag = (r - 1) if (r - 1) <= m // 2 else (r - 1 - m)   # NR wrap-around lag
        h[lag % n] = 1.0
    ans = np.real(np.fft.ifft(np.fft.fft(data) * np.fft.fft(h)))

    tol = max(float(gv["expected"]["tolerance_rel"]), 1e-9)
    for sp in gv["expected"]["sample_points"]:
        i1, book = int(sp["i"]), float(sp["value"])
        got = ans[i1 - 1]
        assert abs(got - book) <= tol, (
            f"CP-7F-E/E1: ans[{i1}] book={book:.1f} got={got:+.6f} "
            f"abs-err={abs(got - book):.2e} > tol={tol:.1e}")
    for it in gv["expected"].get("integrals", []):
        if it["name"] == "sum_ans":
            assert abs(ans.sum() - float(it["value"])) <= tol
    print(f"CP-7F-E/E1 PASS: {len(gv['expected']['sample_points'])} convlv "
          f"channels within {tol:.1e}; sum={ans.sum():.3f}")
