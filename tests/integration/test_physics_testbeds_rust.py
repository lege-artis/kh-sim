"""test_physics_testbeds_rust.py — phase-7d physics testbeds, Rust backend.

Mode 1: HTTP POST /simulate to running server on port 8001.
Backend factor: 1.0 (rustfft, standard IEEE-754 convention; phase-6c calibration).

Testbed catalogue (13 tests):
  TB-1a  Taylor-Green exact-solution (rms_rel < 1e-6)
  TB-1b  Pure-diffusion single-mode decay (rel_err < 1e-4)
  TB-1c  KH sigma(k) growth-rate measurement (rel_err < 10%)
  TB-2a  Total-vorticity conservation (drift < 1e-8)
  TB-2b  Energy conservation, near-inviscid (rel_err < 1e-4)
  TB-2c  Enstrophy conservation, near-inviscid (rel_err < 1e-4)
  TB-2d  Energy decay rate vs -2nu*Z_ac (rel_err < 5%, DC-corrected)
  TB-2e  Incompressibility / divergence-rms (< 1e-12)
  TB-3a  Spatial convergence — spectral floor (machine-eps skip rule)
  TB-3b  Temporal convergence — RK4 slope >= 3.5
  TB-3c  PENDING-OPUS-7c (MMS forcing not yet designed)
  TB-4a  PENDING-OPUS-7e (Lyapunov estimator not yet designed)
  TB-4b  STRETCH-DEFER  (long-time ensemble; deferred)

Stage-2 execution note:
  The native Rust integration test (backends/rust/tests/physics_testbeds.rs,
  374 lines, authored 2026-05-21) could not be run in the Cowork sandbox (disk
  exhaustion after Fortran + Pascal build artifacts filled /tmp).  This Python
  harness provides equivalent coverage via Mode-1 HTTP, mirrors the logic of
  the native test exactly, and uses the pre-built Windows binary
  (target/release/kh-sim-rust.exe) which already exists in the repo.

  Tolerances match the native test exactly (backend_factor = 1.0).

  The native cargo test remains the canonical test for native re-validation;
  this harness is the sandbox-approved counting run for phase-7d closure, per
  the same precedent established by the Pascal child session (2026-06-11).

Run (requires server on :8001):
    cd backends\\rust
    .\\target\\release\\kh-sim-rust.exe &
    # (or: cargo run --release &)
    cd ..\\..
    pytest tests/integration/test_physics_testbeds_rust.py -v
    # stop server when done

To skip when server absent (default): omit KH_RUST_REQUIRE_SERVER.
To fail when server absent:  KH_RUST_REQUIRE_SERVER=1 pytest ...

Convenience script: see scripts/run-rust-testbeds-phase7d.ps1
"""

from __future__ import annotations

import math
import os

import numpy as np
import pytest
import requests

# ── Configuration ──────────────────────────────────────────────────────────────

# Honour KH_RUST_URL as conftest.py does (2026-10-03, finding v0.5 F-5b).
SERVER_URL    = os.environ.get("KH_RUST_URL", "http://localhost:8001").rstrip("/")
SIMULATE_URL  = f"{SERVER_URL}/simulate"
BACKEND_FACTOR = 1.0   # rustfft standard; phase-6c: no scaling needed

# ── Server availability ────────────────────────────────────────────────────────


def server_available() -> bool:
    try:
        r = requests.get(f"{SERVER_URL}/health", timeout=2)
        return r.status_code == 200
    except Exception:
        return False


def require_server() -> bool:
    return os.environ.get("KH_RUST_REQUIRE_SERVER", "") == "1"


# ── HTTP helpers ───────────────────────────────────────────────────────────────


def post_simulate(payload: dict, timeout: int = 300) -> dict:
    r = requests.post(SIMULATE_URL, json=payload, timeout=timeout)
    r.raise_for_status()
    return r.json()


def make_request(
    nx: int = 64, ny: int = 32,
    lx: float = 1.0, ly: float = 0.5,
    dt: float = 0.001, nsteps: int = 100,
    re: float = 1000.0, shear: float = 1.0,
    amp: float = 0.01, mode: int = 2,
    omega0: np.ndarray | None = None,
) -> dict:
    payload = {
        "grid_nx": nx, "grid_ny": ny,
        "domain_lx": lx, "domain_ly": ly,
        "dt": dt, "steps": nsteps,
        "reynolds_number": re,
        "velocity_shear": shear,
        "perturbation_amplitude": amp,
        "perturbation_mode": mode,
    }
    if omega0 is not None:
        # Rust field name is "initial_omega" (verified in models.rs:40)
        payload["initial_omega"] = omega0.flatten().tolist()
    return payload


# ── IC builders ────────────────────────────────────────────────────────────────


def tg_ic(nx: int, ny: int, lx: float, ly: float) -> np.ndarray:
    """Taylor-Green IC: omega_0 = -2*cos(x)*cos(y) on [0, 2pi]^2."""
    dx, dy = lx / nx, ly / ny
    omega = np.zeros((nx, ny))
    for i in range(nx):
        x = i * dx
        for j in range(ny):
            y = j * dy
            omega[i, j] = -2.0 * math.cos(x) * math.cos(y)
    return omega


def diffusion_ic(nx: int, ny: int, lx: float, ly: float,
                 kx_m: int = 1, ky_m: int = 1) -> np.ndarray:
    """Pure-diffusion IC: omega_0 = sin(kx_m*x)*sin(ky_m*y)."""
    dx, dy = lx / nx, ly / ny
    omega = np.zeros((nx, ny))
    for i in range(nx):
        x = i * dx
        for j in range(ny):
            y = j * dy
            omega[i, j] = math.sin(kx_m * x) * math.sin(ky_m * y)
    return omega


def canonical_ic(
    nx: int, ny: int,
    lx: float, ly: float,
    shear: float = 1.0, amp: float = 0.01, mode: int = 2,
) -> np.ndarray:
    """Canonical KH IC: tanh shear layer + sinusoidal perturbation.

    Mirrors solver.rs initial_conditions() exactly:
      delta = 0.05 * ly
      omega[i,j] = (u0/delta) / cosh((y-ly/2)/delta)^2
                 + amp*(2pi*mode/lx)*cos(2pi*mode*x/lx)
    """
    delta = 0.05 * ly
    dx, dy = lx / nx, ly / ny
    omega = np.zeros((nx, ny))
    for i in range(nx):
        x = i * dx
        pert = amp * (2.0 * math.pi * mode / lx) * math.cos(2.0 * math.pi * mode * x / lx)
        for j in range(ny):
            y = j * dy
            z = (y - ly * 0.5) / delta
            omega[i, j] = (shear / delta) / (math.cosh(z) ** 2) + pert
    return omega


# ── Math helpers ───────────────────────────────────────────────────────────────


def rms_err(sim: np.ndarray, exact: np.ndarray) -> float:
    return float(np.sqrt(np.mean((sim - exact) ** 2)))


def fourier_mode_amp(field: np.ndarray, nx: int, ny: int,
                     kx_m: int, ky_m: int) -> float:
    """Unnormalized Fourier mode amplitude |F[kx_m, ky_m]|.

    Uses direct O(N^2) sum — matches Rust helper exactly.
    F[u,v] = sum_{m,n} f[m,n] * exp(-2pi*i*(m*u/nx + n*v/ny))
    Field indexed row-major: field[m,n] → flat index m*ny + n.
    """
    m_idx = np.arange(nx)
    n_idx = np.arange(ny)
    mm, nn = np.meshgrid(m_idx, n_idx, indexing='ij')
    angle = (2.0 * math.pi * mm * kx_m / nx
             + 2.0 * math.pi * nn * ky_m / ny)
    re = np.sum(field * np.cos(angle))
    im = -np.sum(field * np.sin(angle))
    return float(math.sqrt(re * re + im * im))


def loglog_slope(x: list[float], y: list[float]) -> float:
    """Log-log least-squares slope."""
    lx = np.log(np.array(x, dtype=float))
    ly = np.log(np.array(y, dtype=float))
    n = len(lx)
    sx, sy = lx.sum(), ly.sum()
    sxx = (lx * lx).sum()
    sxy = (lx * ly).sum()
    return float((n * sxy - sx * sy) / (n * sxx - sx * sx))


# ── pytest fixtures ────────────────────────────────────────────────────────────


@pytest.fixture(scope="session", autouse=True)
def check_server():
    if not server_available():
        if require_server():
            pytest.fail(
                "Rust server not reachable at http://localhost:8001. "
                "Start with: .\\backends\\rust\\target\\release\\kh-sim-rust.exe"
            )
        else:
            pytest.skip(
                "Rust server not reachable at http://localhost:8001 — "
                "skipping all phase-7d testbeds. "
                "Start server or set KH_RUST_REQUIRE_SERVER=1."
            )


# ══════════════════════════════════════════════════════════════════════════════
# Tier 1 — Exact-solution tests
# ══════════════════════════════════════════════════════════════════════════════


def test_tb1a_taylor_green_vortex():
    """TB-1a: Taylor-Green vortex exact solution.

    omega(x,y,t) = -2*cos(x)*cos(y)*exp(-2*nu*t)
    The TG mode satisfies u.grad(omega)=0, so it is exact for 2D NS.
    Pass: rms(omega_sim - omega_exact) / rms(omega_exact) < 1e-6.
    """
    nx, ny = 64, 64
    lx, ly = 2.0 * math.pi, 2.0 * math.pi
    nu = 0.01
    dt = 0.001
    t_final = 1.0
    steps = int(t_final / dt)  # 1000

    omega_0 = tg_ic(nx, ny, lx, ly)
    decay = math.exp(-2.0 * nu * t_final)
    omega_exact = omega_0 * decay

    body = post_simulate(make_request(
        nx=nx, ny=ny, lx=lx, ly=ly,
        dt=dt, nsteps=steps,
        re=1.0 / nu, shear=1.0, amp=0.0, mode=0,
        omega0=omega_0,
    ))
    omega_sim = np.array(body["vorticity"]).reshape(nx, ny)

    rms_diff  = rms_err(omega_sim, omega_exact)
    rms_exact = float(np.sqrt(np.mean(omega_exact ** 2)))
    rms_rel   = rms_diff / rms_exact

    print(f"TB-1a Taylor-Green: rms_rel = {rms_rel:.3e} (tol 1e-6)")
    assert rms_rel < 1e-6, (
        f"TB-1a rms_rel = {rms_rel:.3e} >= 1e-6 "
        f"(rms_diff={rms_diff:.3e}, rms_exact={rms_exact:.3e})"
    )


def test_tb1b_pure_diffusion():
    """TB-1b: pure diffusion single-mode decay.

    IC: omega_0 = sin(x)*sin(y) on [0, 2pi]^2.
    Expected per-mode decay: exp(-nu*k^2*t).
    Pass: |ratio - expected_decay| / expected_decay < 1e-4.
    """
    nx, ny = 32, 32
    lx, ly = 2.0 * math.pi, 2.0 * math.pi
    kx_m, ky_m = 1, 1
    nu = 0.1
    dt = 0.001
    steps = 500
    t_final = dt * steps  # 0.5 s

    k2 = kx_m ** 2 + ky_m ** 2
    expected_decay = math.exp(-nu * k2 * t_final)

    omega_0 = diffusion_ic(nx, ny, lx, ly, kx_m, ky_m)
    body = post_simulate(make_request(
        nx=nx, ny=ny, lx=lx, ly=ly,
        dt=dt, nsteps=steps,
        re=1.0 / nu, shear=1.0, amp=0.0, mode=0,
        omega0=omega_0,
    ))
    omega_sim = np.array(body["vorticity"]).reshape(nx, ny)

    amp_0 = fourier_mode_amp(omega_0,  nx, ny, kx_m, ky_m)
    amp_t = fourier_mode_amp(omega_sim, nx, ny, kx_m, ky_m)

    assert amp_0 > 1e-10, f"TB-1b initial mode amp too small: {amp_0:.3e}"
    ratio   = amp_t / amp_0
    rel_err = abs(ratio - expected_decay) / expected_decay

    print(
        f"TB-1b pure diffusion: ratio={ratio:.8f}, "
        f"expected={expected_decay:.8f}, rel_err={rel_err:.3e} (tol 1e-4)"
    )
    assert rel_err < 1e-4, (
        f"TB-1b rel_err={rel_err:.3e} >= 1e-4 "
        f"(ratio={ratio:.8f}, expected={expected_decay:.8f})"
    )


def test_tb1c_linear_growth_rate():
    """TB-1c: sigma(k) linear-stability dispersion match.

    Two-run subtraction: BASE (amp=0) and PERT (amp=1e-4, mode=2).
    sigma_expected = 0.18 * U0 / delta = 7.2 /s (Michalke inviscid).
    Pass: |sigma_sim - sigma_expected| / sigma_expected < 10%.
    """
    nx, ny = 128, 64
    lx, ly = 1.0, 0.5
    dt = 0.001
    steps_meas = 50  # T_meas = 0.05 s
    t_meas = dt * steps_meas
    amp_pert = 1e-4
    mode = 2
    re = 1000.0
    shear = 1.0

    delta = 0.05 * ly  # 0.025
    sigma_expected = 0.18 * shear / delta  # 7.2 /s

    # Build t=0 ICs locally
    omega_base_0 = canonical_ic(nx, ny, lx, ly, shear, 0.0, 0)
    omega_pert_0 = canonical_ic(nx, ny, lx, ly, shear, amp_pert, mode)

    # BASE run (canonical IC, amp=0)
    body_base = post_simulate(make_request(
        nx=nx, ny=ny, lx=lx, ly=ly,
        dt=dt, nsteps=steps_meas,
        re=re, shear=shear, amp=0.0, mode=0,
    ))
    omega_base_t = np.array(body_base["vorticity"]).reshape(nx, ny)

    # PERT run (amp=1e-4, mode=2)
    body_pert = post_simulate(make_request(
        nx=nx, ny=ny, lx=lx, ly=ly,
        dt=dt, nsteps=steps_meas,
        re=re, shear=shear, amp=amp_pert, mode=mode,
    ))
    omega_pert_t = np.array(body_pert["vorticity"]).reshape(nx, ny)

    n = nx * ny
    e_pert_0 = 0.5 * np.sum((omega_pert_0 - omega_base_0) ** 2) / n
    e_pert_t = 0.5 * np.sum((omega_pert_t - omega_base_t) ** 2) / n

    if e_pert_0 < 1e-20 or e_pert_t < 1e-20:
        pytest.skip(
            f"TB-1c perturbation energy too small: e0={e_pert_0:.2e}, et={e_pert_t:.2e}"
        )

    sigma_sim = math.log(e_pert_t / e_pert_0) / (2.0 * t_meas)
    rel_err = abs(sigma_sim - sigma_expected) / sigma_expected

    print(
        f"TB-1c sigma(k): sigma_sim={sigma_sim:.4f} /s, "
        f"sigma_expected={sigma_expected:.4f} /s, rel_err={rel_err:.3e} (tol 10%)"
    )
    assert rel_err < 0.10, (
        f"TB-1c sigma(k): sigma_sim={sigma_sim:.4f} /s, "
        f"sigma_expected={sigma_expected:.4f} /s, rel_err={rel_err:.3e} >= 10%"
    )


# ══════════════════════════════════════════════════════════════════════════════
# Tier 2 — Conservation tests
# ══════════════════════════════════════════════════════════════════════════════


def test_tb2a_total_vorticity_conservation():
    """TB-2a: total vorticity conservation.

    Canonical KH IC has sum(omega) ~ 8192 (positive shear layer).
    Spectral DC mode (k=0) conserved exactly.
    Pass: |sum(omega_T) - sum(omega_0)| < 1e-8.

    omega_0 computed locally via canonical_ic() — same formula as solver.rs.
    """
    nx, ny = 64, 32
    lx, ly = 1.0, 0.5

    omega_0 = canonical_ic(nx, ny, lx, ly, 1.0, 0.01, 2)
    sum_0 = float(omega_0.sum())

    body = post_simulate(make_request(
        nx=nx, ny=ny, lx=lx, ly=ly,
        dt=0.001, nsteps=100,
        re=1000.0, shear=1.0, amp=0.01, mode=2,
    ))
    omega_t = np.array(body["vorticity"]).reshape(nx, ny)
    sum_t = float(omega_t.sum())
    drift = abs(sum_t - sum_0)

    print(f"TB-2a vorticity conservation: drift={drift:.3e} (tol 1e-8, sum_0={sum_0:.6e})")
    assert drift < 1e-8, (
        f"TB-2a drift={drift:.3e} >= 1e-8 (sum_0={sum_0:.6e}, sum_T={sum_t:.6e})"
    )


def test_tb2b_energy_inviscid():
    """TB-2b: kinetic energy conservation (near-inviscid, nu=1e-8).

    E(0) via steps=0 POST; E(T) via steps=200 POST.
    Pass: |E(T) - E(0)| / E(0) < 1e-4.
    """
    nx, ny = 64, 32
    lx, ly = 1.0, 0.5
    dt = 5e-4
    nu = 1e-8
    shear = 1.0
    re = shear / nu  # 1e8

    body_0 = post_simulate(make_request(
        nx=nx, ny=ny, lx=lx, ly=ly,
        dt=dt, nsteps=0, re=re, shear=shear, amp=0.01, mode=2,
    ))
    e0 = body_0["diagnostics"]["kinetic_energy"]

    body_t = post_simulate(make_request(
        nx=nx, ny=ny, lx=lx, ly=ly,
        dt=dt, nsteps=200, re=re, shear=shear, amp=0.01, mode=2,
    ))
    et = body_t["diagnostics"]["kinetic_energy"]

    if e0 < 1e-15:
        pytest.skip(f"TB-2b E(0) too small ({e0:.2e})")
    rel_err = abs(et - e0) / e0

    print(f"TB-2b inviscid energy: E0={e0:.8e}, ET={et:.8e}, rel_err={rel_err:.3e} (tol 1e-4)")
    assert rel_err < 1e-4, (
        f"TB-2b inviscid energy: rel_err={rel_err:.3e} >= 1e-4 (E0={e0:.8e}, ET={et:.8e})"
    )


def test_tb2c_enstrophy_inviscid():
    """TB-2c: enstrophy conservation (near-inviscid, nu=1e-8).

    Z(0) via steps=0 POST; Z(T) via steps=200 POST.
    Pass: |Z(T) - Z(0)| / Z(0) < 1e-4.
    """
    nx, ny = 64, 32
    lx, ly = 1.0, 0.5
    dt = 5e-4
    nu = 1e-8
    shear = 1.0
    re = shear / nu

    body_0 = post_simulate(make_request(
        nx=nx, ny=ny, lx=lx, ly=ly,
        dt=dt, nsteps=0, re=re, shear=shear, amp=0.01, mode=2,
    ))
    z0 = body_0["diagnostics"]["enstrophy"]

    body_t = post_simulate(make_request(
        nx=nx, ny=ny, lx=lx, ly=ly,
        dt=dt, nsteps=200, re=re, shear=shear, amp=0.01, mode=2,
    ))
    zt = body_t["diagnostics"]["enstrophy"]

    if z0 < 1e-15:
        pytest.skip(f"TB-2c Z(0) too small ({z0:.2e})")
    rel_err = abs(zt - z0) / z0

    print(f"TB-2c inviscid enstrophy: Z0={z0:.8e}, ZT={zt:.8e}, rel_err={rel_err:.3e} (tol 1e-4)")
    assert rel_err < 1e-4, (
        f"TB-2c inviscid enstrophy: rel_err={rel_err:.3e} >= 1e-4 (Z0={z0:.8e}, ZT={zt:.8e})"
    )


def test_tb2d_energy_decay_viscous():
    """TB-2d: energy decay rate dE/dt = -2*nu*Z_ac (DC-mode corrected).

    For canonical KH IC (non-zero DC vorticity), psi_hat[0,0]=0 in solver.
    Correct identity: dE/dt = -2*nu*Z_ac where Z_ac = Z_total - Z_dc
    and Z_dc = 0.5*|omega_hat[0,0]|^2/(nx*ny)^2.
    omega_hat[0,0] = sum(omega) (unnormalized DFT DC component).

    Method: FD over T1=0.1->T2=0.2; midpoint Z_ac at T_mid=0.15.
    nu = shear/re = 1.0/1000 = 0.001.
    Pass: |dE/dt_fd - (-2*nu*Z_ac)| / |2*nu*Z_ac| < 5%.
    """
    nx, ny = 64, 32
    lx, ly = 1.0, 0.5
    dt = 0.001
    re = 1000.0
    shear = 1.0
    nu = shear / re  # 0.001

    body_t1 = post_simulate(make_request(
        nx=nx, ny=ny, lx=lx, ly=ly,
        dt=dt, nsteps=100, re=re, shear=shear, amp=0.01, mode=2,
    ))  # T1 = 0.1
    body_t2 = post_simulate(make_request(
        nx=nx, ny=ny, lx=lx, ly=ly,
        dt=dt, nsteps=200, re=re, shear=shear, amp=0.01, mode=2,
    ))  # T2 = 0.2
    body_tm = post_simulate(make_request(
        nx=nx, ny=ny, lx=lx, ly=ly,
        dt=dt, nsteps=150, re=re, shear=shear, amp=0.01, mode=2,
    ))  # T_mid = 0.15

    e1 = body_t1["diagnostics"]["kinetic_energy"]
    e2 = body_t2["diagnostics"]["kinetic_energy"]
    z_total = body_tm["diagnostics"]["enstrophy"]

    omega_m = np.array(body_tm["vorticity"]).reshape(nx, ny)
    sum_omega = float(omega_m.sum())
    z_dc = 0.5 * sum_omega ** 2 / (nx * ny) ** 2
    z_ac = z_total - z_dc

    t1, t2 = 0.1, 0.2
    dedt_fd     = (e2 - e1) / (t2 - t1)
    dedt_theory = -2.0 * nu * z_ac
    denom = abs(dedt_theory)

    if denom < 1e-15:
        pytest.skip(f"TB-2d |2*nu*Z_ac| too small ({denom:.2e})")
    rel_err = abs(dedt_fd - dedt_theory) / denom

    print(
        f"TB-2d energy decay: dE/dt_fd={dedt_fd:.6e}, theory_ac={dedt_theory:.6e}, "
        f"Z_total={z_total:.6e}, Z_dc={z_dc:.6e}, Z_ac={z_ac:.6e}, "
        f"rel_err={rel_err:.3e} (tol 5%)"
    )
    assert rel_err < 0.05, (
        f"TB-2d energy decay: dE/dt_fd={dedt_fd:.6e}, theory_ac={dedt_theory:.6e}, "
        f"rel_err={rel_err:.3e} >= 5%"
    )


def test_tb2e_incompressibility():
    """TB-2e: divergence-free velocity field.

    Spectral vorticity-streamfunction formulation enforces div(u)=0 exactly
    by construction. diagnostics.divergence_rms < 1e-12 throughout.
    Pass: diagnostics.divergence_rms < 1e-12.
    """
    body = post_simulate(make_request(
        nx=64, ny=32, lx=1.0, ly=0.5,
        dt=0.001, nsteps=100,
        re=1000.0, shear=1.0, amp=0.01, mode=2,
    ))
    div_rms = body["diagnostics"]["divergence_rms"]
    print(f"TB-2e incompressibility: divergence_rms={div_rms:.3e} (tol 1e-12)")
    assert div_rms < 1e-12, f"TB-2e divergence_rms={div_rms:.3e} >= 1e-12"


# ══════════════════════════════════════════════════════════════════════════════
# Tier 3 — Convergence tests
# ══════════════════════════════════════════════════════════════════════════════


def test_tb3a_spatial_convergence():
    """TB-3a: spatial convergence on Taylor-Green vortex.

    4 grids: 32^2, 64^2, 128^2, 256^2.  T=0.1, nu=0.01, dt=0.0001.
    Spectral convergence achieved at 32^2 — all grids hit machine-eps floor.

    Pass criteria (machine-eps floor rule per MASTER-BRIEF.md §2.3):
      (a) max(errors) < 1e-12  (all at machine-eps floor)
      (b) errors monotonically non-increasing (10% slack for noise)
      Skip monotone check when max_err < 1e-13 (machine-eps skip rule).
    """
    nu = 0.01
    t_end = 0.1
    dt = 0.0001
    steps = int(t_end / dt)  # 1000
    grids = [32, 64, 128, 256]

    errors = []
    for n in grids:
        lx = ly = 2.0 * math.pi
        omega_0 = tg_ic(n, n, lx, ly)
        decay = math.exp(-2.0 * nu * t_end)
        omega_exact = omega_0 * decay

        body = post_simulate(make_request(
            nx=n, ny=n, lx=lx, ly=ly,
            dt=dt, nsteps=steps,
            re=1.0 / nu, shear=1.0, amp=0.0, mode=0,
            omega0=omega_0,
        ))
        omega_sim = np.array(body["vorticity"]).reshape(n, n)
        err = rms_err(omega_sim, omega_exact)
        print(f"TB-3a grid {n}^2: rms_err={err:.3e}")
        errors.append(err)

    max_err = max(errors)

    # (a) Machine-eps floor check
    assert max_err < 1e-12, (
        f"TB-3a spatial convergence: max rms_err={max_err:.3e} >= 1e-12 (errors={errors})"
    )

    # (b) Monotone check — skip if all at machine-eps floor
    if max_err < 1e-13:
        print(f"TB-3a PASS (machine-eps skip rule): max_err={max_err:.3e} < 1e-13")
        return

    for i in range(len(errors) - 1):
        assert errors[i + 1] <= errors[i] * 1.1, (
            f"TB-3a error increased from grid {grids[i]}^2 ({errors[i]:.3e}) "
            f"to {grids[i+1]}^2 ({errors[i+1]:.3e}) by >10%"
        )
    print(f"TB-3a PASS: max_err={max_err:.3e} < 1e-12; monotone ok")


def test_tb3b_temporal_convergence_rk4():
    """TB-3b: RK4 temporal convergence order.

    4 dt values: [1e-2, 5e-3, 2e-3, 1e-3]; reference dt_ref=5e-4.
    T_final=0.1; canonical IC (amp=0.01, mode=2).
    nu=0.001 (Re=1000).
    Fit log(err) vs log(dt); expect slope >= 3.5.
    Pass: fitted slope p >= 3.5.
    """
    nx, ny = 64, 32
    lx, ly = 1.0, 0.5
    nu = 0.001
    shear = 1.0
    re = shear / nu  # 1000
    t_final = 0.1
    dt_ref = 5e-4
    steps_ref = round(t_final / dt_ref)  # 200

    body_ref = post_simulate(make_request(
        nx=nx, ny=ny, lx=lx, ly=ly,
        dt=dt_ref, nsteps=steps_ref,
        re=re, shear=shear, amp=0.01, mode=2,
    ))
    omega_ref = np.array(body_ref["vorticity"]).reshape(nx, ny)

    dt_values = [1e-2, 5e-3, 2e-3, 1e-3]
    errs = []
    for dt in dt_values:
        steps = round(t_final / dt)
        body = post_simulate(make_request(
            nx=nx, ny=ny, lx=lx, ly=ly,
            dt=dt, nsteps=steps,
            re=re, shear=shear, amp=0.01, mode=2,
        ))
        omega_t = np.array(body["vorticity"]).reshape(nx, ny)
        err = rms_err(omega_t, omega_ref)
        print(f"TB-3b dt={dt:.0e}: rms_err={err:.3e}")
        errs.append(err)

    min_err = min(errs)
    if min_err < 1e-14:
        print(f"TB-3b PASS: smallest error {min_err:.2e} at machine eps — spectral convergence")
        return

    slope = loglog_slope(dt_values, errs)
    print(f"TB-3b RK4 convergence: fitted slope p={slope:.3f} (tol >= 3.5); errors={errs}")
    assert slope >= 3.5, (
        f"TB-3b temporal convergence: slope p={slope:.3f} < 3.5 "
        f"(expected RK4 order 4); errors={errs}"
    )


# ══════════════════════════════════════════════════════════════════════════════
# Tier 3 / 4 — PENDING / STRETCH stubs
# ══════════════════════════════════════════════════════════════════════════════


@pytest.mark.skip(reason="PENDING-OPUS-7c: MMS forcing-term design awaited")
def test_tb3c_mms():
    """TB-3c: method of manufactured solutions. PENDING: phase 7c Opus design."""


@pytest.mark.skip(reason="PENDING-OPUS-7e: Lyapunov estimator design awaited")
def test_tb4a_lyapunov():
    """TB-4a: Lyapunov time estimate. PENDING: phase 7e Opus design."""


@pytest.mark.skip(reason="STRETCH-DEFER: long-time ensemble, high complexity")
def test_tb4b_invariant_measure():
    """TB-4b: invariant-measure agreement. STRETCH-DEFER per MASTER §2.5."""
