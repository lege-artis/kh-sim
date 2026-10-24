"""
kh_physics.py — Kelvin-Helmholtz instability reference implementation
======================================================================
2D incompressible Navier-Stokes, vorticity-streamfunction formulation.
Doubly periodic domain, pseudo-spectral Poisson solver, RK4 time stepping.

This is the CANONICAL reference. All backend implementations (Rust, Scala,
C++, Fortran, Pascal) must validate their output against this module.

See: shared/physics/KH-PHYSICS.md for full physics documentation.
See: shared/api/openapi.yaml for REST API contract.

Usage:
    python kh_physics.py                        # run reference smoke test
    python kh_physics.py --nx 128 --ny 64 --steps 200 --re 1000
"""

import warnings
import numpy as np
import json
import math
import time
import argparse
from dataclasses import dataclass, field, asdict
from typing import Optional, List


# ── Data structures (mirrors REST API schema) ──────────────────────────────────

@dataclass
class SimulationRequest:
    grid_nx:               int   = 128
    grid_ny:               int   = 64
    domain_lx:             float = 1.0
    domain_ly:             float = 0.5
    dt:                    float = 0.001
    steps:                 int   = 100
    reynolds_number:       float = 1000.0
    velocity_shear:        float = 1.0
    perturbation_amplitude:float = 0.01
    perturbation_mode:     int   = 2
    initial_omega:         Optional[List[float]] = None  # flattened Nx*Ny
    # Phase-7f steady-forcing extension (OQ-7F-PY-01; see
    # development record, forcing-API design v0.1 for the locked convention).
    forcing_type:          str   = "none"   # "none" | "mode" | "band"
    forcing_amplitude:     float = 0.0
    forcing_kx:            int   = 1        # "mode": integer mode index along x
    forcing_ky:            int   = 1        # "mode": integer mode index along y
    forcing_kmin:          float = 0.0      # "band": inner integer-mode radius
    forcing_kmax:          float = 0.0      # "band": outer integer-mode radius
    # v0.2 (E2): rotating phases + Ekman drag (design doc s7; defaults = v0.1)
    forcing_omega0:        float = 0.0      # "band": phase-rotation rate scale
    linear_drag:           float = 0.0      # Ekman-type linear drag alpha
    # Time of the initial field (2026-10-10, CP-7F-B/D harness): a run continued
    # from a returned field passes the previous t_final, so time-dependent forcing
    # keeps its phase. t_final = t_start + steps*dt.
    t_start:               float = 0.0


@dataclass
class Diagnostics:
    kinetic_energy:  float = 0.0
    enstrophy:       float = 0.0
    max_vorticity:   float = 0.0
    divergence_rms:  float = 0.0
    cfl_peak:        float = 0.0   # max over all steps of dt*(max|u|/dx + max|v|/dy)


@dataclass
class SimulationResult:
    backend:          str   = "python-numpy"
    language:         str   = "Python"
    steps_completed:  int   = 0
    t_final:          float = 0.0
    grid_nx:          int   = 0
    grid_ny:          int   = 0
    u_velocity:       List[float] = field(default_factory=list)
    v_velocity:       List[float] = field(default_factory=list)
    vorticity:        List[float] = field(default_factory=list)
    pressure:         List[float] = field(default_factory=list)  # streamfunction ψ
    diagnostics:      Diagnostics = field(default_factory=Diagnostics)
    compute_time_ms:  float = 0.0
    warnings:         List[str] = field(default_factory=list)   # K5 (modelling rules v0.1)


def request_warnings(dt: float, dx: float, dy: float, nu: float, cfl_peak: float) -> List[str]:
    """K5 warnings (modelling rules v0.1 s2): computed, but possibly not resolved.
    The same three conditions, same order and wording prefix, in every backend."""
    out: List[str] = []
    if cfl_peak > 1.0:
        out.append(f"cfl_advective: peak Courant number {cfl_peak:.3g} exceeds 1")
    diff = 2.0 * nu * dt / min(dx, dy) ** 2
    if diff > 1.0:
        out.append(f"cfl_diffusive: 2*nu*dt/min(dx,dy)^2 = {diff:.3g} exceeds 1")
    if nu == 0.0:
        out.append("inviscid: nu = 0, no physical dissipation; enstrophy accumulates at the 2/3 cutoff")
    return out


# ── Physics kernels ────────────────────────────────────────────────────────────

def _wavenumbers(nx: int, ny: int, dx: float, dy: float):
    """
    Return 2D wavenumber arrays for rfft2 output shape (nx, ny//2+1).
    rfft2 applies full FFT along axis-0 (kx) and real FFT along axis-1 (ky).
    """
    kx = 2.0 * np.pi * np.fft.fftfreq(nx, d=dx)    # full FFT axis-0: shape (nx,)
    ky = 2.0 * np.pi * np.fft.rfftfreq(ny, d=dy)   # real FFT axis-1: shape (ny//2+1,)
    KX, KY = np.meshgrid(kx, ky, indexing='ij')    # shape (nx, ny//2+1) ✓
    return KX, KY


def initial_conditions(nx: int, ny: int, lx: float, ly: float,
                       U0: float, delta: float, amp: float, mode: int
                       ) -> np.ndarray:
    """
    Compute initial vorticity field analytically.

    IC velocity (KH-GE-4; the sign of u corrected 2026-10-12, the field unchanged):
        u(x,y) = -U0 * tanh((y - Ly/2) / delta)    upper stream in -x
        v(x,y) = amp * sin(2*pi*mode*x / Lx)

    Vorticity ω = dv/dx - du/dy:
        = amp*(2*pi*mode/Lx)*cos(2*pi*mode*x/Lx)  [from v]
        + U0/delta * sech²((y-Ly/2)/delta)         [from -du/dy]
    """
    dx = lx / nx
    dy = ly / ny
    x = np.arange(nx) * dx
    y = np.arange(ny) * dy
    X, Y = np.meshgrid(x, y, indexing='ij')

    # Vorticity from shear layer (positive for U0 > 0: upper stream in -x)
    omega_shear = U0 / delta / np.cosh((Y - ly / 2.0) / delta) ** 2

    # Vorticity perturbation from sinusoidal v-perturbation
    omega_pert = amp * (2.0 * np.pi * mode / lx) * np.cos(2.0 * np.pi * mode * X / lx)

    return omega_shear + omega_pert


def solve_poisson(omega_hat: np.ndarray, KX: np.ndarray, KY: np.ndarray,
                  nx: int, ny: int) -> np.ndarray:
    """
    Spectral Poisson solve: ∇²ψ = -ω  →  ψ̂ = ω̂ / (kx²+ky²)

    Returns psi_hat (rfft2 coefficients).
    Mean mode (kx=ky=0) set to zero: fixes ψ_mean = 0.
    """
    K2 = KX ** 2 + KY ** 2
    K2[0, 0] = 1.0          # avoid divide-by-zero; corrected below
    psi_hat = omega_hat / K2
    psi_hat[0, 0] = 0.0     # zero mean streamfunction
    return psi_hat


def velocity_from_psi_hat(psi_hat: np.ndarray, KX: np.ndarray, KY: np.ndarray,
                           nx: int, ny: int) -> tuple:
    """
    Recover velocity from streamfunction via spectral differentiation.
        u =  dψ/dy  →  û = i·ky·ψ̂
        v = -dψ/dx  →  v̂ = -i·kx·ψ̂
    """
    u = np.fft.irfft2(1j * KY * psi_hat, s=(nx, ny))
    v = np.fft.irfft2(-1j * KX * psi_hat, s=(nx, ny))
    return u, v


def dealias_two_thirds(field_hat: np.ndarray, KX: np.ndarray, KY: np.ndarray) -> np.ndarray:
    """KH-NM-4: 2/3-rule dealiasing (Orszag 1971).
    Zeros rfft2 modes where |kx| > 2/3·kx_max or ky > 2/3·ky_max."""
    kx_cut = (2.0 / 3.0) * np.max(np.abs(KX))
    ky_cut = (2.0 / 3.0) * np.max(KY)
    out = field_hat.copy()
    out[(np.abs(KX) > kx_cut) | (KY > ky_cut)] = 0.0
    return out


def cfl_dt_bound(u_max: float, v_max: float, dx: float, dy: float, nu: float) -> float:
    """KH-NM-5: CFL-limited time step. Returns min of advective and diffusive bounds."""
    dt_adv  = min(dx / (u_max + 1e-30), dy / (v_max + 1e-30))
    # nu = 0 (inviscid) is a valid limit of the model: no diffusive bound
    # (modelling rules v0.1 s3, D; ruling 2026-10-07).
    dt_diff = math.inf if nu == 0.0 else min(dx ** 2 / (2.0 * nu), dy ** 2 / (2.0 * nu))
    return min(dt_adv, dt_diff)


def vorticity_rhs(omega: np.ndarray, u: np.ndarray, v: np.ndarray,
                  nu: float, KX: np.ndarray, KY: np.ndarray,
                  nx: int, ny: int) -> np.ndarray:
    """RHS of vorticity transport equation with 2/3-rule dealiasing on nonlinear term."""
    omega_hat = np.fft.rfft2(omega)
    K2 = KX ** 2 + KY ** 2

    domega_dx = np.fft.irfft2(1j * KX * omega_hat, s=(nx, ny))
    domega_dy = np.fft.irfft2(1j * KY * omega_hat, s=(nx, ny))
    laplacian_omega = np.fft.irfft2(-K2 * omega_hat, s=(nx, ny))

    # Nonlinear term with 2/3-rule dealiasing (KH-NM-4, Orszag 1971)
    nl_hat = np.fft.rfft2(-u * domega_dx - v * domega_dy)
    nl_hat = dealias_two_thirds(nl_hat, KX, KY)
    nonlinear = np.fft.irfft2(nl_hat, s=(nx, ny))

    return nonlinear + nu * laplacian_omega


def build_forcing_field(nx: int, ny: int, lx: float, ly: float,
                        ftype: str, amp: float,
                        fkx: int, fky: int,
                        fkmin: float, fkmax: float) -> Optional[np.ndarray]:
    """
    Steady (time-independent) vorticity forcing field F(x,y) -- phase-7f extension.
    Locked cross-backend convention: development record, forcing-API design v0.1 s3.

      "none": returns None (also when amp == 0).
      "mode": F = A * cos(2*pi*fkx*x/lx) * cos(2*pi*fky*y/ly)   (CP-7F-D)
      "band": F = A/sqrt(M) * sum over half-plane integer modes (ix,iy),
              fkmin <= sqrt(ix^2+iy^2) <= fkmax, of
              cos(2*pi*(ix*x/lx + iy*y/ly) + phi(ix,iy)),
              phi(ix,iy) = 2*pi*frac(ix*sqrt(2) + iy*sqrt(3))    (CP-7F-B)
              Requires lx == ly.
    """
    if ftype == "none" or amp == 0.0:
        return None

    x = np.arange(nx) * (lx / nx)
    y = np.arange(ny) * (ly / ny)
    X, Y = np.meshgrid(x, y, indexing="ij")

    if ftype == "mode":
        return amp * np.cos(2.0 * np.pi * fkx * X / lx) \
                   * np.cos(2.0 * np.pi * fky * Y / ly)

    if ftype == "band":
        if lx != ly:
            raise ValueError("forcing_type='band' requires a square domain (lx == ly)")
        kmax_i = int(np.floor(fkmax))
        modes = []
        for ix in range(0, kmax_i + 1):
            for iy in range(-kmax_i, kmax_i + 1):
                if ix == 0 and iy <= 0:
                    continue  # half-plane: avoid double-counting conjugate modes
                r = np.sqrt(float(ix * ix + iy * iy))
                if fkmin <= r <= fkmax:
                    modes.append((ix, iy))
        if not modes:
            return None
        F = np.zeros((nx, ny))
        for (ix, iy) in modes:
            phi = 2.0 * np.pi * ((ix * np.sqrt(2.0) + iy * np.sqrt(3.0)) % 1.0)
            F += np.cos(2.0 * np.pi * (ix * X / lx + iy * Y / ly) + phi)
        return (amp / np.sqrt(float(len(modes)))) * F

    raise ValueError(f"unknown forcing_type: {ftype!r}")


def build_forcing_operator(nx: int, ny: int, lx: float, ly: float,
                           ftype: str, amp: float,
                           fkx: int, fky: int,
                           fkmin: float, fkmax: float,
                           omega0: float = 0.0):
    """
    v0.2 forcing operator (design doc s7). Returns:
      None                      -- forcing off
      ndarray (nx, ny)          -- static field (v0.1 paths: "mode", or "band" with omega0=0)
      callable t -> ndarray     -- rotating-phase band forcing ("band" with omega0 != 0)

    Rotating band: per-mode phase phi(ix,iy,t) = phi0 + Omega*t with
      phi0  = 2*pi*frac(ix*sqrt(2)  + iy*sqrt(3))
      Omega = 2*pi*frac(ix*sqrt(11) + iy*sqrt(13)) * omega0
    F(t) = (A/sqrt(M)) * sum_m [C_m*cos(Omega_m*t) - S_m*sin(Omega_m*t)],
    C_m/S_m precomputed spatial cos/sin fields. RNG-free, bit-reproducible.
    """
    if ftype != "band" or omega0 == 0.0:
        return build_forcing_field(nx, ny, lx, ly, ftype, amp, fkx, fky, fkmin, fkmax)
    if amp == 0.0:
        return None
    if lx != ly:
        raise ValueError("forcing_type='band' requires a square domain (lx == ly)")

    x = np.arange(nx) * (lx / nx)
    y = np.arange(ny) * (ly / ny)
    X, Y = np.meshgrid(x, y, indexing="ij")

    kmax_i = int(np.floor(fkmax))
    C_list, S_list, Om_list = [], [], []
    for ix in range(0, kmax_i + 1):
        for iy in range(-kmax_i, kmax_i + 1):
            if ix == 0 and iy <= 0:
                continue
            r = np.sqrt(float(ix * ix + iy * iy))
            if fkmin <= r <= fkmax:
                phi0 = 2.0 * np.pi * ((ix * np.sqrt(2.0) + iy * np.sqrt(3.0)) % 1.0)
                om   = 2.0 * np.pi * ((ix * np.sqrt(11.0) + iy * np.sqrt(13.0)) % 1.0) * omega0
                psi  = 2.0 * np.pi * (ix * X / lx + iy * Y / ly) + phi0
                C_list.append(np.cos(psi))
                S_list.append(np.sin(psi))
                Om_list.append(om)
    if not C_list:
        return None
    C = np.stack(C_list); S = np.stack(S_list); Om = np.array(Om_list)
    norm = amp / np.sqrt(float(len(Om_list)))

    def forcing_at(t: float) -> np.ndarray:
        ct = np.cos(Om * t); st = np.sin(Om * t)
        return norm * (np.tensordot(ct, C, axes=1) - np.tensordot(st, S, axes=1))

    return forcing_at


def rk4_step(omega: np.ndarray, nu: float,
             KX: np.ndarray, KY: np.ndarray,
             nx: int, ny: int, dt: float,
             forcing=None, t: float = 0.0, drag: float = 0.0) -> np.ndarray:
    """Single RK4 time step for vorticity field.

    forcing: None | ndarray (steady) | callable t->ndarray (v0.2 rotating band).
    Time-dependent forcing evaluated at RK4 stage times (design doc s7);
    drag adds -drag*omega to every stage RHS (Ekman friction, v0.2).
    """
    if forcing is None:
        f1 = f2 = f4 = None
    elif callable(forcing):
        f1, f2, f4 = forcing(t), forcing(t + 0.5 * dt), forcing(t + dt)
    else:
        f1 = f2 = f4 = forcing

    def step_with_omega(w, f):
        w_hat = np.fft.rfft2(w)
        psi_hat = solve_poisson(w_hat, KX, KY, nx, ny)
        u, v = velocity_from_psi_hat(psi_hat, KX, KY, nx, ny)
        rhs = vorticity_rhs(w, u, v, nu, KX, KY, nx, ny)
        if drag != 0.0:
            rhs = rhs - drag * w
        return rhs if f is None else rhs + f

    k1 = step_with_omega(omega, f1)
    k2 = step_with_omega(omega + 0.5 * dt * k1, f2)
    k3 = step_with_omega(omega + 0.5 * dt * k2, f2)
    k4 = step_with_omega(omega + dt * k3, f4)

    return omega + (dt / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4)


# ── Main simulation entry point ────────────────────────────────────────────────

def simulate(req: SimulationRequest) -> SimulationResult:
    """
    Run KH instability simulation. Returns SimulationResult matching REST schema.
    """
    t_start = time.perf_counter()

    nx, ny   = req.grid_nx, req.grid_ny
    lx, ly   = req.domain_lx, req.domain_ly
    dt       = req.dt
    steps    = req.steps
    # nu = |U| / Re: a reversed shear is the mirrored flow, not negative viscosity
    # (modelling rules v0.1 s2, K3; ruling 2026-10-07). Identical for U >= 0.
    nu       = abs(req.velocity_shear) / req.reynolds_number
    delta    = 0.05 * ly

    dx, dy   = lx / nx, ly / ny
    KX, KY   = _wavenumbers(nx, ny, dx, dy)

    # Initialise vorticity
    if req.initial_omega is not None:
        omega = np.array(req.initial_omega, dtype=float).reshape(nx, ny)
    else:
        omega = initial_conditions(nx, ny, lx, ly,
                                   U0=req.velocity_shear,
                                   delta=delta,
                                   amp=req.perturbation_amplitude,
                                   mode=req.perturbation_mode)

    # Phase-7f forcing (None unless explicitly requested; see design doc s5+s7)
    forcing = build_forcing_operator(nx, ny, lx, ly,
                                     req.forcing_type, req.forcing_amplitude,
                                     req.forcing_kx, req.forcing_ky,
                                     req.forcing_kmin, req.forcing_kmax,
                                     omega0=req.forcing_omega0)

    # CFL guard (KH-NM-5): warn once if user-supplied dt exceeds stability bound
    _omega_hat0 = np.fft.rfft2(omega)
    _psi_hat0   = solve_poisson(_omega_hat0, KX, KY, nx, ny)
    _u0, _v0    = velocity_from_psi_hat(_psi_hat0, KX, KY, nx, ny)
    _dt_cfl     = cfl_dt_bound(float(np.max(np.abs(_u0))), float(np.max(np.abs(_v0))), dx, dy, nu)
    if dt > _dt_cfl:
        warnings.warn(
            f"KH-NM-5: dt={dt:.3e} exceeds CFL bound {_dt_cfl:.3e}; simulation may be unstable",
            RuntimeWarning, stacklevel=2,
        )

    # Time loop — track peak CFL across all steps (KH-NM-5 diagnostic).
    # t starts at t_start: the forcing phase of a continued run (CP-7F-B/D harness).
    t        = float(req.t_start)
    cfl_peak = 0.0
    for _ in range(steps):
        w_hat_step   = np.fft.rfft2(omega)
        psi_hat_step = solve_poisson(w_hat_step, KX, KY, nx, ny)
        u_step, v_step = velocity_from_psi_hat(psi_hat_step, KX, KY, nx, ny)
        cfl_step = dt * (float(np.max(np.abs(u_step))) / dx +
                         float(np.max(np.abs(v_step))) / dy)
        if cfl_step > cfl_peak:
            cfl_peak = cfl_step
        omega = rk4_step(omega, nu, KX, KY, nx, ny, dt,
                         forcing=forcing, t=t, drag=req.linear_drag)
        t += dt

    # Final field recovery
    omega_hat = np.fft.rfft2(omega)
    psi_hat   = solve_poisson(omega_hat, KX, KY, nx, ny)
    psi       = np.fft.irfft2(psi_hat, s=(nx, ny))
    u, v      = velocity_from_psi_hat(psi_hat, KX, KY, nx, ny)

    # Diagnostics
    ke         = 0.5 * float(np.mean(u ** 2 + v ** 2))
    enstrophy  = 0.5 * float(np.mean(omega ** 2))
    max_vort   = float(np.max(np.abs(omega)))
    # Divergence should be machine-zero for spectral method
    div        = float(np.sqrt(np.mean(
        (np.fft.irfft2(1j * KX * np.fft.rfft2(u), s=(nx, ny)) +
         np.fft.irfft2(1j * KY * np.fft.rfft2(v), s=(nx, ny))) ** 2
    )))

    elapsed_ms = (time.perf_counter() - t_start) * 1000.0

    return SimulationResult(
        backend          = "python-numpy",
        language         = "Python",
        steps_completed  = steps,
        t_final          = round(t, 10),
        grid_nx          = nx,
        grid_ny          = ny,
        u_velocity       = u.flatten().tolist(),
        v_velocity       = v.flatten().tolist(),
        vorticity        = omega.flatten().tolist(),
        pressure         = psi.flatten().tolist(),
        diagnostics      = Diagnostics(
            kinetic_energy = ke,
            enstrophy      = enstrophy,
            max_vorticity  = max_vort,
            divergence_rms = div,
            cfl_peak       = cfl_peak,
        ),
        compute_time_ms  = round(elapsed_ms, 2),
        warnings         = request_warnings(dt, dx, dy, nu, cfl_peak),
    )


# ── CLI / smoke test ───────────────────────────────────────────────────────────

def _parse_args():
    p = argparse.ArgumentParser(description="KH instability reference solver")
    p.add_argument("--nx",    type=int,   default=64,    help="Grid points x")
    p.add_argument("--ny",    type=int,   default=32,    help="Grid points y")
    p.add_argument("--lx",    type=float, default=1.0,   help="Domain length x")
    p.add_argument("--ly",    type=float, default=0.5,   help="Domain length y")
    p.add_argument("--dt",    type=float, default=0.001, help="Time step")
    p.add_argument("--steps", type=int,   default=100,   help="Steps to run")
    p.add_argument("--re",    type=float, default=1000.0,help="Reynolds number")
    p.add_argument("--amp",   type=float, default=0.01,  help="Perturbation amplitude")
    p.add_argument("--mode",  type=int,   default=2,     help="Perturbation mode")
    p.add_argument("--out",   type=str,   default="kh_reference_output.json",
                   help="Output JSON path")
    p.add_argument("--full",  action="store_true",
                   help="Include full field arrays in output (large)")
    return p.parse_args()


if __name__ == "__main__":
    args = _parse_args()

    req = SimulationRequest(
        grid_nx=args.nx, grid_ny=args.ny,
        domain_lx=args.lx, domain_ly=args.ly,
        dt=args.dt, steps=args.steps,
        reynolds_number=args.re,
        perturbation_amplitude=args.amp,
        perturbation_mode=args.mode,
    )

    print(f"KH instability reference simulation")
    print(f"  Grid: {req.grid_nx}x{req.grid_ny} | dt={req.dt} | steps={req.steps} | Re={req.reynolds_number}")
    print(f"  Domain: [{req.domain_lx} x {req.domain_ly}] | amp={req.perturbation_amplitude} | mode={req.perturbation_mode}")
    print()

    result = simulate(req)

    print(f"  t_final          = {result.t_final:.4f}")
    print(f"  kinetic_energy   = {result.diagnostics.kinetic_energy:.6f}")
    print(f"  enstrophy        = {result.diagnostics.enstrophy:.6f}")
    print(f"  max_vorticity    = {result.diagnostics.max_vorticity:.6f}")
    print(f"  divergence_rms   = {result.diagnostics.divergence_rms:.2e}  (should be < 1e-10)")
    print(f"  compute_time     = {result.compute_time_ms:.1f} ms")

    # Save reference output (strip large arrays unless --full)
    out_dict = asdict(result)
    if not args.full:
        for fld in ("u_velocity", "v_velocity", "vorticity", "pressure"):
            vals = out_dict[fld]
            out_dict[fld] = {
                "length":   len(vals),
                "first_10": vals[:10],
                "last_10":  vals[-10:],
                "min":      float(min(vals)),
                "max":      float(max(vals)),
                "mean":     float(sum(vals) / len(vals)),
            }

    with open(args.out, "w") as f:
        json.dump(out_dict, f, indent=2)

    print(f"\n  Reference output -> {args.out}")
    print("\n  [PASS] Reference simulation complete.")
