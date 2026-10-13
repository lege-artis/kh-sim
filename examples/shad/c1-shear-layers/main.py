"""main.py -- shad tier, chapter 01 (The billows): three figures and one number.

    fig1  the input: the shear profile U(y) of the layer alone and as the periodic
          box realises it, and the small ripple that starts everything;
    fig2  the growth: the ripple's amplitude against time on a log axis (a straight
          line = exponential growth), next to the vorticity in mid roll-up;
    fig3  the takeaway: the late billows with the stream-function contours drawn
          over them -- the "cat's eyes".

and prints the measured early growth rate next to the textbook value for an
unbounded tanh layer at the same wavenumber.

Run from anywhere (Python 3.10+, NumPy, Matplotlib; about a minute):

    python examples/shad/c1-shear-layers/main.py [--outdir DIR]

Default output directory: docs/shad/figures/ (c1-shear-layers-fig{1,2,3}.png).
Everything is computed by the kh-sim Python reference, shared/physics/kh_physics.py.

Things to try: RE (viscosity), MODE (the ripple's wavelength), AMP. Above
k*delta = 1 the ripple stops growing; mode 7 is the first one there in this box.

License: Apache-2.0 (example code), lege-artis/kh-sim.
ASCII-only source.
"""

import argparse
import math
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[3]          # .../kh-sim
sys.path.insert(0, str(ROOT))
from shared.physics.kh_physics import (  # noqa: E402
    SimulationRequest, simulate, initial_conditions, _wavenumbers,
    solve_poisson, velocity_from_psi_hat)

# The canonical box of the request path: 1.0 x 0.5, layer width delta = Ly/20.
NX, NY = 256, 128
LX, LY = 1.0, 0.5
DELTA = 0.05 * LY
U0 = 1.0
RE = 10000.0
DT = 0.001
MODE = 2                      # two wavelengths across the box
AMP_ROLLUP = 0.01             # the canonical ripple (figures 2b and 3)
AMP_LINEAR = 1e-4             # a ripple small enough to stay linear (figure 2a)
T_MID, T_LATE = 1.2, 2.4
FIT_WINDOW = (0.2, 0.6)       # where the small ripple grows exponentially

# Inviscid growth rate of an unbounded tanh layer at k*delta = pi/10 (this
# box, mode 2), sigma*delta/U0: Rayleigh equation, shooting method of
# shared/physics/rayleigh_eigenvalue.py (0.1763). Michalke (1964): the peak
# is 0.1897 at k*delta = 0.4446.
SIGMA_UNBOUNDED = 0.176


def run(amp, t_marks, chunk=None):
    """Continue one run through t_marks; return {t: (omega, psi, diagnostics)}."""
    out, omega, t = {}, None, 0.0
    for t_target in t_marks:
        steps = int(round((t_target - t) / DT))
        r = simulate(SimulationRequest(grid_nx=NX, grid_ny=NY, domain_lx=LX, domain_ly=LY,
                                       dt=DT, steps=steps, reynolds_number=RE,
                                       velocity_shear=U0,
                                       # a continued run sends the field, no ripple (K4b)
                                       perturbation_amplitude=amp if omega is None else 0.0,
                                       perturbation_mode=MODE, initial_omega=omega,
                                       t_start=t))
        omega, t = r.vorticity, r.t_final
        out[round(t, 6)] = (np.array(omega).reshape(NX, NY),
                            np.array(r.pressure).reshape(NX, NY),   # "pressure" carries psi
                            np.array(r.v_velocity).reshape(NX, NY),
                            r.diagnostics)
    return out


def mode_amplitude(v):
    """Largest |v-hat| of the perturbation mode along x, over all y."""
    return float(np.abs(np.fft.rfft(v, axis=0)[MODE]).max())


def realised_profile():
    """U(y) as the doubly periodic box realises the tanh layer (x-average of u)."""
    w0 = initial_conditions(NX, NY, LX, LY, U0, DELTA, 0.0, MODE)
    kx, ky = _wavenumbers(NX, NY, LX / NX, LY / NY)
    u, _ = velocity_from_psi_hat(solve_poisson(np.fft.rfft2(w0), kx, ky, NX, NY),
                                 kx, ky, NX, NY)
    return u.mean(axis=0)


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--outdir", default=str(ROOT / "docs" / "shad" / "figures"))
    args = ap.parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    y = np.arange(NY) * LY / NY
    x = np.arange(NX) * LX / NX
    k = 2.0 * math.pi * MODE / LX

    # -- the small ripple: exponential growth ------------------------------------
    marks = [round(0.02 * i, 2) for i in range(1, 41)]            # t = 0.02 .. 0.80
    lin = run(AMP_LINEAR, marks)
    ts = np.array(sorted(lin))
    amps = np.array([mode_amplitude(lin[t][2]) for t in ts])
    sel = (ts >= FIT_WINDOW[0] - 1e-9) & (ts <= FIT_WINDOW[1] + 1e-9)
    slope, icpt = np.polyfit(ts[sel], np.log(amps[sel]), 1)
    sigma_nd = slope * DELTA / U0

    # -- the canonical ripple: roll-up ------------------------------------------
    big = run(AMP_ROLLUP, [T_MID, T_LATE])
    w_mid, _, _, d_mid = big[T_MID]
    w_late, psi_late, _, d_late = big[T_LATE]

    # -- figure 1: the input -----------------------------------------------------
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(11, 4), constrained_layout=True)
    a1.plot(-U0 * np.tanh((y - LY / 2) / DELTA), y, "k--", lw=1.2,
            label=r"the layer alone: $-U_0\tanh((y-L_y/2)/\delta)$")
    a1.plot(realised_profile(), y, "C3", lw=2, label="as the periodic box realises it")
    a1.axhspan(LY / 2 - DELTA, LY / 2 + DELTA, color="C0", alpha=0.12,
               label=r"the layer, $\pm\delta$")
    a1.set_xlabel("streamwise velocity U(y)")
    a1.set_ylabel("y")
    a1.set_title("the shear: upper stream to the left, lower to the right", fontsize=10)
    a1.legend(fontsize=8, loc="lower right")
    a2.plot(x, AMP_ROLLUP * np.sin(k * x), "C0", lw=2)
    a2.axhline(0.0, color="0.6", lw=0.8)
    a2.set_xlabel("x")
    a2.set_ylabel("v at t = 0")
    a2.set_title(f"the ripple: v = {AMP_ROLLUP} sin(2 pi {MODE} x / Lx), "
                 f"k delta = {k * DELTA:.3f}", fontsize=10)
    fig.savefig(outdir / "c1-shear-layers-fig1.png", dpi=110)
    plt.close(fig)

    # -- figure 2: the growth and the mid roll-up -----------------------------------
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(12, 4.2), constrained_layout=True,
                                 gridspec_kw={"width_ratios": [1, 1.6]})
    a1.semilogy(ts, amps, "o", ms=3, color="C0", label="kh-sim, amplitude of mode 2")
    tf = np.linspace(0.0, 0.8, 50)
    a1.semilogy(tf, np.exp(icpt + slope * tf), "C0", lw=1,
                label=fr"fit {FIT_WINDOW[0]}-{FIT_WINDOW[1]}: $\sigma\delta/U_0$ = {sigma_nd:.3f}")
    a1.semilogy(tf, np.exp(icpt + SIGMA_UNBOUNDED * U0 / DELTA * (tf - 0.4)
                           + slope * 0.4), "k--", lw=1,
                label=fr"unbounded layer: {SIGMA_UNBOUNDED}")
    a1.set_xlabel("t")
    a1.set_title("a small ripple grows exponentially", fontsize=10)
    a1.legend(fontsize=8)
    im = a2.imshow(w_mid.T, origin="lower", extent=(0, LX, 0, LY), cmap="RdBu_r",
                   vmin=-40, vmax=40, interpolation="bilinear")
    a2.set_title(f"the canonical ripple rolls up, t = {T_MID}", fontsize=10)
    a2.set_xlabel("x")
    a2.set_ylabel("y")
    fig.colorbar(im, ax=a2, shrink=0.85, label=r"vorticity $\omega$")
    fig.savefig(outdir / "c1-shear-layers-fig2.png", dpi=110)
    plt.close(fig)

    # -- figure 3: the takeaway ----------------------------------------------------
    fig, ax = plt.subplots(figsize=(9, 4.5), constrained_layout=True)
    im = ax.imshow(w_late.T, origin="lower", extent=(0, LX, 0, LY), cmap="RdBu_r",
                   vmin=-40, vmax=40, interpolation="bilinear")
    ax.contour(x, y, psi_late.T, levels=24, colors="k", linewidths=0.6,
               linestyles="solid")
    ax.set_title(f"t = {T_LATE}: two billows, the stream function's closed "
                 "'cat's eyes' around them", fontsize=10)
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    fig.colorbar(im, ax=ax, shrink=0.8, label=r"vorticity $\omega$")
    fig.savefig(outdir / "c1-shear-layers-fig3.png", dpi=90)
    plt.close(fig)

    print(f"k*delta                      = {k * DELTA:.4f}")
    print(f"measured sigma*delta/U0      = {sigma_nd:.3f}  (t = {FIT_WINDOW[0]} .. {FIT_WINDOW[1]}, Re = {RE:.0f})")
    print(f"unbounded tanh layer         = {SIGMA_UNBOUNDED:.3f}")
    print(f"ratio                        = {sigma_nd / SIGMA_UNBOUNDED:.2f}")
    for t, d in ((T_MID, d_mid), (T_LATE, d_late)):
        print(f"t = {t:<4}  kinetic_energy = {d.kinetic_energy:.4f}  enstrophy = {d.enstrophy:.2f}  "
              f"max_vorticity = {d.max_vorticity:.2f}  divergence_rms < 1e-10: {d.divergence_rms < 1e-10}")
    print(f"wrote {outdir}/c1-shear-layers-fig{{1,2,3}}.png")


if __name__ == "__main__":
    main()
