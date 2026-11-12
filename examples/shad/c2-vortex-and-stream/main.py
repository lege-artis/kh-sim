"""main.py -- shad tier, chapter 02 (The vortex and the stream): three figures, two numbers.

    fig1  one vortex: its vorticity, the stream function's circles around it, and
          the swirl speed across it against the formula (Lamb-Oseen), alone and as the
          periodic box realises it (no net circulation in the box);
    fig2  viscosity at work: the vortex's peak vorticity against time, next to the
          exact answer for a spreading Gaussian vortex;
    fig3  two vortices of the same sign orbit each other: the vorticity is carried,
          the stream function drawn in the frame that turns with them.

It prints the peak vorticity against the exact Lamb-Oseen value at five times, and
the orbiting rate of the pair against two point vortices of the same strength.

Run from anywhere (Python 3.10+, NumPy, Matplotlib; about a minute):

    python examples/shad/c2-vortex-and-stream/main.py [--outdir DIR]

Default output directory: docs/shad/figures/ (c2-vortex-and-stream-fig{1,2,3}.png).
Everything is computed by the kh-sim Python reference, shared/physics/kh_physics.py.

Things to try: RE (viscosity), CORE (vortex size), D_PAIR (pair separation; below
about 3.3 core radii the two merge into one).

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
    SimulationRequest, simulate, _wavenumbers, solve_poisson, velocity_from_psi_hat)

N = 128                       # grid, both directions
L = 2.0 * math.pi             # box, both directions
GAMMA = 2.0 * math.pi         # circulation of each vortex
RE = 1000.0                   # nu = 1 / RE (velocity_shear 1)
NU = 1.0 / RE
DT = 0.002

CORE = 0.3                    # the single vortex's core radius a
T_DECAY = (0.0, 2.5, 5.0, 7.5, 10.0)

CORE_PAIR = 0.2               # each vortex of the pair
D_PAIR = 1.2                  # their separation (6 core radii: they orbit, no merger)
T_PAIR = 4.5                  # about one orbit
T_SNAP = 2.25                 # about half an orbit (figure 3)
SAMPLE = 0.25                 # how often the pair's positions are read

X = np.arange(N) * L / N
XX, YY = np.meshgrid(X, X, indexing="ij")
C = L / 2.0                   # the box centre, a grid point


def gaussian(x0, y0, a):
    """A Gaussian vortex of circulation GAMMA and core radius a at (x0, y0),
    with its periodic images (so the field is smooth across the box edges)."""
    w = np.zeros((N, N))
    for i in (-1, 0, 1):
        for j in (-1, 0, 1):
            r2 = (XX - x0 - i * L) ** 2 + (YY - y0 - j * L) ** 2
            w += GAMMA / (math.pi * a * a) * np.exp(-r2 / (a * a))
    return w


def advance(omega, t, t_end):
    r = simulate(SimulationRequest(grid_nx=N, grid_ny=N, domain_lx=L, domain_ly=L, dt=DT,
                                   steps=int(round((t_end - t) / DT)), reynolds_number=RE,
                                   velocity_shear=1.0, perturbation_amplitude=0.0,
                                   initial_omega=omega.ravel().tolist(), t_start=t))
    return (np.array(r.vorticity).reshape(N, N), np.array(r.pressure).reshape(N, N),
            r.t_final, r.diagnostics)


def stream_and_velocity(omega):
    kx, ky = _wavenumbers(N, N, L / N, L / N)
    psi_hat = solve_poisson(np.fft.rfft2(omega), kx, ky, N, N)
    u, v = velocity_from_psi_hat(psi_hat, kx, ky, N, N)
    return np.fft.irfft2(psi_hat, s=(N, N)), u, v


def peaks(omega):
    """The two largest local maxima of omega (periodic), as (x, y) positions refined
    by a parabola through the neighbours."""
    nb = [np.roll(np.roll(omega, i, 0), j, 1) for i in (-1, 0, 1) for j in (-1, 0, 1)
          if (i, j) != (0, 0)]
    is_max = np.all([omega > m for m in nb], axis=0)
    idx = np.argwhere(is_max)
    idx = idx[np.argsort(-omega[is_max])][:2]
    out = []
    for i, j in idx:
        f = lambda a, b: omega[(i + a) % N, (j + b) % N]  # noqa: E731
        dx = 0.5 * (f(-1, 0) - f(1, 0)) / (f(-1, 0) - 2 * f(0, 0) + f(1, 0))
        dy = 0.5 * (f(0, -1) - f(0, 1)) / (f(0, -1) - 2 * f(0, 0) + f(0, 1))
        out.append(((i + dx) * L / N, (j + dy) * L / N))
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--outdir", default=str(ROOT / "docs" / "shad" / "figures"))
    args = ap.parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    # -- one vortex: spreading by viscosity --------------------------------------------
    w = gaussian(C, C, CORE)
    psi0, u0, v0 = stream_and_velocity(w)
    tt, pk, t = [0.0], [float(w.max())], 0.0
    for t_mark in np.arange(0.5, T_DECAY[-1] + 1e-9, 0.5):
        w, _, t, _ = advance(w, t, float(t_mark))
        tt.append(float(t_mark))
        pk.append(float(w.max()))
    series = (np.array(tt), np.array(pk))
    rows = [(t_mark, pk[tt.index(t_mark)], GAMMA / (math.pi * (CORE ** 2 + 4.0 * NU * t_mark)))
            for t_mark in T_DECAY]

    # -- the pair: carried by each other's flow ------------------------------------------
    h = D_PAIR / 2.0
    wp = gaussian(C - h, C, CORE_PAIR) + gaussian(C + h, C, CORE_PAIR)
    t, times, angles, seps, snap = 0.0, [], [], [], None
    for t_mark in np.arange(0.0, T_PAIR + 1e-9, SAMPLE):
        if t_mark > t:
            wp, psi_p, t, d_pair = advance(wp, t, float(t_mark))
        (x1, y1), (x2, y2) = peaks(wp)
        times.append(t_mark)
        angles.append(math.atan2(y2 - y1, x2 - x1))
        seps.append(math.hypot(x2 - x1, y2 - y1))
        if abs(t_mark - T_SNAP) < 1e-9:
            snap = (wp.copy(), [(x1, y1), (x2, y2)])
    # The two vortices are identical, so the line through them is defined up to a half
    # turn: unwrap twice the angle.
    omega_rate = float(np.polyfit(times, np.unwrap(2.0 * np.array(angles)) / 2.0, 1)[0])
    omega_points = GAMMA / (math.pi * D_PAIR ** 2)
    # The box holds no net circulation: the solver removes the mean vorticity, so the
    # pair sits in a uniform background of -2 GAMMA / L^2, which turns it back at
    # GAMMA / L^2. With the separation as measured (viscosity draws the pair in slowly):
    omega_box = float(np.mean([GAMMA / (math.pi * s * s) for s in seps])) - GAMMA / L ** 2

    # -- figure 1: one vortex ------------------------------------------------------------
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(11, 4.4), constrained_layout=True)
    im = a1.imshow(gaussian(C, C, CORE).T, origin="lower", extent=(0, L, 0, L),
                   cmap="RdBu_r", vmin=-25, vmax=25, interpolation="bilinear")
    a1.contour(X, X, psi0.T, levels=16, colors="k", linewidths=0.6, linestyles="solid")
    a1.set_xlim(C - 1.5, C + 1.5)
    a1.set_ylim(C - 1.5, C + 1.5)
    a1.set_title("vorticity (colour) and stream function (lines)", fontsize=10)
    a1.set_xlabel("x")
    a1.set_ylabel("y")
    fig.colorbar(im, ax=a1, shrink=0.85, label=r"vorticity $\omega$")
    j0 = int(round(C * N / L))
    r = X[j0:] - C
    a2.plot(r, v0[j0:, j0], "o", ms=3, color="C0", label="kh-sim, from the stream function")
    rr = np.linspace(1e-3, r.max(), 300)
    a2.plot(rr, GAMMA / (2 * math.pi * rr) * (1 - np.exp(-rr ** 2 / CORE ** 2)), "k-", lw=1,
            label=r"Lamb-Oseen: $\Gamma/(2\pi r)\,(1-e^{-r^2/a^2})$")
    a2.plot(rr, GAMMA / (2 * math.pi * rr) * (1 - np.exp(-rr ** 2 / CORE ** 2))
            - GAMMA * rr / (2 * L * L), "C3--", lw=1,
            label=r"the same in the box: $-\,\Gamma r/(2L^2)$")
    a2.set_xlabel("distance from the centre r")
    a2.set_ylabel("swirl speed")
    a2.set_title("the swirl speed: solid-body inside, 1/r outside, less the box's share", fontsize=10)
    a2.legend(fontsize=8)
    fig.savefig(outdir / "c2-vortex-and-stream-fig1.png", dpi=100)
    plt.close(fig)

    # -- figure 2: spreading ---------------------------------------------------------------
    fig, ax = plt.subplots(figsize=(7, 4.2), constrained_layout=True)
    ts, pks = series
    tf = np.linspace(0, 10, 200)
    ax.plot(tf, GAMMA / (math.pi * (CORE ** 2 + 4 * NU * tf)), "k-", lw=1,
            label=r"exact: $\Gamma/(\pi(a^2+4\nu t))$")
    ax.plot(ts, pks, "o", ms=4, color="C3", label="kh-sim, peak vorticity")
    ax.set_xlabel("t")
    ax.set_ylabel("peak vorticity")
    ax.set_title(f"viscosity spreads the vortex (Re = {RE:.0f})", fontsize=10)
    ax.legend(fontsize=8)
    fig.savefig(outdir / "c2-vortex-and-stream-fig2.png", dpi=100)
    plt.close(fig)

    # -- figure 3: the pair ------------------------------------------------------------------
    w_s, pos = snap
    psi_s, _, _ = stream_and_velocity(w_s)
    psi_rot = psi_s + 0.5 * omega_rate * ((XX - C) ** 2 + (YY - C) ** 2)
    fig, ax = plt.subplots(figsize=(6.4, 5.4), constrained_layout=True)
    im = ax.imshow(w_s.T, origin="lower", extent=(0, L, 0, L), cmap="RdBu_r",
                   vmin=-40, vmax=40, interpolation="bilinear")
    ax.contour(X, X, psi_rot.T, levels=40, colors="k", linewidths=0.6, linestyles="solid")
    ax.plot([p[0] for p in pos], [p[1] for p in pos], "k+", ms=10)
    ax.set_xlim(C - 1.8, C + 1.8)
    ax.set_ylim(C - 1.8, C + 1.8)
    ax.set_title(f"t = {T_SNAP}: the pair has turned half way;\nstream function in the "
                 "frame that turns with it", fontsize=10)
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    fig.colorbar(im, ax=ax, shrink=0.85, label=r"vorticity $\omega$")
    fig.savefig(outdir / "c2-vortex-and-stream-fig3.png", dpi=90)
    plt.close(fig)

    print(f"one vortex: circulation {GAMMA:.4f}, core a = {CORE}, Re = {RE:.0f} (nu = {NU})")
    for t_mark, got, exact in rows:
        print(f"t = {t_mark:4.1f}  peak vorticity = {got:8.4f}  Lamb-Oseen = {exact:8.4f}  "
              f"ratio = {got / exact:.4f}")
    print(f"the pair: separation {D_PAIR}, cores {CORE_PAIR}, t = 0 .. {T_PAIR}")
    print(f"separation                   = {seps[0]:.3f} .. {seps[-1]:.3f}")
    print(f"two point vortices, alone    = {omega_points:.3f}  (Gamma / (pi d^2), d = {D_PAIR})")
    print(f"the same in the box          = {omega_box:.3f}  (minus Gamma / L^2, d as measured)")
    print(f"measured orbit rate          = {omega_rate:.3f}")
    print(f"ratio measured / box         = {omega_rate / omega_box:.3f}")
    print(f"divergence_rms < 1e-10: {d_pair.divergence_rms < 1e-10}")
    print(f"wrote {outdir}/c2-vortex-and-stream-fig{{1,2,3}}.png")


if __name__ == "__main__":
    main()
