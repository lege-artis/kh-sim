"""main.py -- shad tier, chapter 03 (Reading the growth rate): three figures, one curve.

    fig1  the set-up: the two shear layers a periodic box needs, and the seed that
          starts a wavelength (the box's own unstable mode, made small);
    fig2  the measurement: the amplitude of each seeded wavelength against time on a
          log axis -- each a straight line, with its own slope;
    fig3  the curve: growth rate against wavenumber, kh-sim's dots on the Rayleigh
          equation's line (Michalke 1964), and the viscous correction removed by
          measuring one wavelength at three Reynolds numbers.

It prints the growth rate measured at k*delta = 0.1 .. 1.1 next to the inviscid
theory, and the inviscid limit extrapolated at k*delta = 0.4.

Run from anywhere (Python 3.10+, NumPy, SciPy, Matplotlib; about two minutes):

    python examples/shad/c3-growth-rate/main.py [--outdir DIR]

Default output directory: docs/shad/figures/ (c3-growth-rate-fig{1,2,3}.png).
Everything is computed by the kh-sim Python reference, shared/physics/kh_physics.py,
and the Rayleigh solver next to it, shared/physics/rayleigh_eigenvalue.py.

Things to try: RE (viscosity), NY (resolution across the layers), T_END.

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
sys.path.insert(0, str(ROOT / "shared" / "physics"))
from shared.physics.kh_physics import SimulationRequest, simulate  # noqa: E402
import rayleigh_eigenvalue  # noqa: E402

# The box of kh-sim's testbed A (engineer tier, chapter 01 s6), at a quarter of its
# resolution: Lx = pi^2 makes mode m sit at k*delta = m/10 exactly.
LX, LY = math.pi ** 2, 2.0 * math.pi
NX, NY = 128, 256              # y finer: the layers are thin, and near the cut-off the
                              # weakly unstable modes need their critical layer resolved
DELTA = LY / 40.0
RE = 40000.0
DT = 0.005
T_END = 4.0
CHUNKS = 8                    # amplitude read every 0.5 time units
MODES = list(range(1, 12))    # k*delta = 0.1 .. 1.1
SEED = 1e-6
NEUTRAL = 1e-6                # growth rate * delta below which the box has no growing mode
EXTRAP_MODE = 4
EXTRAP_RE = (10000.0, 20000.0, 40000.0)
SHOWN = (2, 4, 8, 10)         # the growth histories drawn in figure 2

X = np.arange(NX) * LX / NX
Y = np.arange(NY) * LY / NY
XX, YY = np.meshgrid(X, Y, indexing="ij")
ETA1 = (YY - LY / 4) / DELTA
ETA2 = (YY - 3 * LY / 4) / DELTA


def base_flow():
    """Two tanh layers, the second the mirror of the first, so the velocity is periodic
    in y: U = tanh(eta1) - tanh(eta2) - 1. Vorticity omega = -dU/dy."""
    u = np.tanh(ETA1[0]) - np.tanh(ETA2[0]) - 1.0
    w = -(1.0 / np.cosh(ETA1) ** 2 - 1.0 / np.cosh(ETA2) ** 2) / DELTA
    return u, w


def own_mode(mode):
    """The box's own inviscid linear modes for this wavelength (the Rayleigh operator
    (U - c)(D^2 - k^2) psi = U'' psi on the NY grid points, by Fourier collocation):
    the fastest growth rate times delta, and its vorticity profile w(y), max |w| = 1.
    The same construction as kh-sim's testbed A."""
    k = 2.0 * math.pi * mode / LX
    u, _ = base_flow()
    ky = np.fft.fftfreq(NY, 1.0 / NY) * 2.0 * math.pi / LY
    d2 = np.fft.ifft(np.diag(-(ky ** 2)) @ np.fft.fft(np.eye(NY), axis=0), axis=0)
    lap = d2 - k * k * np.eye(NY)
    a = -1j * k * np.diag(u) @ lap + 1j * k * np.diag(np.real(d2 @ u))
    lam, vec = np.linalg.eig(np.linalg.solve(lap, a))
    i = int(np.argmax(lam.real))
    if lam[i].real * DELTA < NEUTRAL:
        # At and past the cut-off every mode is neutral: the "fastest" one is round-off,
        # and which neutral mode wins differs from one LAPACK build to the next.
        return 0.0, None
    w = lap @ vec[:, i]
    return float(lam[i].real) * DELTA, w / w[np.argmax(np.abs(w))]


def wobble():
    """The seed where the box has no growing shape: vorticity on the two layers,
    sech^2 across each, the same for every wavelength."""
    return 1.0 / np.cosh(ETA1[0]) ** 2 + 1.0 / np.cosh(ETA2[0]) ** 2


def measure(mode, re):
    """Seed one wavelength, run to T_END, return (times, amplitudes, growth rate)."""
    _, w = base_flow()
    k = 2.0 * math.pi * mode / LX
    shape = own_mode(mode)[1]
    if shape is None:
        shape = wobble()
    w = w + SEED * np.real(shape[None, :] * np.exp(1j * k * XX))
    t, ts, amps = 0.0, [], []
    for c in range(1, CHUNKS + 1):
        t_end = c * T_END / CHUNKS
        r = simulate(SimulationRequest(grid_nx=NX, grid_ny=NY, domain_lx=LX, domain_ly=LY,
                                       dt=DT, steps=int(round((t_end - t) / DT)),
                                       reynolds_number=re, velocity_shear=1.0,
                                       perturbation_amplitude=0.0,
                                       initial_omega=w.ravel().tolist(), t_start=t))
        w, t = np.array(r.vorticity).reshape(NX, NY), r.t_final
        ts.append(t)
        amps.append(float(np.abs(np.fft.fft2(w)[mode, :]).max()))
    ts, amps = np.array(ts), np.array(amps)
    rate = float(np.polyfit(ts, np.log(amps), 1)[0])
    return ts, amps, rate


def rayleigh(kd):
    """Inviscid growth rate sigma*delta/U of a single tanh layer (Rayleigh equation,
    Chebyshev collocation), walls far enough away not to matter at this k*delta."""
    if kd >= 1.0:
        return 0.0                               # the cut-off: no unstable mode
    # Long waves feel walls from afar; short ones need the critical layer resolved.
    ly, n = (40.0, 200) if kd < 0.6 else (16.0, 400)
    return rayleigh_eigenvalue.rayleigh_cheb(kd, ly, ly / 2, 1.0, N=n)[0]


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--outdir", default=str(ROOT / "docs" / "shad" / "figures"))
    args = ap.parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    runs = {m: measure(m, RE) for m in MODES}
    sig = {m: runs[m][2] * DELTA for m in MODES}
    theory = {m: rayleigh(m / 10.0) for m in MODES}
    box = {m: own_mode(m)[0] for m in MODES}
    ext = [(re, measure(EXTRAP_MODE, re)[2] * DELTA) for re in EXTRAP_RE if re != RE]
    ext.append((RE, sig[EXTRAP_MODE]))
    ext.sort()
    slope, sigma_inf = np.polyfit([1.0 / re for re, _ in ext], [s for _, s in ext], 1)

    # -- figure 1: the set-up ----------------------------------------------------------
    u, w = base_flow()
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(11, 4), constrained_layout=True)
    a1.plot(u, Y, "C3", lw=2)
    for yc in (LY / 4, 3 * LY / 4):
        a1.axhspan(yc - DELTA, yc + DELTA, color="C0", alpha=0.12)
    a1.set_xlabel("streamwise velocity U(y)")
    a1.set_ylabel("y")
    a1.set_title("two layers, so that the flow can repeat in y", fontsize=10)
    m = EXTRAP_MODE
    seed = np.real(own_mode(m)[1][None, :] * np.exp(2j * math.pi * m * XX / LX))
    im = a2.imshow(seed.T, origin="lower", extent=(0, LX, 0, LY), cmap="RdBu_r",
                   aspect="auto", interpolation="bilinear", vmin=-1, vmax=1)
    a2.set_title(f"the seed for k delta = {m / 10:.1f}: the box's own unstable mode", fontsize=10)
    a2.set_xlabel("x")
    fig.colorbar(im, ax=a2, shrink=0.85, label=r"vorticity added (x %g)" % SEED)
    fig.savefig(outdir / "c3-growth-rate-fig1.png", dpi=100)
    plt.close(fig)

    # -- figure 2: the measurement --------------------------------------------------------
    fig, ax = plt.subplots(figsize=(7.5, 4.4), constrained_layout=True)
    for i, m in enumerate(SHOWN):
        ts, amps, rate = runs[m]
        ax.semilogy(ts, amps / amps[0], "o-", ms=3, lw=1, color=f"C{i}",
                    label=fr"k$\delta$ = {m / 10:.1f}: $\sigma\delta/U$ = {rate * DELTA:.3f}")
    ax.set_xlabel("t")
    ax.set_ylabel("amplitude of the seeded wavelength (relative)")
    ax.set_title(f"each wavelength grows at its own rate (Re = {RE:.0f})", fontsize=10)
    ax.legend(fontsize=8)
    fig.savefig(outdir / "c3-growth-rate-fig2.png", dpi=100)
    plt.close(fig)

    # -- figure 3: the curve ----------------------------------------------------------------
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(11.5, 4.3), constrained_layout=True,
                                 gridspec_kw={"width_ratios": [1.6, 1]})
    kd = np.r_[np.linspace(0.05, 0.95, 37), 1.0, 1.1]
    a1.plot(kd, [rayleigh(x) for x in kd], "k-", lw=1,
            label="Rayleigh equation, inviscid (Michalke 1964)")
    a1.plot([m / 10 for m in MODES], [box[m] for m in MODES], "s", ms=7, mfc="none",
            color="k", label="the box's own modes, inviscid")
    a1.plot([m / 10 for m in MODES], [sig[m] for m in MODES], "o", ms=5, color="C3",
            label=f"kh-sim, Re = {RE:.0f}")
    a1.axhline(0.0, color="0.6", lw=0.8)
    a1.axvline(1.0, color="0.6", lw=0.8, ls=":")
    a1.annotate("peak 0.1897\nat 0.4446", (0.4446, 0.1897), (0.08, 0.17), fontsize=8,
                arrowprops=dict(arrowstyle="-", lw=0.6))
    a1.annotate("cut-off:\nk delta = 1", (1.0, 0.0), (1.02, 0.06), fontsize=8,
                arrowprops=dict(arrowstyle="-", lw=0.6))
    a1.set_xlim(0.0, 1.2)
    a1.set_xlabel(r"k$\delta$")
    a1.set_ylabel(r"growth rate $\sigma\delta/U$")
    a1.set_title("the dispersion relation", fontsize=10)
    a1.legend(fontsize=8)
    inv = np.array([1.0 / re for re, _ in ext])
    a2.plot(inv * 1e5, [s for _, s in ext], "o", color="C3", ms=5, label="kh-sim")
    xx = np.linspace(0, inv.max() * 1.05, 20)
    a2.plot(xx * 1e5, sigma_inf + slope * xx, "C3-", lw=1,
            label=fr"fit: $\sigma_\infty\delta/U$ = {sigma_inf:.4f}")
    a2.axhline(box[EXTRAP_MODE], color="k", lw=1, ls="--",
               label=f"the box's own mode, inviscid {box[EXTRAP_MODE]:.4f}")
    a2.set_xlabel(r"1/Re  ($\times 10^{-5}$)")
    a2.set_title(fr"k$\delta$ = {EXTRAP_MODE / 10:.1f}: viscosity removed by extrapolation",
                 fontsize=10)
    a2.legend(fontsize=8)
    fig.savefig(outdir / "c3-growth-rate-fig3.png", dpi=100)
    plt.close(fig)

    print(f"box {LX:.4f} x {LY:.4f}, {NX} x {NY}, delta = Ly/40 = {DELTA:.4f}, Re = {RE:.0f}")
    print("k*delta   kh-sim   box, inviscid   single layer, inviscid")
    for m in MODES:
        print(f"  {m / 10:.1f}     {sig[m]:7.4f}     {box[m]:7.4f}        {theory[m]:7.4f}")
    print(f"k*delta = {EXTRAP_MODE / 10:.1f} at Re = " + ", ".join(f"{re:.0f}" for re, _ in ext)
          + ": " + ", ".join(f"{s:.4f}" for _, s in ext))
    print(f"extrapolated to 1/Re = 0       = {sigma_inf:.4f}")
    print(f"the box's own mode, inviscid   = {box[EXTRAP_MODE]:.4f}")
    print(f"single layer (Michalke)        = {theory[EXTRAP_MODE]:.4f}")
    print(f"wrote {outdir}/c3-growth-rate-fig{{1,2,3}}.png")


if __name__ == "__main__":
    main()
