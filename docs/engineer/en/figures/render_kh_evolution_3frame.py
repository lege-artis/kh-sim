"""render_kh_evolution_3frame.py -- figure for engineer tier, chapter 01.

Three vorticity snapshots of the canonical Kelvin-Helmholtz run of the Python
reference (shared/physics/kh_physics.py): the shear layer at t = 0, the roll-up
at t = 0.8 and the saturated billows at t = 1.6.

Run from anywhere (Python 3.10+, NumPy, Matplotlib; about 30 s):

    python docs/engineer/en/figures/render_kh_evolution_3frame.py [--out FILE.png]

Default output: kh-evolution-3frame.png next to this script.
The run is continued in chunks: each call to simulate() starts from the field
the previous one returned (initial_omega + t_start), which is how any client of
kh-sim takes snapshots of one run.

License: Apache-2.0 (example code), lege-artis/kh-sim.
ASCII-only source.
"""

import argparse
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[4]          # .../kh-sim
sys.path.insert(0, str(ROOT))
from shared.physics.kh_physics import SimulationRequest, simulate  # noqa: E402

NX, NY = 256, 128            # grid; the domain is the default 1.0 x 0.5
RE = 10000.0                 # nu = U0 / Re = 1e-4
DT = 0.001
AMP, MODE = 0.01, 2          # the canonical ripple: v = AMP sin(2 pi MODE x / Lx)
FRAMES = (0.0, 0.8, 1.6)     # snapshot times


def run_frames():
    """Return [(t, omega[nx, ny], diagnostics)] at the FRAMES times."""
    out = []
    omega, t = None, 0.0
    for t_target in FRAMES:
        steps = int(round((t_target - t) / DT))
        # A continued run sends the field and no ripple: the HTTP backends refuse
        # perturbation_amplitude != 0 together with initial_omega (K4b).
        amp = AMP if omega is None else 0.0
        r = simulate(SimulationRequest(grid_nx=NX, grid_ny=NY, dt=DT, steps=steps,
                                       reynolds_number=RE, perturbation_amplitude=amp,
                                       perturbation_mode=MODE, initial_omega=omega,
                                       t_start=t))
        omega, t = r.vorticity, r.t_final
        out.append((t, np.array(omega).reshape(NX, NY), r.diagnostics))
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--out", default=str(Path(__file__).with_name("kh-evolution-3frame.png")))
    args = ap.parse_args()

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    frames = run_frames()
    fig, axes = plt.subplots(1, 3, figsize=(12, 4), constrained_layout=True)
    vmax = 40.0                                      # U0 / delta of the initial layer
    for ax, (t, w, d) in zip(axes, frames):
        im = ax.imshow(w.T, origin="lower", extent=(0.0, 1.0, 0.0, 0.5),
                       cmap="RdBu_r", vmin=-vmax, vmax=vmax, interpolation="bilinear")
        ax.set_title(f"t = {t:.1f}   max|$\\omega$| = {d.max_vorticity:.1f}", fontsize=10)
        ax.set_xlabel("x")
        print(f"t = {t:.1f}  kinetic_energy = {d.kinetic_energy:.4f}  "
              f"enstrophy = {d.enstrophy:.2f}  max_vorticity = {d.max_vorticity:.2f}")
    axes[0].set_ylabel("y")
    fig.colorbar(im, ax=axes, shrink=0.8, label=r"vorticity $\omega$")
    fig.suptitle(f"kh-sim Python reference: {NX} x {NY}, Re = {RE:.0f}, "
                 "perturbation mode 2, amplitude 0.01", fontsize=11)
    fig.savefig(args.out, dpi=120)
    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
