"""main.py -- shad tier, chapter 04 (Why six languages): two figures, two bugs.

Runs chapter 1's billows three times with the kh-sim Python reference:

    correct         as shipped;
    conjugate FFT   the Fourier transform with the opposite sign convention, the bug
                    the C++, Fortran and Pascal backends once shared;
    no de-aliasing  the 2/3 rule left out of the nonlinear term, the bug the Rust
                    backend once had;

and compares each wrong run with the right one twice: by the summary numbers a
solver reports (kinetic energy, enstrophy, peak vorticity, peak Courant number) and
point by point, field against field.

    fig1  the velocity of the right run and of the conjugate-FFT run, side by side;
    fig2  how far each wrong run is from the right one, by summary and by field.

Run from anywhere (Python 3.10+, NumPy, Matplotlib; about twenty seconds):

    python examples/shad/c4-six-languages/main.py [--outdir DIR]

Default output directory: docs/shad/figures/ (c4-six-languages-fig{1,2}.png).
The bugs are introduced into this script's copy of the reference only, by
replacing two of its building blocks for the length of one run.

License: Apache-2.0 (example code), lege-artis/kh-sim.
ASCII-only source.
"""

import argparse
import sys
from contextlib import contextmanager
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[3]          # .../kh-sim
sys.path.insert(0, str(ROOT))
import shared.physics.kh_physics as kh  # noqa: E402

NX, NY = 128, 64
LX, LY = 1.0, 0.5
STEPS = 1200                  # t = 1.2, chapter 1's mid roll-up
# divergence_rms is left out: it is round-off (1e-16) in every run, and the relative
# difference of two round-off numbers says nothing.
SUMMARY = ("kinetic_energy", "enstrophy", "max_vorticity", "cfl_peak")
FIELDS = ("vorticity", "u_velocity", "v_velocity")


class _ConjugateFFT:
    """numpy.fft with the opposite sign convention: forward exp(+2 pi i kn/N), inverse
    exp(-2 pi i kn/N). A consistent pair, so a transform followed by its inverse is
    still the identity -- which is why the mistake survives every round-trip test."""
    def __getattr__(self, name):
        return getattr(np.fft, name)

    def rfft2(self, a, *args, **kw):
        return np.conj(np.fft.rfft2(a, *args, **kw))

    def irfft2(self, a, *args, **kw):
        return np.fft.irfft2(np.conj(a), *args, **kw)


class _NumpyWithConjugateFFT:
    fft = _ConjugateFFT()

    def __getattr__(self, name):
        return getattr(np, name)


@contextmanager
def conjugate_fft():
    saved = kh.np
    kh.np = _NumpyWithConjugateFFT()
    try:
        yield
    finally:
        kh.np = saved


@contextmanager
def no_dealiasing():
    saved = kh.dealias_two_thirds
    kh.dealias_two_thirds = lambda field_hat, kx, ky: field_hat
    try:
        yield
    finally:
        kh.dealias_two_thirds = saved


@contextmanager
def as_shipped():
    yield


def run(variant):
    with variant():
        r = kh.simulate(kh.SimulationRequest(grid_nx=NX, grid_ny=NY, domain_lx=LX, domain_ly=LY,
                                             dt=0.001, steps=STEPS, reynolds_number=10000.0,
                                             velocity_shear=1.0, perturbation_amplitude=0.01,
                                             perturbation_mode=2))
    out = {k: getattr(r.diagnostics, k) for k in SUMMARY}
    out.update({f: np.array(getattr(r, f)).reshape(NX, NY) for f in FIELDS})
    return out


def rel(a, b):
    """Relative difference of two numbers, or of two fields (largest point difference
    over the largest value)."""
    if isinstance(a, np.ndarray):
        return float(np.abs(a - b).max() / np.abs(a).max())
    return abs(a - b) / abs(a)


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--outdir", default=str(ROOT / "docs" / "shad" / "figures"))
    args = ap.parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    good = run(as_shipped)
    wrong = {"conjugate FFT": run(conjugate_fft), "no de-aliasing": run(no_dealiasing)}

    # -- figure 1: the same energy, the opposite flow --------------------------------------
    fig, axes = plt.subplots(1, 2, figsize=(12, 3.6), constrained_layout=True)
    for ax, (name, res) in zip(axes, (("correct", good), ("conjugate FFT", wrong["conjugate FFT"]))):
        im = ax.imshow(res["u_velocity"].T, origin="lower", extent=(0, LX, 0, LY), cmap="PuOr",
                       vmin=-1.2, vmax=1.2, interpolation="bilinear")
        ax.set_title(f"{name}: streamwise velocity u at t = 1.2\n"
                     f"kinetic energy {res['kinetic_energy']:.10f}", fontsize=9)
        ax.set_xlabel("x")
    axes[0].set_ylabel("y")
    fig.colorbar(im, ax=axes, shrink=0.85, label="u")
    fig.savefig(outdir / "c4-six-languages-fig1.png", dpi=90)
    plt.close(fig)

    # -- figure 2: what each comparison sees --------------------------------------------
    labels = list(SUMMARY) + list(FIELDS)
    fig, ax = plt.subplots(figsize=(10, 4.2), constrained_layout=True)
    floor = 1e-17
    width = 0.38
    for i, (name, res) in enumerate(wrong.items()):
        vals = [max(rel(good[k], res[k]), floor) for k in labels]
        ax.bar(np.arange(len(labels)) + (i - 0.5) * width, vals, width, label=name,
               color=f"C{i}")
    ax.set_yscale("log")
    ax.set_ylim(floor / 3, 10)
    ax.axhline(1e-3, color="0.4", lw=0.8, ls="--")
    ax.text(len(labels) - 0.5, 1.4e-3, "a typical summary tolerance, 1e-3", ha="right", fontsize=8)
    ax.axvline(len(SUMMARY) - 0.5, color="0.6", lw=0.8)
    ax.set_xticks(np.arange(len(labels)))
    ax.set_xticklabels([lbl.replace("_", "\n") for lbl in labels], fontsize=8)
    ax.set_ylabel("relative difference from the correct run")
    ax.set_title("left: the summary numbers a solver reports;  right: the fields, point by point\n"
                 f"(a difference of exactly zero is drawn at {floor:g})", fontsize=9)
    ax.legend(fontsize=8)
    fig.savefig(outdir / "c4-six-languages-fig2.png", dpi=90)
    plt.close(fig)

    print(f"chapter 1's billows, {NX} x {NY}, t = {STEPS / 1000:.1f}; relative difference from the correct run")
    print(f"{'':16}" + "".join(f"{n:>16}" for n in wrong))
    for k in labels:
        print(f"{k:16}" + "".join(f"{rel(good[k], res[k]):16.1e}" for res in wrong.values()))
    flip = wrong["conjugate FFT"]
    print("conjugate FFT: u, v exactly reversed:",
          bool(np.array_equal(flip["u_velocity"], -good["u_velocity"])
               and np.array_equal(flip["v_velocity"], -good["v_velocity"])))
    print(f"wrote {outdir}/c4-six-languages-fig{{1,2}}.png")


if __name__ == "__main__":
    main()
