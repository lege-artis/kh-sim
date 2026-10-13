# Just Shad's Guide to the Kelvin–Helmholtz Galaxy

The story-first documentation tier of lege-artis/kh-sim: the physics told as a story, for readers who want to see why moving fluids curl and how a simulation earns the right to be believed, without first learning vector calculus. Every "the flow does X" is one link away from its rigorous statement in the canonical tier.

## Chapters

| # | Title | Status | Adds |
|---|---|---|---|
| 00 | [Prologue](00-prologue.md) | v0.1.0 | What this is and is not, how to read it |
| 01 | [The billows](01-the-billows.md) | v0.1.0 | Where you have seen KH; the pressure feedback that grows a ripple; a simulation that rolls up the billows, measures their growth rate, and explains why it is three quarters of the textbook's |
| 02 | The vortex and the stream | queued, v0.1.x | Vorticity and stream function, without pressure |
| 03 | Reading the growth rate | queued, v0.1.x | The dispersion relation as the code's self-test |
| 04 | Why six languages | queued, v0.1.x | Element-wise agreement and the sign error that energy totals never caught |
| 05+ | Shear layers in the wild | queued, v0.1.x | Clouds and the Richardson number, ocean mixing layers, jets |

## Example scripts

The figures of each chapter are rendered by a script in [`../../examples/shad/`](../../examples/shad/):

- [`c1-shear-layers/main.py`](../../examples/shad/c1-shear-layers/main.py) → `figures/c1-shear-layers-fig{1,2,3}.png`

Each script runs on its own with `python main.py` (Python 3.10+, NumPy, Matplotlib; no network) and computes everything with the kh-sim Python reference, [`shared/physics/kh_physics.py`](../../shared/physics/kh_physics.py).

## Iconography

The tier's icon is not drawn yet; the prologue carries a placeholder.

## License

This tier is **CC-BY-SA-4.0**, like the rest of `docs/`. The example scripts in `examples/shad/` are **Apache-2.0**, like the code.

## Cross-references

- Engineer tier (how to run it, what it computes, how): [`../engineer/en/`](../engineer/en/00-quick-start.md)
- Canonical tier (the equations): [`../../shared/canonical-equations/`](../../shared/canonical-equations/)
