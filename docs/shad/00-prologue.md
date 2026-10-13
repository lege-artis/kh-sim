# Just Shad's Guide to the Kelvin–Helmholtz Galaxy

> *placeholder for the kh-sim icon — one great curling billow, mid-roll, looking as though it knows something about the weather that you do not*

*License: CC-BY-SA-4.0 (text); the example scripts are Apache-2.0.*

## What this is

Most books about fluid instability begin by writing the Navier–Stokes equations on the page. They write them in full vector form, with the pressure gradient and the viscous Laplacian and the convective term that is nonlinear in a way the authors clearly find thrilling and the reader clearly finds alarming, and they apologise for none of it. The reader, having opened the book hoping to understand why clouds sometimes look like a row of breaking waves, closes it again.

This guide does not do that. This guide starts with the clouds.

You have almost certainly seen the Kelvin–Helmholtz instability. It is the row of curling wave-crests that sometimes appears along the top of a cloud layer, each one rolled over in the same direction like a line of frozen surf. It is the swirl where milk meets coffee before you stir. It is the ripple that runs along a flag, and the roll-up at the edge of a jet of smoke. Two things are sliding past each other — fast air over slow air, cream over coffee — and the sliding, given the smallest excuse, curls the boundary into vortices. **Don't panic** if "vorticity" and "shear layer" are not yet words you own; by the end of this guide you will have watched the billows roll up in a simulation, measured how fast they grow, and seen why the same physics had to be written in six programming languages before anyone was allowed to call it correct.

What you need to bring:

- a curiosity about why moving fluids curl
- the willingness to read a plot of a swirling field
- a Python interpreter with NumPy and Matplotlib, if you want to run the code (optional, but it is more fun)
- about 20 minutes per chapter
- no fluid-dynamics course, no remembered vector calculus

The equations will show up, briefly, and we will deal with them when they do — and always with the rigorous version one link away.

## What this is *not*

This is not the canonical reference. Anywhere this guide says "the flow does X", there is a more careful statement, with the derivation attached, in [`../../shared/canonical-equations/`](../../shared/canonical-equations/). Wanting to know *why* rather than *that* is exactly the signal to follow that link.

It is also not the engineer's manual. If you want to build and run the solver, or understand the pseudo-spectral method and the tolerance gates, that is the [engineer tier](../engineer/en/00-quick-start.md). This guide is the story; that one is the recipe.

And it is not about real weather, quite. kh-sim's fluid is two-dimensional, has one density throughout, and lives in a box that repeats itself in both directions. Real billow clouds have none of those luxuries. Chapter 1 is honest about what that costs.

## How to use this guide

Chapters are self-contained but written to be read in order; each adds one idea.

| Chapter | Status | The new thing |
|---|---|---|
| [1 — The billows](01-the-billows.md) | v0.1.0 | Where you have seen KH, why a smooth shear will not stay smooth, and a simulation that rolls up the billows and measures how fast they grow |
| 2 — The vortex and the stream | queued, v0.1.x | Vorticity and stream function — the two-quantity language that makes two-dimensional flow tractable without ever mentioning pressure |
| 3 — Reading the growth rate | queued, v0.1.x | The dispersion relation: how a 1964 theorem lets the code test itself, and what it takes to make that test fair |
| 4 — Why six languages | queued, v0.1.x | What "canonical reference" costs: element-wise agreement, and the sign error that matching energy totals never caught |
| 5 onward — Shear layers in the wild | queued, v0.1.x | Billow clouds and the Richardson number, mixing layers in the ocean, the edge of a jet — each its own chapter |

Each chapter has three figures — the input, the thing that happens, the takeaway — and a runnable script under [`../../examples/shad/`](../../examples/shad/) that renders them. Run a script and you get the same figures back, modulo your Matplotlib's fonts and its views on colour.

## Tone

Dry. Technical. Slightly amused. Occasionally wandering off after a related observation for a paragraph before returning to the point. The narrator is allowed to find things interesting on the reader's behalf, partly because fluid instability genuinely is interesting once you can see what it is doing, and partly because the narrator is the only one available to keep the reader company while they work it out.

If a sentence in here makes you think *"the author is being clever for the sake of being clever"*, please file an issue. The point is to make the instability feel obvious, not to perform cleverness at you.

## One idea to carry through

The whole point of kh-sim is a small act of stubbornness: refusing to trust a swirling picture just because it *looks* like the real thing. A simulation can produce beautiful billows for entirely wrong reasons. Everything in this guide, and everything in the tiers below it, exists to replace "it looks right" with "it hit the number theory said it would". That is the difference between a pretty animation and a reference.

## License

This guide is CC-BY-SA-4.0; the example scripts are Apache-2.0. The instability itself belongs to Hermann von Helmholtz (1868) and William Thomson, Lord Kelvin (1871), who described it a century and a half before anyone could simulate it, and who are in no position to object to the use made of their names.

---

**Next:** [Chapter 1 — The billows](01-the-billows.md)
