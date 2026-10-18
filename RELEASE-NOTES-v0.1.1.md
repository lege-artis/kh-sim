# lege-artis/kh-sim — v0.1.1

**Tag:** `v0.1.1`
**Type:** CI fix. No change to the solver, the tests or the documentation chapters.

## What changed

- **Scala workflow:** `sbt` is not preinstalled on the Ubuntu 24.04 GitHub runner image, so the v0.1.0 Scala job stopped at *Compile* with exit code 127. The workflow now installs it with `sbt/setup-sbt@v1` after setting up Java 17.
- **All six workflows** run on `ubuntu-24.04` instead of `ubuntu-latest`, which moves to Ubuntu 26 from 2026-10-19; the toolchains the jobs install (gfortran, FPC, CMake, Python 3.13) stay the ones they were tested with.
- `README.md` status line.

## CI on v0.1.0, for the record

Fortran, Rust, C++, Pascal and Python/NumPy green on the first run; Scala failed as above. With v0.1.1 all six are expected green.

See `RELEASE-NOTES-v0.1.0.md` for the content and validation of the release.
