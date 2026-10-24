#!/usr/bin/env python3
"""
nr1986_json_to_namelist.py  --  CP-7F-E Fortran loader (JSON -> namelist)

The Fortran backend has no JSON parser; its IO idiom is Fortran namelist
(see backends/fortran/src/kh_io.f90). Python (oracle) and Pascal (fpjson) read
the golden vector JSON directly; only Fortran needs this transcode.

CANONICAL SCHEMA: nr1986-pascal-ch19/v1, as defined verbatim in
  development record, NR 1986 example lookup.
That doc is the single source of truth for what Pete transcribes. This tool
parses THAT schema (problem.*, expected.sample_points[], expected.integrals[]).

Input : shared/golden-vectors/nr1986-pascal-ch19.json   (schema nr1986-pascal-ch19/v1)
Output: backends/fortran/tests/nr1986_pascal_ch19.nml    (namelist /NR1986_CP7FE/)

Behaviour:
  * Refuses to emit if any required value is still an unfilled placeholder
    (contains '<' ... '>' or 'FILL') -> exit 2, listing the offenders. The
    Fortran test then finds no .nml and SKIPs -- never a false green.
  * Emits sample points as parallel arrays (sample_i, sample_j, sample_val) plus
    n_samples, and the first integral check if present.
  * On success exit 0.

Usage:
  python3 tools/nr1986_json_to_namelist.py \
      [shared/golden-vectors/nr1986-pascal-ch19.json] \
      [-o backends/fortran/tests/nr1986_pascal_ch19.nml]
"""
import argparse
import json
import sys
from pathlib import Path

DEFAULT_IN = "shared/golden-vectors/nr1986-pascal-ch19.json"
DEFAULT_OUT = "backends/fortran/tests/nr1986_pascal_ch19.nml"
MAX_SAMPLES = 64  # fixed Fortran array bound; bump in both places if exceeded


def _is_placeholder(v):
    """A still-unfilled template value: angle-bracket stub or FILL marker."""
    if isinstance(v, str):
        s = v.strip()
        return ("<" in s and ">" in s) or s.upper().startswith("FILL") or s == ""
    return False


def _f_str(s):
    """Fortran namelist string literal: single-quoted, quotes doubled, ASCII."""
    s = "" if s is None else str(s)
    s = s.encode("ascii", "replace").decode("ascii")
    return "'" + s.replace("'", "''") + "'"


def main():
    ap = argparse.ArgumentParser(description="CP-7F-E JSON -> Fortran namelist")
    ap.add_argument("infile", nargs="?", default=DEFAULT_IN)
    ap.add_argument("-o", "--outfile", default=DEFAULT_OUT)
    args = ap.parse_args()

    inp = Path(args.infile)
    if not inp.exists():
        print(f"[nr1986] no input JSON at {inp} -- nothing to do "
              f"(Fortran test will SKIP). Expected pre-transcription.")
        return 0

    data = json.loads(inp.read_text())
    if data.get("schema") not in (None, "nr1986-pascal-ch19/v1"):
        print(f"[nr1986] WARNING: unexpected schema '{data.get('schema')}' "
              f"(expected nr1986-pascal-ch19/v1). Parsing best-effort.",
              file=sys.stderr)

    src = data.get("source", {})
    prob = data.get("problem", {})
    exp = data.get("expected", {})
    samples = exp.get("sample_points", []) or []
    integrals = exp.get("integrals", []) or []

    # --- completeness gate -----------------------------------------------------
    bad = []
    for label, obj, keys in (
        ("source", src, ("chapter", "section")),
        ("problem", prob, ("type", "grid_nx", "grid_ny", "boundary")),
    ):
        for k in keys:
            if k not in obj or _is_placeholder(obj.get(k)):
                bad.append(f"{label}.{k}")
    if not samples:
        bad.append("expected.sample_points (empty)")
    for n, sp in enumerate(samples):
        for k in ("i", "j", "value"):
            if k not in sp or _is_placeholder(sp.get(k)):
                bad.append(f"expected.sample_points[{n}].{k}")
    if bad:
        print("[nr1986] transcription INCOMPLETE -- unfilled/placeholder fields:",
              file=sys.stderr)
        for f in bad:
            print(f"          - {f}", file=sys.stderr)
        print("[nr1986] no namelist emitted; Fortran CP-7F-E will SKIP.",
              file=sys.stderr)
        return 2

    if len(samples) > MAX_SAMPLES:
        print(f"[nr1986] {len(samples)} sample points > MAX_SAMPLES={MAX_SAMPLES}; "
              f"bump the bound in this tool AND run_cp7fe.", file=sys.stderr)
        return 4

    # --- coerce ----------------------------------------------------------------
    try:
        nx = int(prob["grid_nx"]); ny = int(prob["grid_ny"])
        si = [int(s["i"]) for s in samples]
        sj = [int(s["j"]) for s in samples]
        sv = [float(s["value"]) for s in samples]
        tol_rel = float(exp.get("tolerance_rel", 1.0e-3))
        integ_name = str(integrals[0]["name"]) if integrals else ""
        integ_val = float(integrals[0]["value"]) if integrals else 0.0
        params = prob.get("params", {}) or {}
        src_i = int(params.get("source_i_1based", 0))
        src_j = int(params.get("source_j_1based", 0))
        src_val = float(params.get("source_value", 0.0))
        m_resp = int(params.get("m_response", 0))
    except (KeyError, TypeError, ValueError) as e:
        print(f"[nr1986] could not coerce a required numeric field: {e}",
              file=sys.stderr)
        return 3

    def arr(name, vals, fmt):
        body = ", ".join(fmt(v) for v in vals)
        return f"  {name:<14}= {body}"

    out = Path(args.outfile)
    out.parent.mkdir(parents=True, exist_ok=True)
    L = []
    L.append("! Auto-generated by tools/nr1986_json_to_namelist.py")
    L.append(f"! Source JSON: {inp.as_posix()}  (schema nr1986-pascal-ch19/v1)")
    L.append("! DO NOT EDIT BY HAND -- regenerate from the golden vector.")
    L.append("&NR1986_CP7FE")
    L.append("  data_present  = .true.")
    L.append(f"  chapter       = {_f_str(src.get('chapter'))}")
    L.append(f"  section       = {_f_str(src.get('section'))}")
    L.append(f"  example_name  = {_f_str(src.get('example_name'))}")
    L.append(f"  prob_type     = {_f_str(prob.get('type'))}")
    L.append(f"  boundary      = {_f_str(prob.get('boundary'))}")
    L.append(f"  grid_nx       = {nx}")
    L.append(f"  grid_ny       = {ny}")
    L.append(f"  n_samples     = {len(samples)}")
    L.append(arr("sample_i", si, lambda v: str(v)))
    L.append(arr("sample_j", sj, lambda v: str(v)))
    L.append(arr("sample_val", sv, lambda v: repr(v)))
    L.append(f"  integ_name    = {_f_str(integ_name)}")
    L.append(f"  integ_val     = {integ_val!r}")
    L.append(f"  tolerance_rel = {tol_rel!r}")
    # problem.params (2026-10-11, CP-7F-E asserts): the E2 source and the E1
    # response length; 0 when the golden vector has none.
    L.append(f"  src_i         = {src_i}")
    L.append(f"  src_j         = {src_j}")
    L.append(f"  src_val       = {src_val!r}")
    L.append(f"  m_resp        = {m_resp}")
    L.append("/")
    out.write_text("\n".join(L) + "\n")
    print(f"[nr1986] wrote {out}  (n_samples={len(samples)}, "
          f"tol_rel={tol_rel!r}, type={prob.get('type')!r})")
    return 0


if __name__ == "__main__":
    sys.exit(main())
