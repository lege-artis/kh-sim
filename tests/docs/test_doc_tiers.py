"""test_doc_tiers.py -- the engineer and shad documentation tiers (wave 4).

What a reader of docs/ relies on, checked:
  * the chapters exist and carry their license line;
  * every relative link and image resolves inside the kh-sim tree, and none
    points into the workspace folders that the public mirror leaves out;
  * figures are present and at most 300 KB;
  * no forbidden identifier of the sanitisation policy appears in docs/ or
    examples/; the scripts are ASCII-only;
  * the code a chapter shows runs, and prints what the chapter says it prints
    (the Hello-World and Taylor-Green snippets, the figure script of engineer 01,
    the example script of shad 01).

The two script runs take about a minute each:  python -m pytest tests/docs -q
Needs NumPy and Matplotlib.

ASCII-only source.
"""

import re
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]          # .../kh-sim
DOCS = ROOT / "docs"
EXAMPLES = ROOT / "examples"

CHAPTERS = [
    DOCS / "engineer" / "en" / "00-quick-start.md",
    DOCS / "engineer" / "en" / "01-what-kh-sim-computes.md",
    DOCS / "engineer" / "en" / "02-numerical-methods.md",
    DOCS / "shad" / "00-prologue.md",
    DOCS / "shad" / "01-the-billows.md",
]
MARKDOWN = CHAPTERS + [DOCS / "shad" / "README.md"]
ENGINEER_FIGURE_SCRIPT = DOCS / "engineer" / "en" / "figures" / "render_kh_evolution_3frame.py"
SHAD_EXAMPLE = EXAMPLES / "shad" / "c1-shear-layers" / "main.py"

MAX_FIGURE_BYTES = 300 * 1024
# The sanitisation policy's forbidden identifiers (public-mirror context).
FORBIDDEN = re.compile("SUPIN|Boura(?:\u010d|c)ka|MI-M-T|MIM2000|Yamyang|Improwave", re.I)
# Workspace folders that are not part of the public mirror.
PRIVATE_DIRS = ("_handoffs", "_audit", "_config", "_specs")
LINK = re.compile(r"!?\[[^\]]*\]\(([^)\s]+)\)")
TIMEOUT = 900


def _text(p):
    return p.read_text(encoding="utf-8")


def _links(md):
    """(target, is_image) for every inline markdown link outside code blocks."""
    out, in_code = [], False
    for line in _text(md).splitlines():
        if line.lstrip().startswith("```"):
            in_code = not in_code
            continue
        if in_code:
            continue
        for m in LINK.finditer(line):
            out.append((m.group(1), m.group(0).startswith("!")))
    return out


def _block_after(md, marker):
    """The fenced block that follows the HTML comment <!-- marker -->."""
    lines = _text(md).splitlines()
    i = lines.index(f"<!-- {marker} -->")
    assert lines[i + 1].startswith("```"), f"{md.name}: no fenced block after {marker}"
    j = i + 2
    while not lines[j].startswith("```"):
        j += 1
    return lines[i + 2:j]


def _run(args, stdin=None):
    r = subprocess.run([sys.executable] + args, input=stdin, capture_output=True, text=True,
                       cwd=ROOT, timeout=TIMEOUT)
    assert r.returncode == 0, f"exit {r.returncode}\n{r.stdout}\n{r.stderr}"
    return [ln.rstrip() for ln in r.stdout.splitlines()]


# -- structure ---------------------------------------------------------------------

def test_chapters_present_with_license():
    for md in MARKDOWN:
        assert md.is_file(), md
        head = "\n".join(_text(md).splitlines()[:12])
        assert "CC-BY-SA-4.0" in head or md.name == "README.md", f"{md.name}: no license line"
    assert "CC-BY-SA-4.0" in _text(DOCS / "shad" / "README.md")


def test_links_resolve_inside_the_public_tree():
    bad = []
    for md in MARKDOWN:
        for target, _ in _links(md):
            if re.match(r"[a-z]+:", target) or target.startswith("#"):
                continue
            path = target.split("#", 1)[0]
            dest = (md.parent / path).resolve()
            rel = dest.relative_to(ROOT) if dest.is_relative_to(ROOT) else None
            if rel is None:
                bad.append(f"{md.name}: {target} leaves the kh-sim tree")
            elif any(part in PRIVATE_DIRS for part in rel.parts):
                bad.append(f"{md.name}: {target} points into a workspace folder")
            elif not dest.exists():
                bad.append(f"{md.name}: {target} does not exist")
    assert not bad, "\n".join(bad)


def test_figures_present_and_small():
    images = [(md.parent / t).resolve() for md in MARKDOWN for t, img in _links(md) if img]
    assert len(images) == 4, images                    # engineer 01: 1, shad 01: 3
    for p in images:
        assert p.is_file(), p
        assert p.stat().st_size <= MAX_FIGURE_BYTES, f"{p.name}: {p.stat().st_size} bytes"


def test_no_forbidden_identifiers():
    hits = []
    for base in (DOCS, EXAMPLES):
        for p in base.rglob("*"):
            if p.is_file() and p.suffix in (".md", ".py", ".txt", ".tex"):
                for n, line in enumerate(_text(p).splitlines(), 1):
                    if FORBIDDEN.search(line):
                        hits.append(f"{p.relative_to(ROOT)}:{n}")
    assert not hits, hits


def test_scripts_are_ascii():
    for p in (ENGINEER_FIGURE_SCRIPT, SHAD_EXAMPLE, Path(__file__)):
        p.read_bytes().decode("ascii")


# -- the code the chapters show ------------------------------------------------------

@pytest.mark.parametrize("md, code, output", [
    (CHAPTERS[0], "hello-world", "hello-world-output"),
    (CHAPTERS[1], "taylor-green", "taylor-green-output"),
])
def test_snippet_prints_what_the_chapter_says(md, code, output):
    src = "\n".join(_block_after(md, code)) + "\n"
    assert _run(["-"], stdin=src) == [ln.rstrip() for ln in _block_after(md, output)]


@pytest.mark.slow
def test_engineer_figure_script(tmp_path):
    out = tmp_path / "kh-evolution-3frame.png"
    got = _run([str(ENGINEER_FIGURE_SCRIPT), "--out", str(out)])
    assert got[:3] == [
        "t = 0.0  kinetic_energy = 0.1208  enstrophy = 53.34  max_vorticity = 40.13",
        "t = 0.8  kinetic_energy = 0.1143  enstrophy = 44.90  max_vorticity = 35.25",
        "t = 1.6  kinetic_energy = 0.1091  enstrophy = 35.59  max_vorticity = 34.46",
    ]
    assert 0 < out.stat().st_size <= MAX_FIGURE_BYTES
    # The chapter quotes these numbers (01 s5).
    text = _text(CHAPTERS[1])
    for quoted in ("40.1 \u2192 35.2 \u2192 34.5", "53.3 \u2192 44.9 \u2192 35.6",
                   "0.121 \u2192 0.114 \u2192 0.109"):
        assert quoted in text, quoted


@pytest.mark.slow
def test_shad_example(tmp_path):
    got = _run([str(SHAD_EXAMPLE), "--outdir", str(tmp_path)])
    expected = [ln.rstrip() for ln in _block_after(CHAPTERS[4], "c1-output")]
    strip = lambda lines: [ln for ln in lines if not ln.startswith("wrote ")]  # noqa: E731
    assert strip(got) == strip(expected)
    for i in (1, 2, 3):
        png = tmp_path / f"c1-shear-layers-fig{i}.png"
        assert 0 < png.stat().st_size <= MAX_FIGURE_BYTES, png
