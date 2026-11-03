"""Checks that tools/compare_pdf.py fails on real changes and passes on harmless ones.

    python tools/test_compare_pdf.py

Draws a small PDF with matplotlib (DejaVu Sans, so no Arial is needed), then compares variants with it:
  must pass: the same drawing again; the same file with its objects renumbered and reordered;
  must fail: a changed crop box; a label moved by 1 pt; a missing label; the regular and bold fonts swapped in
  the resource dictionary; two overlapping shapes drawn in the other order.
The "sorted streams" column shows what a weaker rule (all decoded streams, sorted) would say: it passes the
changed crop box and the swapped fonts. Exit status 1 if any case does not behave as expected.
"""
from __future__ import annotations

import shutil
import sys
import tempfile
from pathlib import Path

import matplotlib

matplotlib.use("pdf")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.patches import Rectangle  # noqa: E402
from pypdf import PdfReader, PdfWriter  # noqa: E402
from pypdf.generic import NameObject, RectangleObject  # noqa: E402

import compare_pdf  # noqa: E402

plt.rcParams.update({"pdf.fonttype": 42, "font.family": "DejaVu Sans"})


def draw(path: Path, dx_pt=0.0, beta=True, red_on_top=False):
    fig = plt.figure(figsize=(3, 2))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, 3)
    ax.set_ylim(0, 2)
    ax.axis("off")
    ax.add_patch(Rectangle((0.3, 0.3), 1.2, 0.8, color="#D55E00", zorder=3 if red_on_top else 1))
    ax.add_patch(Rectangle((0.9, 0.6), 1.2, 0.8, color="#0072B2", zorder=2))
    ax.text(0.6 + dx_pt / 72, 1.6, "alpha", fontsize=10)
    if beta:
        ax.text(2.0, 0.3, "beta", fontsize=10, fontweight="bold")
    fig.savefig(path, metadata={"CreationDate": None})
    plt.close(fig)


def rewrite(src: Path, dst: Path, change=None):
    """Copy a PDF through pypdf, which renumbers and reorders its objects; optionally change the page."""
    w = PdfWriter()
    page = w.add_page(PdfReader(str(src)).pages[0])
    if change:
        change(page)
    with open(dst, "wb") as f:
        w.write(f)


def crop(page):
    b = page.mediabox
    page.cropbox = RectangleObject([float(b.left) + 6, float(b.bottom) + 6, float(b.right) - 6, float(b.top) - 6])


def swap_fonts(page):
    fonts = page["/Resources"]["/Font"]
    a, b = sorted(fonts.keys())[:2]
    fonts[NameObject(a)], fonts[NameObject(b)] = fonts[b], fonts[a]


def font_numbers(path: Path) -> list[int]:
    fonts = PdfReader(str(path)).pages[0]["/Resources"]["/Font"]
    return [fonts.raw_get(k).idnum for k in sorted(fonts.keys())]


def sorted_streams_equal(a: Path, b: Path) -> bool:
    """A weaker rule, for comparison: every decoded stream, sorted, ignoring everything outside streams."""
    import re
    import zlib

    def streams(p):
        out = []
        for m in re.finditer(rb"<<(.*?)>>\s*stream\r?\n(.*?)\r?\nendstream", p.read_bytes(), re.S):
            if b"/XRef" in m.group(1):
                continue
            try:
                out.append(zlib.decompress(m.group(2)))
            except zlib.error:
                out.append(m.group(2))
        return sorted(out)
    return streams(a) == streams(b)


def main() -> int:
    tmp = Path(tempfile.mkdtemp())
    try:
        base = tmp / "base.pdf"
        draw(base)
        cases = []
        draw(tmp / "same.pdf")
        cases.append(("same drawing again", tmp / "same.pdf", True))
        rewrite(base, tmp / "reordered.pdf")
        if font_numbers(tmp / "reordered.pdf") == font_numbers(base):
            raise SystemExit("the rewrite did not renumber the font objects; the reordering case proves nothing")
        cases.append(("objects renumbered and reordered", tmp / "reordered.pdf", True))
        rewrite(base, tmp / "cropped.pdf", crop)
        cases.append(("crop box 6 pt smaller on each side", tmp / "cropped.pdf", False))
        draw(tmp / "moved.pdf", dx_pt=1.0)
        cases.append(("label moved by 1 pt", tmp / "moved.pdf", False))
        draw(tmp / "missing.pdf", beta=False)
        cases.append(("label missing", tmp / "missing.pdf", False))
        rewrite(base, tmp / "swapped.pdf", swap_fonts)
        cases.append(("regular and bold fonts swapped", tmp / "swapped.pdf", False))
        draw(tmp / "order.pdf", red_on_top=True)
        cases.append(("overlapping shapes in the other order", tmp / "order.pdf", False))

        bad = 0
        print(f"{'case':<40}{'expected':<10}{'got':<10}{'sorted streams':<15}reason")
        for name, path, want in cases:
            got, why = compare_pdf.compare_full(path, base, dpi=144, do_render=True)
            ok = got == want
            bad += not ok
            prev = "pass" if sorted_streams_equal(path, base) else "fail"
            print(f"{name:<40}{'pass' if want else 'fail':<10}{'pass' if got else 'fail':<10}{prev:<15}{why}"
                  + ("" if ok else "   <-- UNEXPECTED"))
        print(f"\n{len(cases) - bad} of {len(cases)} cases behave as expected")
        return 1 if bad else 0
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


if __name__ == "__main__":
    sys.exit(main())
