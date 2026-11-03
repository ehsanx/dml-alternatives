"""Compare figure PDFs with the manuscript's reference copies.

    python tools/compare_pdf.py output/figures figures/reference              # full check
    python tools/compare_pdf.py --dpi 200 output/figures figures/reference    # render at another resolution
    python tools/compare_pdf.py a.pdf b.pdf                                   # one pair
    python tools/compare_pdf.py --text-only output/figures figures/reference  # limited fallback

Each PDF in the reference folder is compared with the file of the same name in the first folder. The full
check passes a pair only if all three parts agree:

1. Geometry: the page count, and each page's MediaBox, CropBox, BleedBox, TrimBox, ArtBox, /Rotate and
   /UserUnit.
2. Structure: each page dictionary, with every indirect reference replaced by the object it points to and every
   stream by its decoded data. The content streams are compared byte for byte, so drawing order counts, and each
   resource name (/F1, /GS0, ...) must point to the same font program, widths, encoding and graphics state. Object
   numbers and the order of objects in the file are ignored: matplotlib writes embedded fonts sorted by their file
   path, so a copy installed elsewhere numbers the same objects differently. The document information dictionary
   (creation date, producer) is not part of a page and is ignored.
3. Rendering: every page rendered with the same renderer (PDFium, through pypdfium2) at the same resolution
   (144 dpi unless --dpi says otherwise), compared pixel by pixel.

Expect a full pass only with Arial and the package versions in renv.lock and requirements.txt.

--text-only is a limited fallback for machines without Arial, where the full check cannot pass: it compares the
extracted words only, so it misses changes to drawing, positions, sizes, colours and fonts. It uses Poppler's
pdftotext if it is on the PATH, or else pypdf, and prints which extractor ran.
--no-render skips part 3 when pypdfium2 is unavailable, and says so in the result.

Exit status 1 if any figure differs or is missing; 2 for a usage error.
"""
from __future__ import annotations

import argparse
import shutil
import subprocess
import sys
from pathlib import Path

from pypdf import PdfReader
from pypdf.generic import (ArrayObject, BooleanObject, DictionaryObject, FloatObject, IndirectObject,
                           NameObject, NullObject, NumberObject, StreamObject)

BOXES = ("/MediaBox", "/CropBox", "/BleedBox", "/TrimBox", "/ArtBox")
STREAM_CODING = {"/Length", "/Filter", "/DecodeParms", "/DL"}     # how a stream is stored, not what it holds


# ---- 1. geometry ---------------------------------------------------------------------------------
def geometry(reader: PdfReader) -> list[dict]:
    pages = []
    for page in reader.pages:
        g = {box: [round(float(v), 4) for v in page[box]] if box in page else None for box in BOXES[1:]}
        g["/MediaBox"] = [round(float(v), 4) for v in page.mediabox]
        g["/CropBox (effective)"] = [round(float(v), 4) for v in page.cropbox]
        g["/Rotate"] = int(page.get("/Rotate", 0))
        g["/UserUnit"] = float(page.get("/UserUnit", 1))
        pages.append(g)
    return pages


# ---- 2. structure ----------------------------------------------------------------------------------
class Canon:
    """Canonical form of a PDF object graph: references resolved, streams decoded, dictionaries sorted."""

    def __init__(self):
        self.memo, self.stack = {}, set()

    def __call__(self, obj):
        if isinstance(obj, IndirectObject):
            key = (obj.idnum, obj.generation)
            if key in self.memo:
                return self.memo[key]
            if key in self.stack:                       # a cycle: compare it by position, not by content
                return ("cycle",)
            self.stack.add(key)
            out = self(obj.get_object())
            self.stack.discard(key)
            self.memo[key] = out
            return out
        if isinstance(obj, StreamObject):
            head = tuple((k, self(v)) for k, v in sorted(obj.items()) if k not in STREAM_CODING)
            return ("stream", head, obj.get_data())
        if isinstance(obj, DictionaryObject):
            return ("dict", tuple((k, self(v)) for k, v in sorted(obj.items()) if k != "/Parent"))
        if isinstance(obj, ArrayObject):
            return ("array", tuple(self(v) for v in obj))
        if isinstance(obj, NameObject):
            return ("name", str(obj))
        if isinstance(obj, BooleanObject):
            return ("bool", bool(obj))
        if isinstance(obj, (NumberObject, FloatObject)):
            return ("num", round(float(obj), 6))
        if isinstance(obj, NullObject) or obj is None:
            return ("null",)
        if isinstance(obj, (bytes, str)):
            return ("str", obj.encode("latin-1", "replace") if isinstance(obj, str) else bytes(obj))
        return ("other", type(obj).__name__, repr(obj))


def structure(reader: PdfReader) -> list:
    canon = Canon()
    return [canon(page.indirect_reference or page) for page in reader.pages]


def first_difference(a, b, path="page") -> str:
    """Where two canonical forms first differ, for the report."""
    if type(a) is not type(b) or not isinstance(a, tuple):
        return path
    if a[:1] != b[:1]:
        return f"{path}: {a[0]} vs {b[0]}"
    kind = a[0]
    if kind == "dict":
        ka, kb = [k for k, _ in a[1]], [k for k, _ in b[1]]
        if ka != kb:
            return f"{path}: keys {sorted(set(ka) ^ set(kb))}"
        for (k, va), (_, vb) in zip(a[1], b[1]):
            if va != vb:
                return first_difference(va, vb, f"{path}{k}")
    if kind == "array":
        if len(a[1]) != len(b[1]):
            return f"{path}: array length {len(a[1])} vs {len(b[1])}"
        for i, (va, vb) in enumerate(zip(a[1], b[1])):
            if va != vb:
                return first_difference(va, vb, f"{path}[{i}]")
    if kind == "stream":
        if a[1] != b[1]:
            return first_difference(("dict", a[1]), ("dict", b[1]), path)
        return f"{path}: stream data ({len(a[2])} vs {len(b[2])} bytes)"
    return f"{path}: {a!r:.60} vs {b!r:.60}"


# ---- 3. rendering ----------------------------------------------------------------------------------
def render(path: Path, dpi: int) -> list:
    import numpy as np
    import pypdfium2 as pdfium
    doc = pdfium.PdfDocument(str(path))
    try:
        return [np.asarray(doc[i].render(scale=dpi / 72).to_pil().convert("RGB")) for i in range(len(doc))]
    finally:
        doc.close()


def compare_full(new: Path, ref: Path, dpi: int, do_render: bool) -> tuple[bool, str]:
    ra, rb = PdfReader(str(new)), PdfReader(str(ref))
    ga, gb = geometry(ra), geometry(rb)
    if len(ga) != len(gb):
        return False, f"page count {len(ga)} vs {len(gb)}"
    for i, (pa, pb) in enumerate(zip(ga, gb), 1):
        for k in pa:
            if pa[k] != pb[k]:
                return False, f"page {i} {k} {pa[k]} vs {pb[k]}"
    sa, sb = structure(ra), structure(rb)
    for i, (pa, pb) in enumerate(zip(sa, sb), 1):
        if pa != pb:
            return False, f"page {i} structure: {first_difference(pa, pb)}"
    note = f"{len(ga)} page(s), boxes, content and resources equal"
    if not do_render:
        return True, note + "; NOT rendered (--no-render)"
    import numpy as np
    for i, (pa, pb) in enumerate(zip(render(new, dpi), render(ref, dpi)), 1):
        if pa.shape != pb.shape:
            return False, f"page {i} rendered size {pa.shape[1]}x{pa.shape[0]} vs {pb.shape[1]}x{pb.shape[0]} px"
        n = int((np.abs(pa.astype(int) - pb.astype(int)).sum(axis=2) > 0).sum())
        if n:
            return False, f"page {i} rendered at {dpi} dpi: {n} pixels differ"
    return True, note + f"; rendered at {dpi} dpi: 0 pixels differ"


# ---- text-only fallback ----------------------------------------------------------------------------
def extractor():
    exe = shutil.which("pdftotext")
    if exe:
        ver = subprocess.run([exe, "-v"], capture_output=True, text=True)
        name = ((ver.stderr or ver.stdout).splitlines() or ["pdftotext"])[0]
        return name, lambda p: subprocess.run([exe, "-enc", "UTF-8", str(p), "-"], capture_output=True,
                                              check=True).stdout.decode("utf-8")
    return "pypdf", lambda p: "\n".join(pg.extract_text() for pg in PdfReader(str(p)).pages)


def main(argv: list[str]) -> int:
    ap = argparse.ArgumentParser(description="Compare figure PDFs with reference copies.", add_help=True)
    ap.add_argument("new", help="folder of drawn PDFs, or one PDF")
    ap.add_argument("ref", help="folder of reference PDFs, or one PDF")
    ap.add_argument("--dpi", type=int, default=144, help="rendering resolution (default 144)")
    ap.add_argument("--text-only", action="store_true", help="limited fallback: compare extracted words only")
    ap.add_argument("--no-render", action="store_true", help="skip the rendered comparison")
    a = ap.parse_args(argv)
    new, ref = Path(a.new), Path(a.ref)
    if ref.is_file():
        pairs = [(new, ref)]
    else:
        refs = sorted(ref.glob("*.pdf"))
        if not refs:
            print(f"no PDFs in {ref}")
            return 2
        pairs = [(new / r.name, r) for r in refs]
    do_render = not a.no_render
    if do_render and not a.text_only:
        try:
            import pypdfium2  # noqa: F401
        except ImportError:
            print("the rendered comparison needs pypdfium2 (requirements.txt); --no-render skips it")
            return 2
    if a.text_only:
        name, get = extractor()
        print(f"TEXT-ONLY FALLBACK (limited: extracted words only; drawing, positions, sizes, colours and fonts"
              f" are not compared). Extractor: {name}")
    bad = 0
    for n, r in pairs:
        if not n.exists():
            print(f"MISSING     {r.name}")
            bad += 1
            continue
        if a.text_only:
            same, why = get(n).split() == get(r).split(), "extracted words"
            label = "TEXT-SAME" if same else "TEXT-DIFFERS"
        else:
            same, why = compare_full(n, r, a.dpi, do_render)
            label = "IDENTICAL" if same else "DIFFERS"
        print(f"{label:<12}{r.name}: {why}")
        bad += not same
    kind = "identical" if not a.text_only else "with the same words (text-only fallback)"
    print(f"{len(pairs) - bad} of {len(pairs)} {kind}")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
