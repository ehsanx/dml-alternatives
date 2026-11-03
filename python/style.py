"""Shared style for the Python-drawn manuscript figures (Figures 3, 4, 9 and 10).

One encoding across the estimator figures:
  * hue    = nuisance-fitting scheme (full sample / single cross-fit / double cross-fit),
             grey for the singly-robust baselines (g-computation, IPTW);
  * marker = combination rule (circle = one-step: AIPW, DML, DC-AIPW;
             diamond = targeting: TMLE, CV-TMLE, DC-TMLE).
Colours are Okabe-Ito hues. The figures are drawn at the manuscript's text width (5.147 in) with
8-9 pt lettering, in Arial or, where Arial is not installed, a metric-compatible substitute
(Liberation Sans or Arimo): same text and layout, possibly not byte-identical output.
"""
from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("pdf")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib import font_manager as fm  # noqa: E402

# ---- colour -------------------------------------------------------------------------
C = dict(
    full="#D55E00",        # full-sample (in-sample) nuisance fits      Okabe-Ito vermillion
    cf="#0072B2",          # single cross-fit (out of fold)             Okabe-Ito blue
    dc="#009E73",          # double cross-fit (disjoint g / Q folds)    Okabe-Ito bluish green
    sr="#8A8A8A",          # singly-robust baselines (de-emphasis grey)
    treated="#AA3377",     # overlap plots: treated / quit (A = 1)
    control="#B8860B",     # overlap plots: comparison / did not quit (A = 0)
    ink="#1A1A1A", ink2="#444444", muted="#6E6E6E",
    grid="#E3E3E3", axis="#9A9A9A",
    navy="#002145",        # the median band of Figure 10 and the target box of Figure 3
    accent="#C8102E",      # sparing emphasis (one mark per figure at most)
)
SCHEME_COL = {"full": C["full"], "cf": C["cf"], "dc": C["dc"], "sr": C["sr"]}
RULE_MARK = {"onestep": "o", "targeting": "D", "plugin": "o"}

# ---- size: the manuscript's text width with 8-9 pt lettering ------------------------
P = dict(fs=8, fs_small=8, fs_label=9, fs_title=9, lw=1.1, lw_thin=0.6, ms=5.2, ms_small=2.6, width=5.147)
FONTS = ("Arial", "Liberation Sans", "Arimo")


def pick_font() -> str:
    """The first installed of Arial, Liberation Sans and Arimo."""
    for name in FONTS:
        try:
            fm.findfont(fm.FontProperties(family=name), fallback_to_default=False)
            return name
        except ValueError:
            continue
    raise SystemExit("install Arial, Liberation Sans or Arimo")


def setup() -> dict:
    """Set rcParams for the figures; the same font for text and mathtext."""
    fam = pick_font()
    P["family"] = fam
    plt.rcParams.update({
        "pdf.fonttype": 42, "ps.fonttype": 42,
        "font.family": fam, "font.size": P["fs"],
        "mathtext.fontset": "custom", "mathtext.rm": fam,
        "mathtext.it": f"{fam}:italic", "mathtext.bf": f"{fam}:bold",
        "axes.edgecolor": C["axis"], "axes.linewidth": P["lw_thin"] * 0.6,
        "axes.labelcolor": C["ink"], "axes.labelsize": P["fs_label"],
        "axes.titlesize": P["fs_title"], "axes.titleweight": "semibold",
        "axes.titlecolor": C["ink"], "axes.spines.top": False, "axes.spines.right": False,
        "xtick.color": C["ink2"], "ytick.color": C["ink2"],
        "xtick.labelsize": P["fs"], "ytick.labelsize": P["fs"],
        "xtick.major.width": P["lw_thin"] * 0.6, "ytick.major.width": P["lw_thin"] * 0.6,
        "xtick.major.size": P["fs"] * 0.25, "ytick.major.size": P["fs"] * 0.25,
        "legend.frameon": False, "legend.fontsize": P["fs"],
        "text.color": C["ink"], "figure.dpi": 72, "savefig.dpi": 72,
        "lines.solid_capstyle": "butt", "lines.dash_capstyle": "butt",
        "svg.fonttype": "none",
    })
    return P


def save(fig, path: Path) -> Path:
    """Write a vector PDF with embedded fonts and no creation date (so reruns compare equal)."""
    fig.savefig(path, metadata={"CreationDate": None})
    plt.close(fig)
    return path


def legend_handles(entries):
    """entries: list of (label, color, marker, filled) -> proxy artists for fig.legend."""
    from matplotlib.lines import Line2D
    hs = []
    for lab, col, mk, filled in entries:
        hs.append(Line2D([], [], ls="none", marker=mk, ms=P["ms"] * (0.85 if mk == "D" else 1.0),
                         mfc=col if filled else "white", mec=col, mew=P["lw_thin"] * 1.2, label=lab))
    return hs
