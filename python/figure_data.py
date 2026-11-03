"""Data readers and derived numbers for the Python-drawn figures (Figures 9 and 10).

Reads the frozen CSVs in docs/data/, or in the folder named by the DML_OUT environment variable:
  nhefs_panel_table.csv       NHEFS estimates (Figure 9, lower left)
  results_all3.csv            LaLonde/NSW-PSID estimates at one split (Figure 9, lower right)
  nhefs_overlap_data.csv      the NHEFS overlap diagnostic's scores (Figure 9, upper left)
  positivity_overlap.csv      the seeded LaLonde/NSW-PSID overlap diagnostic (Figure 9, upper right)
  nhefs_software.csv          single runs of five implementations (Figure 10)
  nhefs_matched_repeated.csv  20 seeded reruns of three of them (Figure 10)
"""
from __future__ import annotations

import os
from pathlib import Path

import numpy as np
import pandas as pd


def repo_root() -> Path:
    """The nearest folder above this file that holds R/estimators.R."""
    p = Path(__file__).resolve().parent
    while not (p / "R" / "estimators.R").exists() and p.parent != p:
        p = p.parent
    return p


DATA = Path(os.environ.get("DML_OUT", repo_root() / "docs" / "data"))

NH_SPLIT = ["DML (cross-fit AIPW)", "CV-TMLE", "DC-AIPW", "DC-TMLE"]
# display name, scheme, rule  (rows of every estimator panel, in this order)
ROWS = [
    ("g-computation", "sr", "plugin"), ("IPTW", "sr", "plugin"),
    ("AIPW", "full", "onestep"), ("TMLE", "full", "targeting"),
    ("DML", "cf", "onestep"), ("CV-TMLE", "cf", "targeting"),
    ("DC-AIPW", "dc", "onestep"), ("DC-TMLE", "dc", "targeting"),
]


def _one(d: pd.DataFrame, what: str) -> pd.Series:
    if len(d) != 1:
        raise ValueError(f"expected one row for {what}, got {len(d)}")
    return d.iloc[0]


def nhefs_panel() -> pd.DataFrame:
    """NHEFS rows: full-sample fits, and the split-based routes at r = 1 and r = 100."""
    p = pd.read_csv(DATA / "nhefs_panel_table.csv")
    full = p[p.r.isna()]
    out = []
    for disp, key in [("g-computation", "g-computation"), ("IPTW", "IPTW"),
                      ("AIPW", "AIPW (one-step)"), ("TMLE", "TMLE")]:
        z = _one(full[full.estimator == key], key)
        out.append(dict(name=disp, est=z.estimate, lo=z.ci_lo, hi=z.ci_hi, r=np.nan))
    for disp, key in zip(["DML", "CV-TMLE", "DC-AIPW", "DC-TMLE"], NH_SPLIT):
        for rr in (1, 100):
            z = _one(p[(p.estimator == key) & (p.r == rr)], f"{key} r={rr}")
            out.append(dict(name=disp, est=z.estimate, lo=z.ci_lo, hi=z.ci_hi, r=rr,
                            partition_sd=z.partition_sd))
    return pd.DataFrame(out)


def lalonde() -> pd.DataFrame:
    """LaLonde/NSW-PSID rows of results_all3.csv, with the display names of ROWS."""
    a = pd.read_csv(DATA / "results_all3.csv")
    a = a[a.dataset.str.startswith("LaLonde")].copy()
    a["name"] = a.estimator.map({
        "g-computation": "g-computation", "IPTW": "IPTW", "AIPW": "AIPW", "TMLE": "TMLE",
        "AIPW = DML": "DML", "TMLE (CV-TMLE)": "CV-TMLE",
        "AIPW (DC-AIPW)": "DC-AIPW", "TMLE (DC-TMLE)": "DC-TMLE"})
    return a


def software() -> dict:
    return dict(table=pd.read_csv(DATA / "nhefs_software.csv"),
                seeds=pd.read_csv(DATA / "nhefs_matched_repeated.csv"))


def overlap() -> dict:
    """NHEFS: the overlap diagnostic's scores (not seeded; also Figure 2). LaLonde: the seeded diagnostic."""
    return dict(nhefs=pd.read_csv(DATA / "nhefs_overlap_data.csv"),
                diag=pd.read_csv(DATA / "positivity_overlap.csv"))


def r_density(x, n=2 ** 14, cut=3):
    """Gaussian kernel density with R's default bandwidth (bw.nrd0), evaluated exactly (not binned) on a fine grid
    from min-3bw to max+3bw. R's density() bins the data and, at its default n = 512, understates narrow peaks:
    for the LaLonde comparison group (bw about 0.0014, below that grid's step) it gives about 123, against about 154
    here."""
    x = np.asarray(x, float)
    sd = x.std(ddof=1)
    iqr = np.subtract(*np.percentile(x, [75, 25]))
    lo = min(sd, iqr / 1.34) if iqr > 0 else sd
    bw = 0.9 * lo * len(x) ** (-0.2)
    grid = np.linspace(x.min() - cut * bw, x.max() + cut * bw, n)
    d = np.exp(-0.5 * ((grid[:, None] - x[None, :]) / bw) ** 2).sum(1) / (len(x) * bw * np.sqrt(2 * np.pi))
    return grid, d


def numbers() -> dict:
    """The numbers the figures print, unrounded. R/verify_numbers.R audits the manuscript's rounded values."""
    F = {}
    nh = nhefs_panel()
    full = nh[nh.r.isna()].set_index("name")
    r100 = nh[nh.r == 100].set_index("name")
    dr6 = np.r_[full.loc[["AIPW", "TMLE"], "est"], r100.est]
    F["nh_dr6_min"], F["nh_dr6_max"] = dr6.min(), dr6.max()
    s = software()
    v = s["table"][["AIPW", "TMLE", "DML_crossfit"]].to_numpy().ravel()
    v = v[~np.isnan(v)]
    F["sw_min"], F["sw_max"], F["sw_span"] = v.min(), v.max(), v.max() - v.min()
    med = s["seeds"][["hand_DML", "AIPW_pkg", "tmle_pkg"]].median()
    F["sw_med_gap"] = med.max() - med.min()
    F["sw_med"] = med.mean()
    lf = lalonde().set_index("name")
    F["la_gcomp"], F["la_iptw"] = lf.estimate["g-computation"], lf.estimate["IPTW"]
    dg = overlap()["diag"]
    lc = dg[dg.dataset.str.startswith("LaLonde") & (dg.group == "control")].pi
    F["la_ctrl_peak"] = r_density(lc)[1].max()          # clipped LaLonde comparison-group peak
    return F
