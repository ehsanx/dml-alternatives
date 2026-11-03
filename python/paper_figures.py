"""Figures 3, 4, 9 and 10 of the manuscript, drawn with matplotlib at the manuscript's text width.

    python python/paper_figures.py [output folder]        (default: output/figures)

writes concept_grid.pdf (Figure 3), fitting_schemes.pdf (Figure 4), overlap_estimates.pdf (Figure 9) and
software_seeds.pdf (Figure 10). Figures 3 and 4 are schematics; Figures 9 and 10 read the frozen CSVs listed in
figure_data.py. Encoding (style.py): hue = nuisance-fitting scheme, marker = combination rule (circle one-step,
diamond targeting), grey hollow = singly robust.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
from matplotlib.patches import Rectangle

import figure_data as fd
import style as st
from style import C, P, RULE_MARK, SCHEME_COL, plt

OUT = Path("output/figures")      # replaced by the command-line argument


# ---- helpers -------------------------------------------------------------------------
_RENDER = {}


def tw(s: str, size: float, weight: str = "normal", style: str = "normal") -> float:
    """Rendered width of s in inches, measured with an Agg renderer (as drawn)."""
    if "fig" not in _RENDER:
        from matplotlib.backends.backend_agg import FigureCanvasAgg
        f = plt.figure(figsize=(4, 4))
        _RENDER["fig"], _RENDER["r"] = f, FigureCanvasAgg(f).get_renderer()
    f = _RENDER["fig"]
    t = f.text(0, 0, s, fontsize=size, fontweight=weight, style=style)
    w = t.get_window_extent(renderer=_RENDER["r"]).width / f.dpi
    t.remove()
    return w


PAD = 0.22   # inches of white kept at each figure edge


def minus(s: str) -> str:
    return s.replace("-", "−")


def mark(ax, x, y, scheme, rule, size=None, filled=True, z=4, alpha=1.0):
    col = SCHEME_COL[scheme]
    mk = RULE_MARK[rule]
    ms = size or P["ms"]
    if mk == "D":
        ms *= 0.86
    hollow = (scheme == "sr") or not filled
    ax.plot([x], [y], ls="none", marker=mk, ms=ms, mfc="white" if hollow else col, mec=col,
            mew=P["lw_thin"] * (1.3 if hollow else 0.9), zorder=z, alpha=alpha, clip_on=False)


def key_handles(schemes=("full", "cf", "dc"), rules=True, singly=True):
    h = []
    lab = {"full": "full sample", "cf": "cross-fit", "dc": "double cross-fit"}
    for s in schemes:
        h.append((lab[s], SCHEME_COL[s], "s", True))
    if singly:
        h.append(("singly robust", C["sr"], "o", False))
    if rules:
        h.append(("one-step", C["ink2"], "o", True))
        h.append(("targeting", C["ink2"], "D", True))
    return st.legend_handles(h)


def add_legend(fig, y=0.985, schemes=("full", "cf", "dc"), rules=True, singly=True, x=0.5):
    hs = key_handles(schemes, rules, singly)
    fig.legend(handles=hs, loc="upper center", bbox_to_anchor=(x, y), ncol=len(hs),
               handlelength=0.8, handletextpad=0.35, columnspacing=1.2, borderaxespad=0,
               fontsize=P["fs_small"])


def style_x(ax, ticks=None, fmt=None):
    ax.spines["left"].set_visible(False)
    ax.tick_params(axis="y", length=0, labelleft=False)
    if ticks is not None:
        ax.set_xticks(ticks)
    if fmt:
        ax.set_xticklabels([fmt(t) for t in ticks])
    ax.grid(axis="x", color=C["grid"], lw=P["lw_thin"] * 0.6, zorder=0)
    ax.set_axisbelow(True)


# ---- Figure 3: one score, two choices, one target ------------------------------------
def fig_concept_grid():
    """One shared nuisance engine and one EIF score, the 2 x 3 grid of combination rule (rows,
    marker shape) by fitting scheme (columns, colour), and one target."""
    from matplotlib.patches import FancyBboxPatch
    W, fs = P["width"], P["fs"]
    H = 2.78
    fig = plt.figure(figsize=(W, H))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, W)
    ax.set_ylim(0, H)
    ax.axis("off")
    ink, ink2 = C["ink"], C["ink2"]
    arrow = dict(arrowstyle="-|>", color=ink2, lw=1.0, mutation_scale=7, shrinkA=0, shrinkB=0)

    def box(x, y, w, h, face, edge):
        ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0,rounding_size=0.04",
                                    facecolor=face, edgecolor=edge, lw=0.7))

    # ---- top: one nuisance engine -> one score ----
    top, bh = H - 0.03, 0.60
    box(0.03, top - bh, 1.85, bh, "white", "#C9C9C9")
    ax.text(0.10, top - 0.13, "One Super Learner", fontsize=fs + 0.5, fontweight="bold", va="center", color=ink)
    ax.text(0.10, top - 0.31, r"$\hat g(X)=\hat P(A=1\,|\,X)$", fontsize=fs, va="center", color=ink)
    ax.text(0.10, top - 0.47, r"$\hat Q(a,X)=\hat E[Y\,|\,A=a,X]$", fontsize=fs, va="center", color=ink)
    bx = 2.12
    box(bx, top - bh, W - 0.03 - bx, bh, "#F3F5F8", "#F3F5F8")
    ax.text(bx + 0.07, top - 0.13, "One shared score (the EIF)", fontsize=fs + 0.5, fontweight="bold",
            va="center", color=ink)
    ax.text(bx + 0.07, top - 0.31, r"$D(O)=H(A,X)\,(Y-Q_A)+(Q_1-Q_0)-\psi$", fontsize=fs, va="center", color=ink)
    ax.text(bx + 0.07, top - 0.47, r"every route solves $\mathbb{P}_n\,D=0$", fontsize=fs, va="center", color=ink)
    ax.annotate("", xy=(bx - 0.02, top - bh / 2), xytext=(1.90, top - bh / 2), arrowprops=arrow)
    # ---- grid ----
    gx, cw = 1.62, (W - 0.03 - 1.62) / 3
    hT = top - bh - 0.22
    hB = hT - 0.17
    rT, rh = hB - 0.13, 0.40
    ax.annotate("", xy=(gx + 1.5 * cw, hT + 0.005), xytext=(gx + 1.5 * cw, top - bh - 0.01), arrowprops=arrow)
    ax.text(gx - 0.06, (hT + hB) / 2, "fitting scheme →", fontsize=fs - 0.5, style="italic", ha="right",
            va="center", color=ink2)
    names = {("onestep", 0): "AIPW", ("onestep", 1): "DML", ("onestep", 2): "DC-AIPW",
             ("targeting", 0): "TMLE", ("targeting", 1): "CV-TMLE", ("targeting", 2): "DC-TMLE"}
    for i, (sc, lab) in enumerate([("full", "full sample"), ("cf", "single cross-fit"), ("dc", "double cross-fit")]):
        x = gx + i * cw
        ax.add_patch(Rectangle((x + 0.015, hB), cw - 0.03, hT - hB, color=SCHEME_COL[sc], lw=0))
        ax.text(x + cw / 2, (hT + hB) / 2, lab, color="white", fontweight="bold", fontsize=fs - 0.5,
                ha="center", va="center")
        for j, rule in enumerate(["onestep", "targeting"]):
            y0 = rT - (j + 1) * rh
            ax.add_patch(FancyBboxPatch((x + 0.015, y0 + 0.02), cw - 0.03, rh - 0.04,
                                        boxstyle="round,pad=0,rounding_size=0.03",
                                        facecolor=_shade(SCHEME_COL[sc], 0.88), edgecolor=SCHEME_COL[sc], lw=0.8))
            mark(ax, x + 0.16, y0 + rh / 2, sc, rule, size=P["ms"] * 0.95)
            ax.text(x + 0.27, y0 + rh / 2, names[(rule, i)], fontsize=fs + 0.5, fontweight="bold",
                    ha="left", va="center", color=ink)
    # row labels
    for j, (rule, title, sub) in enumerate([("onestep", "One-step", "plug-in + mean correction"),
                                            ("targeting", "Targeting", "update $\\hat Q$ until the mean\nweighted residual is 0")]):
        y0 = rT - (j + 1) * rh
        mark(ax, 0.09, y0 + rh - 0.10, "full", rule, size=P["ms"] * 0.9)
        ax.lines[-1].set_color(ink2); ax.lines[-1].set_markerfacecolor(ink2); ax.lines[-1].set_markeredgecolor(ink2)
        ax.text(0.17, y0 + rh - 0.10, title, fontsize=fs + 0.5, fontweight="bold", ha="left", va="center", color=ink)
        ax.text(0.04, y0 + rh - 0.20, sub, fontsize=fs - 0.5, ha="left", va="top", color=ink2, linespacing=1.15)
    # DML = AIPW moved one column (over the one-step row)
    ya = rT - 0.005
    ax.plot([gx + 0.5 * cw, gx + 0.5 * cw, gx + 1.5 * cw], [ya, ya + 0.075, ya + 0.075], color=ink, lw=1.0,
            solid_capstyle="butt")
    ax.annotate("", xy=(gx + 1.5 * cw, ya - 0.012), xytext=(gx + 1.5 * cw, ya + 0.08),
                arrowprops=dict(arrowstyle="-|>", color=ink, lw=1.0, mutation_scale=6, shrinkA=0, shrinkB=0))
    # ---- one target ----
    gB = rT - 2 * rh
    for i, sc in enumerate(["full", "cf", "dc"]):
        ax.annotate("", xy=(gx + (i + 0.5) * cw, gB - 0.13), xytext=(gx + (i + 0.5) * cw, gB + 0.005),
                    arrowprops=dict(arrowstyle="-|>", color=SCHEME_COL[sc], lw=1.0, mutation_scale=6,
                                    shrinkA=0, shrinkB=0))
    aT = gB - 0.13
    box(gx + 0.015, aT - 0.25, 3 * cw - 0.03, 0.25, "white", C["navy"])
    ax.text(gx + 1.5 * cw, aT - 0.125, r"one target: the ATE $\psi=E[Y^1-Y^0]$", fontsize=fs + 0.5,
            ha="center", va="center", color=ink)
    ax.text(0.04, aT - 0.02, "For the ATE, DML is AIPW with\n$\\hat g$, $\\hat Q$ fit out of fold: same\nrow, next column. TMLE is\nthe other row.",
            fontsize=fs - 0.5, ha="left", va="top", color=ink2, linespacing=1.15)
    return st.save(fig, OUT / "concept_grid.pdf")


# ---- Figure 4: which rows fit the initial nuisances and which rows they predict --------
def _shade(hex_col, mix):
    """mix > 0: blend toward white by `mix`; mix < 0: darken toward black by -mix."""
    import matplotlib.colors as mc
    r, g, b = mc.to_rgb(hex_col)
    if mix >= 0:
        return (r + (1 - r) * mix, g + (1 - g) * mix, b + (1 - b) * mix)
    return (r * (1 + mix), g * (1 + mix), b * (1 + mix))


def fig_fitting_schemes():
    """Which rows train the initial nuisance fits and which rows they predict, for the full-sample,
    single cross-fit (V = 5) and double cross-fit (V = 3) schemes. Double cross-fitting uses the
    engine's fixed rotation: fold v is predicted with g-hat from fold v+1 and Q-hat from fold v+2
    (mod 3); each fold is drawn once as the held-out fold."""
    W = P["width"]
    fs = P["fs"]
    bar_h, gap = 0.215, 0.025                      # inches
    H = 2.72
    fig = plt.figure(figsize=(W, H))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, W)
    ax.set_ylim(0, H)
    ax.axis("off")
    x0, x1 = 0.03, W - 0.03
    white = dict(color="white", fontweight="bold", fontsize=fs - 0.5, ha="center", va="center")
    dark = dict(color=C["ink"], fontsize=fs - 0.5, ha="center", va="center")
    note = dict(color=C["ink2"], fontsize=fs - 0.5, ha="left", va="top")

    def header(y, col, bold, rest):
        ax.add_patch(Rectangle((x0, y - 0.045), 0.10, 0.09, color=col, lw=0))
        t = ax.text(x0 + 0.15, y, bold, fontsize=fs + 0.5, fontweight="bold", ha="left", va="center", color=C["ink"])
        w = tw(bold, fs + 0.5, "bold")
        ax.text(x0 + 0.15 + w + 0.05, y, rest, fontsize=fs + 0.5, ha="left", va="center", color=C["ink2"])
        return t

    def cell(x, y, w, col, text, style):
        ax.add_patch(Rectangle((x, y), w, bar_h, facecolor=col, edgecolor="white", lw=1.0))
        ax.text(x + w / 2, y + bar_h / 2, text, **style)

    y = H - 0.13
    # ---- full sample ----
    header(y, C["full"], "Full sample", "(AIPW, TMLE)")
    yb = y - 0.13 - bar_h
    cell(x0, yb, x1 - x0, _shade(C["full"], -0.15),
         r"fit initial $\hat g$, $\hat Q$ on all $n$ rows; predict on those same rows", white)
    ax.text(x0, yb - 0.035, "each row helps fit the models that predict it", **note)
    # ---- single cross-fit ----
    y = yb - 0.36
    header(y, C["cf"], "Single cross-fit, $V=5$", "(DML, CV-TMLE)")
    yb = y - 0.13 - bar_h
    fw = (x1 - x0) / 5
    cell(x0, yb, fw, C["cf"], "predict fold $k$", white)
    for i in range(1, 5):
        cell(x0 + i * fw, yb, fw, _shade(C["cf"], 0.6), r"fit $\hat g$, $\hat Q$", dark)
    ax.text(x0, yb - 0.035, r"rotate $k=1,\ldots,5$: initial predictions are out of fold", **note)
    # ---- double cross-fit ----
    y = yb - 0.36
    header(y, C["dc"], "Double cross-fit, $V=3$", "(DC-AIPW, DC-TMLE)")
    lab = 0.34
    fw = (x1 - x0 - lab) / 3
    roles = {"S": (_shade(C["dc"], -0.2), "predict fold $v$", white),
             "G": (_shade(C["dc"], 0.45), r"fit $\hat g$ (fold $v{+}1$)", dark),
             "Q": (_shade(C["dc"], 0.75), r"fit $\hat Q$ (fold $v{+}2$)", dark)}
    yb = y - 0.13 - bar_h
    for v, order in enumerate(["SGQ", "QSG", "GQS"]):
        yy = yb - v * (bar_h + gap)
        ax.text(x0, yy + bar_h / 2, f"$v={v + 1}$", fontsize=fs - 0.5, ha="left", va="center", color=C["ink2"])
        for i, role in enumerate(order):
            col, text, style = roles[role]
            cell(x0 + lab + i * fw, yy, fw, col, text, style)
    yy = yb - 2 * (bar_h + gap)
    ax.text(x0, yy - 0.035, r"fold numbers mod 3: $\hat g$ and $\hat Q$ come from disjoint folds, each a third of the data",
            **note)
    return st.save(fig, OUT / "fitting_schemes.pdf")


# ---- Figure 9: overlap and the spread of the eight estimates ---------------------------
def _dens(ax, x, y, col, lab=None):
    ax.fill_between(x, 0, y, color=col, alpha=0.30, lw=0, zorder=2)
    ax.plot(x, y, color=col, lw=P["lw"] * 0.8, zorder=3, label=lab)


def fig_overlap_contrast():
    F = fd.numbers()
    ov = fd.overlap()
    nh = ov["nhefs"]
    la_ps = ov["diag"][ov["diag"].dataset.str.startswith("LaLonde")]
    W, fs, fsm = P["width"], P["fs"], P["fs_small"]
    H = W * 0.86
    fig = plt.figure(figsize=(W, H))
    names = [n_ for n_, _, _ in fd.ROWS]
    lab_in = max(max(tw(n_, fsm) for n_ in names), tw("one 95% CI (DML)", fsm, style="italic")) + 0.2 + PAD
    x0 = lab_in / W
    gap = 0.75 / W
    wcol = (1 - x0 - PAD / W - gap) / 2
    xs = (x0, x0 + wcol + gap)
    # ---- top row: estimated-propensity densities (NHEFS diagnostic, not seeded; LaLonde diagnostic, seeded) ----
    dy, dh = 0.595, 0.225
    axA = fig.add_axes([xs[0], dy, wcol, dh])
    axB = fig.add_axes([xs[1], dy, wcol, dh])
    for grp, col in [("did not quit", C["control"]), ("quit smoking", C["treated"])]:
        gx, gy = fd.r_density(nh[nh.group == grp].pi)
        _dens(axA, gx, gy, col)
    for grp, col in [("control", C["control"]), ("treated", C["treated"])]:
        gx, gy = fd.r_density(la_ps[la_ps.group == grp].pi)
        _dens(axB, gx, gy, col)
    axB.set_ylim(0, 14)                                  # clips the comparison-group peak (printed below)
    for ax in (axA, axB):
        ax.set_xlim(0, 1)
        ax.set_xticks([0, 0.25, 0.5, 0.75, 1.0])
        ax.set_xticklabels(["0", "0.25", "0.5", "0.75", "1"])
        ax.set_yticks([])
        ax.spines["left"].set_visible(False)
        ax.set_xlabel("estimated propensity score", fontsize=fsm)
    ya = axA.get_ylim()[1]
    axA.text(0.56, ya * 0.50, "did not quit", color=C["control"], ha="left", va="top", fontsize=fsm,
             fontweight="semibold")
    axA.text(0.56, ya * 0.86, "quit", color=C["treated"], ha="left", va="top", fontsize=fsm,
             fontweight="semibold")
    axB.text(0.045, 13.2, "comparison (PSID)", color=C["control"], ha="left", va="top", fontsize=fsm,
             fontweight="semibold")
    axB.text(0.045, 13.2 - fsm / 72 * 14 / (dh * H) * 1.25, f"peak density \u2248 {F['la_ctrl_peak']:.0f}, clipped", color=C["muted"],
             ha="left", va="top", fontsize=fsm * 0.92)
    axB.text(0.72, 5.0, "treated (NSW)", color=C["treated"], ha="center", va="bottom", fontsize=fsm,
             fontweight="semibold")
    for x_, nm, sub, col in [(xs[0], "NHEFS", "good overlap", C["ink2"]),
                             (xs[1], "LaLonde/NSW\u2013PSID", "poor overlap", C["accent"])]:
        fig.text(x_, 0.955, nm, ha="left", va="top", fontsize=P["fs_label"], fontweight="semibold")
        fig.text(x_, 0.955 - P["fs_label"] / 72 * 1.25 / H, sub, ha="left", va="top", fontsize=fs, color=col)
    # ---- bottom row: the eight point estimates on each dataset ----
    axC = fig.add_axes([xs[0], 0.098, wcol, 0.352])
    axD = fig.add_axes([xs[1], 0.098, wcol, 0.352])
    nhp = fd.nhefs_panel()
    lal = fd.lalonde().set_index("name")
    yrow = {n_: 9 - i for i, n_ in enumerate(names)}
    YCI, YLO, YHI = 0.7, -0.1, 9.6
    for ax in (axC, axD):
        ax.set_ylim(YLO, YHI)
    axC.set_xlim(2.3, 4.5)
    style_x(axC, [2.5, 3.0, 3.5, 4.0, 4.5], lambda t: f"{t:.1f}")
    axD.set_xlim(-16000, 2000)
    style_x(axD, [-15000, -10000, -5000, 0], lambda t: minus(f"{t / 1000:,.0f}k") if t else "0")
    axD.axvline(0, color=C["muted"], lw=P["lw_thin"], ls=(0, (4, 3)), zorder=1)
    dml = nhp[(nhp.name == "DML") & (nhp.r == 100)].iloc[0]
    dl = lal.loc["DML"]
    for ax, (a, b) in [(axC, (dml.lo, dml.hi)), (axD, (dl.ci_lo, dl.ci_hi))]:
        ax.plot([a, b], [YCI, YCI], color=C["muted"], lw=P["lw"] * 1.2, solid_capstyle="butt", zorder=3)
        ax.plot([a, a], [YCI - 0.22, YCI + 0.22], color=C["muted"], lw=P["lw_thin"], zorder=3)
        ax.plot([b, b], [YCI - 0.22, YCI + 0.22], color=C["muted"], lw=P["lw_thin"], zorder=3)
    yfig = lambda y: 0.098 + 0.352 * (y - YLO) / (YHI - YLO)
    fig.text(x0 - 0.12 / W, yfig(YCI), "one 95% CI (DML)", ha="right", va="center", fontsize=fsm,
             color=C["muted"], style="italic")
    for n_, sc, ru in fd.ROWS:
        d = nhp[nhp.name == n_]
        v = d[d.r == 100].est.iloc[0] if d.r.notna().any() else d.est.iloc[0]
        mark(axC, v, yrow[n_], sc, ru, size=P["ms"] * 0.85)
        mark(axD, lal.estimate[n_], yrow[n_], sc, ru, size=P["ms"] * 0.85)
        fig.text(x0 - 0.12 / W, yfig(yrow[n_]), n_, ha="right", va="center",
                 fontsize=fsm, color=C["muted"] if sc == "sr" else C["ink"])
    axC.set_xlabel("estimate (kg)", fontsize=fsm)
    axD.set_xlabel("estimate (US$)", fontsize=fsm)
    axC.set_title(f"doubly-robust: {F['nh_dr6_min']:.2f}\u2013{F['nh_dr6_max']:.2f} kg", loc="left",
                  fontsize=fsm, pad=8, fontweight="normal")
    axD.set_title(f"all eight: \u2212\\${-F['la_gcomp']:,.0f} to \u2212\\${-F['la_iptw']:,.0f}", loc="left",
                  fontsize=fsm, pad=8, fontweight="normal")
    add_legend(fig, y=0.995)
    return st.save(fig, OUT / "overlap_estimates.pdf")


# ---- Figure 10: one estimand across implementations and seeds -------------------------
def fig_software():
    F = fd.numbers()
    s = fd.software()
    t = s["table"].set_index("source")
    seeds = s["seeds"]
    hand = t.loc["hand-coded (this paper)"]
    rows = [
        ("dr3: one-step, full sample", hand.AIPW, "full", "onestep", None),
        ("dr3: targeting, full sample", hand.TMLE, "full", "targeting", None),
        ("dr3: one-step, cross-fit", hand.DML_crossfit, "cf", "onestep", seeds.hand_DML),
        ("tmle package (ĝ in sample)", t.loc["tmle (Gruber & van der Laan)"].TMLE, "cf",
         "targeting", seeds.tmle_pkg),
        ("AIPW package", t.loc["AIPW (Zhong et al.)"].AIPW, "cf", "onestep", seeds.AIPW_pkg),
        ("DoubleML", t.loc["DoubleML (Bach et al.)"].DML_crossfit, "cf", "onestep", None),
        ("tmle3", t.loc["tmle3 (tlverse)"].TMLE, "cf", "targeting", None),
    ]
    W, fs = P["width"], P["fs"]
    lab_w = max(tw(r[0], fs) for r in rows) + 0.3 + PAD
    row_in = fs / 72 * 1.6
    top_in, bot_in = 2.45 * fs / 27, 1.3 * fs / 27
    H = len(rows) * row_in + top_in + bot_in + 0.3
    fig = plt.figure(figsize=(W, H))
    ax = fig.add_axes([lab_w / W, bot_in / H, 1 - (lab_w + 0.3 + PAD) / W, 1 - (top_in + bot_in) / H])
    ax.set_xlim(3.20, 3.60)
    n = len(rows)
    ax.set_ylim(-0.6, n - 0.4)
    style_x(ax, [3.2, 3.3, 3.4, 3.5, 3.6], lambda t: f"{t:.1f}")
    ax.set_xlabel("estimate of the NHEFS ATE (kg)")
    med_lo, med_hi = seeds[["hand_DML", "AIPW_pkg", "tmle_pkg"]].median().agg(["min", "max"])
    ax.axvspan(med_lo - 0.002, med_hi + 0.002, color=C["navy"], alpha=0.12, lw=0, zorder=0.5)
    rng = np.random.default_rng(7)
    for i, (lab, v, scheme, rule, sd) in enumerate(rows):
        y = n - 1 - i
        ax.text(3.20 - 0.006, y, lab, ha="right", va="center", fontsize=fs, color=C["ink"])
        if sd is not None:
            jit = rng.uniform(-0.17, 0.17, len(sd))
            ax.scatter(sd, y + jit, s=(P["ms_small"]) ** 2, color=SCHEME_COL[scheme], alpha=0.45,
                       lw=0, zorder=2, marker=RULE_MARK[rule])
            m = np.median(sd)
            ax.plot([m, m], [y - 0.34, y + 0.34], color=C["ink"], lw=P["lw"], zorder=5,
                    solid_capstyle="butt")
        mark(ax, v, y, scheme, rule, size=P["ms"] * 1.1)
    # annotations
    ytop = n - 0.4
    ax.annotate("", xy=(F["sw_min"], ytop + 0.55), xytext=(F["sw_max"], ytop + 0.55),
                arrowprops=dict(arrowstyle="|-|", lw=P["lw_thin"], color=C["ink2"],
                                mutation_scale=fs * 0.45), annotation_clip=False)
    ax.text((F["sw_min"] + F["sw_max"]) / 2, ytop + 0.75,
            f"single runs: {F['sw_min']:.2f}–{F['sw_max']:.2f} kg (range {F['sw_span']:.2f})",
            ha="center", va="bottom", fontsize=P["fs_small"], color=C["ink2"])
    ax.text(med_hi + 0.008, n - 1, f"medians: {F['sw_med']:.2f} kg\n(gap {F['sw_med_gap']:.3f} kg)",
            ha="left", va="center", fontsize=P["fs_small"], color=C["navy"], linespacing=1.05)
    from matplotlib.lines import Line2D
    hs1 = st.legend_handles([("in-sample nuisance fits", C["full"], "s", True),
                                 ("out-of-fold initial outcome fits", C["cf"], "s", True)])
    hs2 = [Line2D([], [], ls="none", marker="o", ms=P["ms"], mfc=C["ink2"], mec=C["ink2"],
                  label="single run"),
           Line2D([], [], ls="none", marker="o", ms=P["ms_small"], mfc=C["ink2"], mec="none", alpha=0.5,
                  label="all 20 seeds"),
           Line2D([], [], ls="none", marker="|", ms=P["ms"] * 1.4, mew=P["lw"], color=C["ink"],
                  label="20-seed median")]
    kw = dict(loc="upper center", handletextpad=0.35, columnspacing=1.3, fontsize=P["fs_small"],
              borderaxespad=0)
    fig.legend(handles=hs1, bbox_to_anchor=(0.5, 0.995), ncol=2, handlelength=0.8, **kw)
    fig.legend(handles=hs2, bbox_to_anchor=(0.5, 0.995 - 0.55 * fs / 27 / H), ncol=3, handlelength=1.0, **kw)
    return st.save(fig, OUT / "software_seeds.pdf")


FIGURES = {"concept_grid": fig_concept_grid, "fitting_schemes": fig_fitting_schemes,
           "overlap_estimates": fig_overlap_contrast, "software_seeds": fig_software}

if __name__ == "__main__":
    ap = argparse.ArgumentParser(description="Draw Figures 3, 4, 9 and 10 as PDFs.")
    ap.add_argument("out", nargs="?", default=str(fd.repo_root() / "output" / "figures"), help="output folder")
    OUT = Path(ap.parse_args().out)
    OUT.mkdir(parents=True, exist_ok=True)
    st.setup()
    for f in FIGURES.values():
        print("wrote", f())
