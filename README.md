# dml-alternatives

Companion code for the tutorial **"Three Routes to One Answer: Reconciling AIPW, TMLE, and Double Machine
Learning for Applied Researchers"** (M. Ehsan Karim; [arXiv:2609.26142](https://arxiv.org/abs/2609.26142)).
AIPW/AIPTW, TMLE and DML are three routes to the *same* efficient influence function for the average treatment
effect. This code builds all of them by hand on one shared Super Learner nuisance library, in full-sample,
cross-fit and double-cross-fit variants, and compares the results with four R packages. The companion R package
[`dr3`](https://github.com/ehsanx/dr3) packages the estimators of `R/estimators.R` behind a single interface.

## Correction

The double cross-fit (DC-AIPW, DC-TMLE) released with arXiv version 1 averaged the predictions of both role
assignments, so the propensity and outcome models were not trained on disjoint folds. `R/estimators.R` now uses
one fixed rotation, as in Zivich and Breskin (2021): fold *v* is scored with the propensity model fit on fold
*v*+1 alone and the outcome model fit on fold *v*+2 alone (mod 3). Every double-cross-fit result was rerun with
the original fold seeds; the other estimators are unchanged. The main corrected values:

- NHEFS, median over 100 splits: DC-AIPW 3.37 kg and DC-TMLE 3.36 kg (were 3.40 and 3.42); the six doubly robust
  estimates now span 3.32–3.39 kg.
- A single double-cross-fit split varies about 1.8 times as much as a single cross-fit split (SD across splits:
  NHEFS 0.14 vs 0.08 kg; RHC 0.0030 / 0.0037 vs 0.0017 / 0.0021). The earlier finding that the two vary about
  equally came from the averaging.
- IHDP: DC-AIPW and DC-TMLE coverage 0.96 and 0.97, RMSE 0.74 and 0.72.
- LaLonde/NSW–PSID, seed 1: DC-AIPW −$881, DC-TMLE −$5,644.

## Contents

| Path | What it holds |
|---|---|
| `R/estimators.R` | The engine: nuisance fitting (`nuisance_full`, `nuisance_crossfit`, `nuisance_double`) and the estimators `est_gcomp`, `est_iptw`, `est_aipw` (one-step) and `est_tmle`. DML is `est_aipw` on cross-fit nuisances; CV-TMLE is `est_tmle` on cross-fit nuisances. |
| `R/tests/test_double_split.R` | Checks the double cross-fit's splitting rule without fitting a model: scored rows are excluded from training, the propensity and outcome models train on disjoint folds, and every row is scored once. |
| `R/utils_data.R` | Data loaders (`load_nhefs`, `load_rhc`, `load_lalonde`, `load_ihdp`, `load_twins`) returning `Y`, `A`, `X` and, for the benchmarks, the true effect. |
| `R/*.R` (analysis) | The analyses behind the manuscript (see the exhibit map) and `R/dcfix/run_dc_fix.R`, which reruns the double-cross-fit rows. |
| `R/figures/`, `python/` | Draw the 10 figures from the frozen CSVs: Figures 1, 2, 5–8 in R; Figures 3, 4, 9 and 10 in Python (`python/paper_figures.py`). |
| `R/make_figures.R`, `R/make_tables.R` | Draw all figures into `output/figures/`; write the numeric tables (1, 3, 4, 6, 8, 9, 11 and Table 10's numeric entries) into `output/tables/`, with the manuscript's rounding and row order. |
| `R/verify_numbers.R` | Numeric audit: re-derives the reported numbers from the frozen CSVs and compares them with the printed values (exit status 1 on any mismatch). Settings the text states, such as the propensity truncation, are read from the parsed scripts. |
| `R/render_page.R`, `docs/` | Build and hold the companion page, <https://ehsanx.github.io/dml-alternatives/>. |
| `docs/data/` | The frozen result CSVs (see "Results and provenance"). |
| `results/dcfix1_213b1138cc/` | The corrected double-cross-fit run: per-study CSVs with run diagnostics (TMLE convergence and score residual, propensity clipping) and `MANIFEST.txt` (engine hash, R and package versions, settings, seeds, timing). |
| `figures/reference/` | The 10 figure PDFs exactly as they appear in the manuscript. |
| `tools/compare_pdf.py` | Compares drawn figure PDFs with `figures/reference/`: page boxes, page structure (drawing order and resources) and a 144-dpi rendering; `--text-only` is a limited, words-only fallback. `tools/test_compare_pdf.py` checks that it catches real changes. |
| `data/checksums.csv` | MD5 checksums, rows and columns of the raw data files (the data themselves are not included). |
| `extras/` | Supplementary analyses, not in the manuscript: `package_crosscheck.R` (the hand-coded estimators against `tmle`, `AIPW` and `DoubleML` on NHEFS, IHDP and Twins), `tmle3_crosscheck.R` (a `tmle3` arm against the hand-coded CV-TMLE), `walkthrough_ihdp.R` (every step on one IHDP realization) and `dr3_equivalence.R` (`dr3` against this engine; see "Checking dr3 against the engine"). |

## Installation

- **R 4.5.1** with the packages in [`renv.lock`](renv.lock): `renv::restore()`. R must have cairo support
  (`capabilities("cairo")`), because the figure scripts use `cairo_pdf`. On Windows, restore into a short folder
  path: under a 150-character path, glmnet's source build stopped (the deepest header paths of renv's staging
  library reached 257–262 characters, against Windows' usual 260-character limit) and renv then installed nothing.
  On Linux, `pdftools` needs the Poppler C++ library (`libpoppler-cpp-dev` on Debian and Ubuntu), and its
  dependencies `qpdf` and `curl` need the development files of libjpeg, zlib and libcurl.
- **Python 3.10** with the packages pinned in [`requirements.txt`](requirements.txt), for Figures 3, 4, 9 and 10
  and the figure comparison: `pip install -r requirements.txt` in a virtual environment. The comparison renders
  with PDFium through `pypdfium2`, whose wheels include the library.
- **Quarto 1.7.31** (the command line, which bundles Pandoc), for the companion page.
- **Arial.** The figures are set in Arial. Without it, the R scripts get a metric-compatible substitute from
  fontconfig (Liberation Sans or Arimo) and `python/paper_figures.py` picks the first of Arial, Liberation Sans and
  Arimo: the text and layout are the same, but the PDFs differ from `figures/reference/`, and only
  `tools/compare_pdf.py --text-only` applies.
- **Locale (Windows).** The reference figures were drawn under the English UTF-8 locale recorded in
  `sessionInfo.txt`. Windows R cannot use `C.UTF-8`: with `LC_ALL=C.UTF-8`, which some shells set, three figures
  differed from the references in their text encoding. Unset `LC_ALL` or set a Windows UTF-8 locale.

## Reproduce the figures, tables and numbers (no raw data)

The figures, the tables, the numeric audit and the companion page read only the frozen CSVs in `docs/data/`; none
of them refits a model. From the repository root:

```bash
Rscript R/tests/test_double_split.R    # splitting rule of the double cross-fit (no model fits)
Rscript R/verify_numbers.R             # numeric audit of the reported numbers
Rscript R/make_figures.R               # the 10 figures, into output/figures/ (R, then Python)
Rscript R/make_tables.R                # the numeric tables, into output/tables/
python tools/compare_pdf.py output/figures figures/reference
Rscript R/render_page.R                # draws and compares the figures, then renders docs/index.html
```

`DML_OUT=output Rscript R/verify_numbers.R` audits freshly generated CSVs instead of the frozen copies; the figure
scripts read the same variable.

### Which script produces which manuscript exhibit

Exhibit numbers follow the manuscript prepared for arXiv version 2 of
[arXiv:2609.26142](https://arxiv.org/abs/2609.26142) (v2 forthcoming). Each row gives the LaTeX label and, where
it changed, the arXiv v1 number. Drawing scripts are in `R/figures/` and
`python/`; the other scripts are in `R/`.

| Manuscript item | Frozen CSV (`docs/data/`) | Data script | Drawn by |
|---|---|---|---|
| Figure 1 — NHEFS DAG (`fig:dag`) | none (schematic) | | `figures/fig_dag.R` |
| Figure 2 — NHEFS propensity overlap (`fig:nhefs-overlap`) | `nhefs_overlap_data.csv` | `nhefs_analysis.R` | `figures/fig_overlap.R` |
| Table 4 — what moves the NHEFS estimate (`tab:knobs`; it now also holds v1 Table 11, the sensitivity summary) | `nhefs_panel_table.csv`, `nhefs_software.csv` | `nhefs_panel_reps.R`, `nhefs_analysis.R` | |
| Figure 3 — one score, two choices, one target (`fig:concept`; v1 Figure 3 and Table 7) | none (schematic) | | `python/paper_figures.py` |
| Figure 4 — fitting schemes (`fig:fitting-schemes`; new) | none (schematic) | | `python/paper_figures.py` |
| Table 6 — NHEFS estimator panel (`tab:eif-repair`) | `nhefs_panel_table.csv` | `nhefs_panel_reps.R` | |
| Figure 5 — NHEFS estimator forest (`fig:nhefs-panel`; v1 Figure 4) | `nhefs_panel_table.csv` | `nhefs_panel_reps.R` | `figures/fig_forests.R` |
| NHEFS split-to-split variation (Step 4) | `nhefs_reps_raw.csv` | `nhefs_panel_reps.R` | |
| Figure 6 — RHC estimator forest (`fig:rhc-panel`; v1 Figure 5) | `real_application.csv`, `rhc_reps_table.csv` | `real_application.R`, `reps_rhc.R` | `figures/fig_forests.R` |
| Figure 7 — IHDP calibration (`fig:ihdp-grading`; v1 Figure 6) | `grade_table.csv`, `grade_raw.csv` | `grade_ihdp.R` | `figures/fig_ihdp.R` |
| IHDP paired contrasts and *r* = 20 check (Step 5) | `ihdp_contrasts.csv`, `ihdp_aggregated_pilot.csv` | `ihdp_contrasts.R` | |
| RHC residual diagnostic (Step 4) | `rhc_residual_diagnostic.csv` | `rhc_residual_diagnostic.R` | |
| Figure 8 — RHC repetition study (`fig:rhc-reps`; v1 Figure 7) | `rhc_reps_table.csv`, `rhc_reps_raw.csv` | `reps_rhc.R` | `figures/fig_reps.R` |
| Figure 9 — overlap and the spread of the estimates (`fig:positivity`; v1 Figure 8, redrawn) | `nhefs_overlap_data.csv`, `positivity_overlap.csv`, `nhefs_panel_table.csv`, `results_all3.csv` | `nhefs_analysis.R`, `positivity_diagnostics.R`, `nhefs_panel_reps.R`, `comparison_all3.R` | `python/paper_figures.py` |
| Figure 10 — one estimand across implementations and seeds (`fig:software-seeds`; new) | `nhefs_software.csv`, `nhefs_matched_repeated.csv` | `nhefs_analysis.R`, `nhefs_matched_repeated.R` | `python/paper_figures.py` |
| Table 8 — LaLonde/NSW–PSID, every route (`tab:positivity-full`; v1 Table 9) | `results_all3.csv`, `lalonde_seed_sweep.csv` | `comparison_all3.R`, `lalonde_seed_sweep.R` | |
| Table 9 — LaLonde/NSW–PSID trimming (`tab:positivity`; v1 Table 10) | `positivity_trimming_lalonde.csv` | `positivity_diagnostics.R` | |
| Truncation audit (Step 1) | `nhefs_truncation_audit.csv` | `dcfix/run_dc_fix.R` (diagnostics) | |
| Table 10 — reporting checklist, numeric entries (`tab:checklist`; v1 Table 12) | `data_summary.csv`, `nhefs_panel_table.csv`, `nhefs_reps_raw.csv`, `nhefs_truncation_audit.csv`, `positivity_trimming_lalonde.csv`, and the engine's settings in `R/estimators.R` | `data_summary.R`, `nhefs_panel_reps.R`, `dcfix/run_dc_fix.R`, `positivity_diagnostics.R` | |
| Table 11 — software comparison (`tab:software`; v1 Table 13) | `nhefs_software.csv` | `nhefs_analysis.R` | |
| Software reconciliation over 20 seeds (Step 7) | `nhefs_matched_repeated.csv` | `nhefs_matched_repeated.R` | |
| Sample sizes and design facts (Tables 1 and 3, Steps 4 to 6) | `data_summary.csv` | `data_summary.R` (raw data) | |
| Numeric audit | all of the above | `verify_numbers.R` | |

## Results and provenance

The CSVs in `docs/data/` are the reference results behind the manuscript. They come from three runs, and two
later files were written on the workstation:

| Run | When | Environment | Engine (`R/estimators.R` md5) | Notes |
|---|---|---|---|---|
| Original analyses | August 2026 | Windows workstation, R 4.5.1 | `bd6c17335e037cb1e9d094f8fec100c7`, the engine released with arXiv version 1 (superseded; not distributed here) | Every row except the double-cross-fit rows and the two later files. |
| Repeated-split and grading jobs (`nhefs_panel_reps.R`, `reps_rhc.R`, `grade_ihdp.R`) | August 2026 | Linux cluster, R 4.4.0 | `bd6c17335e037cb1e9d094f8fec100c7` | As above. |
| Corrected double cross-fitting (`R/dcfix/run_dc_fix.R`) | 2026-10-31 | Windows workstation, R 4.5.1 | `213b1138cc04454b37be4c76108942e0`, as run | The double-cross-fit rows of every file, with the original seeds; outputs and manifest in `results/dcfix1_213b1138cc/`. |
| RHC residual diagnostic (`rhc_residual_diagnostic.R`) | 2026-10-31 | Windows workstation, R 4.5.1 | full-sample and cross-fit fits only, which the correction did not change | `rhc_residual_diagnostic.csv`: the seed-1 RHC fits of `real_application.R`, refitted. |
| Sample sizes and design facts (`data_summary.R`) | 2026-11-01 | Windows workstation, R 4.5.1 | no fits | `data_summary.csv`, from the raw data. |

The correction changed only `nuisance_double` among the estimating code; the engine also gained diagnostics (TMLE
convergence and score residual, propensity clipping) that do not change estimates. In each CSV, the
double-cross-fit (DC-AIPW, DC-TMLE) rows are the corrected run's and every other row is from the original runs.
Twins is supplementary and was not rerun, so `results_all3.csv` has no Twins double-cross-fit rows.
`nhefs_truncation_audit.csv` summarizes the corrected run's clipping diagnostics (`nhefs_trunc.csv` and the
double-cross-fit rows of `nhefs_r1.csv` and `nhefs_reps.csv`).

`R/estimators.R` differs from the engine that made the corrected run only in comments: the run's engine file had
md5 `213b1138cc04454b37be4c76108942e0` (recorded in `MANIFEST.txt`), and both files give the same code hash,
`d93b660d6d0563b944c746df17a48665`, when parsed without comments and deparsed expression by expression. A rerun of
`R/dcfix/run_dc_fix.R` gets a new run id, because the id carries the engine file's md5.

All fits use the Super Learner library SL.mean, SL.glm, SL.glmnet, SL.ranger and SL.earth, with propensities
truncated to [0.025, 0.975], unless noted. The Run column names a run of the table above (original, cluster or
corrected) or the date of a later file.

| File | Contents | Written by | Run | Seeds |
|---|---|---|---|---|
| `nhefs_overlap_data.csv` | NHEFS estimated propensity scores and treatment group, from the overlap diagnostic's own three-learner Super Learner (SL.glm, SL.glmnet, SL.ranger) | `nhefs_analysis.R` | original | **not seeded**: reference file only |
| `nhefs_panel_table.csv`, `nhefs_reps_raw.csv` | NHEFS estimator panel: full-sample fits, and each split-based estimator's running median and interval at *r* = 1 to 100; the single-split estimates | `nhefs_panel_reps.R` | cluster; DC rows corrected | full-sample fits seed 1; splits 100 + *j*, *j* = 1 to 100 |
| `nhefs_software.csv` | One run of each implementation (hand-coded engine, `tmle`, `AIPW`, `DoubleML`, `tmle3`) on NHEFS | `nhefs_analysis.R` | original | seed 1 |
| `nhefs_matched_repeated.csv` | The hand-coded engine, `tmle` and `AIPW` with 20 seeds each | `nhefs_matched_repeated.R` | original | seeds 1 to 20 |
| `nhefs_truncation_audit.csv` | Pre-truncation propensity ranges and clipped units of the NHEFS analysis fits | from `results/dcfix1_213b1138cc/` | corrected | seed 1; splits 100 + *j* |
| `real_application.csv` | RHC and NHEFS single-run panels | `real_application.R` | original; DC rows corrected | seed 1 |
| `rhc_reps_table.csv`, `rhc_reps_raw.csv` | RHC repetition study: running medians and intervals at *r* = 1 to 100; the single-split estimates | `reps_rhc.R` | cluster; DC rows corrected | splits 100 + *j* |
| `rhc_residual_diagnostic.csv` | RHC outcome-model residuals in and out of fold, the g-computation plug-in and AIPW's correction, from the seed-1 fits | `rhc_residual_diagnostic.R` | workstation, 2026-10-31 | seed 1 |
| `grade_raw.csv`, `grade_table.csv` | IHDP grading: eight estimators over 50 realizations and 20 fold partitions each; the table summarizes the raw rows with realization-level Monte Carlo standard errors | `grade_ihdp.R` | cluster; DC rows corrected | fold seeds 101 to 120; **full-sample fits not seeded per unit** |
| `ihdp_contrasts.csv`, `ihdp_aggregated_pilot.csv` | Paired realization-level contrasts between routes, and the exploratory *r* = 20 check | `ihdp_contrasts.R` (no fits) | from `grade_raw.csv` | none |
| `results_all3.csv` | One run of the full panel on IHDP realization 1, LaLonde/NSW–PSID and Twins; `truth` holds each benchmark's known effect, except for LaLonde, where it holds the $1,794 experimental estimate for a different population: LaLonde's `bias` and `covers` are relative to that reference and are not used | `comparison_all3.R` | original; DC rows corrected for IHDP and LaLonde | seed 1 |
| `lalonde_seed_sweep.csv` | LaLonde/NSW–PSID panel over 8 seeds | `lalonde_seed_sweep.R` | original; DC rows corrected | seeds 1 to 8 |
| `positivity_overlap.csv` | Estimated propensity scores and treatment group for IHDP realization 1 and LaLonde/NSW–PSID, from a three-learner Super Learner (SL.glm, SL.glmnet, SL.ranger) | `positivity_diagnostics.R`, part (1) | original | `set.seed(1)`, then the Twins loader (see "Raw data") |
| `positivity_trimming_lalonde.csv` | Cross-fit AIPW (DML) on LaLonde/NSW–PSID after trimming to two common-support regions | `positivity_diagnostics.R`, part (2) | original | seed 1 |
| `data_summary.csv` | Sample sizes, treated units, outcome means and covariate columns of the four datasets; whether the IHDP design is fixed over its 50 realizations | `data_summary.R` (no fits) | from the raw data, 2026-11-01 | none |

The per-unit files (`nhefs_overlap_data.csv`, `positivity_overlap.csv`) hold only fitted propensity scores and
treatment groups: no covariates, outcomes or identifiers.

**Regenerating.** Rerunning a writer script from the raw data reproduces its rows exactly only for seeded fits on
the same package versions and platform. `nhefs_overlap_data.csv` and the full-sample rows of `grade_raw.csv` come
from fits that were not seeded per unit, so a rerun gives slightly different values: rerunning `grade_ihdp.R` for
realization 1 and fold seed 1 on the workstation reproduced both double-cross-fit rows exactly and the other six
rows to within 0.04. The frozen CSVs are therefore the reference.

## Raw data

The analyses read four public datasets and one supplementary benchmark. This section and "Computing
environments" hold the details that the manuscript's Data and Code availability statements refer to. The data
are not redistributed with the manuscript or here: obtain each file as below, place it under `data/` (or the folder in the `DML_DATA_DIR` environment variable) in the
subfolder shown, and run `Rscript R/check_data.R`, which compares every file with
[`data/checksums.csv`](data/checksums.csv) (MD5, rows and columns). The source's terms of use apply to each
dataset; the licence of this repository covers the code only.

| Dataset | Original source | How to obtain it | Files | Loader |
|---|---|---|---|---|
| NHEFS | NHANES I Epidemiologic Followup Study (National Center for Health Statistics), as prepared by Hernán and Robins, *Causal Inference: What If* (2020) | the `nhefs` data of the `causaldata` R package (version 0.1.4 on CRAN); no file needed | none | `load_nhefs()` |
| Right-heart catheterization (RHC) | Connors et al. (1996), *JAMA* 276:889–897 | Vanderbilt's CSV export of the RHC dataset, `rhc.csv`, from the [Right Heart Catheterization Dataset page](https://hbiostat.org/data/repo/rhc) of hbiostat.org (<https://hbiostat.org/data/repo/rhc.csv>). The file analysed is byte-identical to that download (2,320,412 bytes; checked 2026-11-02). The `dr3` loader instead retrieves the R-format export through `Hmisc::getHdata()`, which differs slightly (see `?dr3::dr_rhc`). | `rhc/rhc.csv` (5,735 rows; 63 columns, the first a row index) | `load_rhc()` |
| LaLonde/NSW–PSID | LaLonde (1986), *American Economic Review* 76:604–620; Dehejia and Wahba (1999), *JASA* 94:1053–1062 | the Dehejia–Wahba experimental sample, `nsw_dw.dta`, and the PSID comparison sample, `psid_controls.dta`, from Dehejia's page, <https://users.nber.org/~rdehejia/nswdata2.html> | `lalonde_nsw/nsw_dw.dta` (445 rows), `lalonde_nsw/psid_controls.dta` (2,490 rows) | `load_lalonde("psid")` |
| IHDP (semi-synthetic) | Hill (2011), *Journal of Computational and Graphical Statistics* 20:217–240 | realizations 1–50, `ihdp_npci_1.csv` to `ihdp_npci_50.csv`, from the DragonNet repository, <https://github.com/claudiashi57/dragonnet>, folder `dat/ihdp/csv/` at commit `2e51b24788` (realizations 1–10 are also identical to the CEVAE repository's `datasets/IHDP/csv/`) | `ihdp/ihdp_npci_1.csv` to `ihdp/ihdp_npci_50.csv` (747 rows, 30 columns each, no header: treatment, observed and counterfactual outcomes, the two true conditional means, then 25 covariates) | `load_ihdp()` |
| Twins (supplementary) | Louizos et al. (2017), from the US linked birth and infant death records (Almond, Chay and Lee, 2005) | `data/Twin_data.csv.gz` from the GANITE repository, <https://github.com/jsyoon0823/GANITE>, at commit `e2997317d5` | `twins/Twin_data.csv.gz` (11,400 twin pairs, 32 columns) | `load_twins()` |

**Twins is needed even though the manuscript does not use it.** Part (1) of `R/positivity_diagnostics.R` builds
the Twins sample before its propensity fits, and that loader sets the random-number state the fits use. Without
the Twins file the script cannot regenerate `docs/data/positivity_overlap.csv` exactly; with it, on the
workstation, the script regenerated that file exactly (checked 2026-11-01).

`R/data_summary.R` checks that the 50 IHDP files share one fixed design (the same covariates and treatment in
every realization) and records the sample sizes the manuscript quotes.

## Refitting from the raw data

With the raw data in place, run the analysis scripts from the repository root; they write to `output/`:

```bash
Rscript R/check_data.R               # compare the raw files with data/checksums.csv
Rscript R/data_summary.R             # sample sizes and design facts
Rscript R/nhefs_analysis.R           # NHEFS single-run panel, software table, overlap diagnostic
Rscript R/real_application.R         # RHC and NHEFS full-sample panels
Rscript R/positivity_diagnostics.R   # overlap diagnostic and LaLonde/NSW-PSID trimming
Rscript R/comparison_all3.R          # IHDP realization 1, LaLonde/NSW-PSID and Twins, one run each
Rscript R/lalonde_seed_sweep.R       # LaLonde/NSW-PSID over 8 seeds
Rscript R/nhefs_matched_repeated.R   # three implementations, 20 seeds each
Rscript R/rhc_residual_diagnostic.R  # RHC residual diagnostic
```

The repeated-split and grading jobs (`nhefs_panel_reps.R`, `reps_rhc.R`, `grade_ihdp.R`, `dcfix/run_dc_fix.R`)
are compute-heavy: they run in parallel over `NCORES` workers and resume from per-unit part files.
`R/ihdp_contrasts.R` then computes the IHDP contrasts from `output/grade_raw.csv`.

## Checking dr3 against the engine

[`dr3`](https://github.com/ehsanx/dr3) packages the estimators of `R/estimators.R`. `extras/dr3_equivalence.R`
checks that the two compute the same estimates. It runs both side by side on NHEFS (from `causaldata`, so no raw
data are needed), with the same data, Super Learner library, folds, propensity truncation and random-number
initialization. It compares the full-sample AIPW and TMLE, the seed-1 cross-fit DML and CV-TMLE, and the seed-1
DC-AIPW and DC-TMLE. With `dr3` installed (`remotes::install_github("ehsanx/dr3")`), from the repository root:

```bash
Rscript extras/dr3_equivalence.R
```

The script prints each pair of estimates, their absolute difference and PASS or FAIL against a tolerance of
1e-10, and exits with status 1 if any comparison fails; it takes a few minutes (36 Super Learner fits).
`--seeds 20` adds the cross-fit DML for seeds 2 to 20, Figure 10's seeds (about half an hour). RHC is not used:
the engine reads Vanderbilt's CSV export and `dr3` its R-format export, which differ slightly (see "Raw data").
On the Windows workstation of "Computing environments" (R 4.5.1, SuperLearner 2.0.29, dr3 0.0.0.9001;
2026-11-02), the six pairs of estimates and the nuisance predictions behind them were identical (absolute
difference 0), and the engine reproduced the archived values to 3.6e-15. That is one platform's observation, not
a guarantee of agreement to machine precision elsewhere.

## Computing environments

| Environment | Used for | Recorded in |
|---|---|---|
| Windows workstation, R 4.5.1 | The original analyses other than the repeated-split and grading jobs; the corrected double-cross-fit run; the figures, the page and the audit | [`renv.lock`](renv.lock); [`sessionInfo.txt`](sessionInfo.txt), an edited record without the operating-system build number and time zone |
| Linux cluster, R 4.4.0 | The repeated-split and grading jobs, apart from their double-cross-fit rows | not recorded: the package versions of the cluster's R library were not saved |
| Python 3.10.6 | Figures 3, 4, 9 and 10 and the figure comparison | [`requirements.txt`](requirements.txt) |

Super Learner ensemble weights depend on learner-package versions and on the platform, so exact regeneration of
the CSVs needs the recorded versions on the same platform. Rows from the cluster jobs cannot be guaranteed to
regenerate exactly, because that environment's package versions are unknown. Two fits were not seeded per unit
(see "Regenerating" above): for them, the frozen CSVs are the reference. The numeric audit reads the frozen CSVs
and does not depend on any of this.

## Licence and citation

The code is released under the [MIT licence](LICENSE); the licence does not cover the datasets (see "Raw data").
To cite this work, cite the paper, [arXiv:2609.26142](https://arxiv.org/abs/2609.26142), and this repository,
<https://github.com/ehsanx/dml-alternatives> ([`CITATION.cff`](CITATION.cff)).
