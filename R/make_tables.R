## =============================================================================
## make_tables.R -- the manuscript's numeric tables from the frozen CSVs, with its rounding and row order:
## Tables 1, 3, 4, 6, 8, 9 and 11, and the numeric entries of Table 10 (the reporting checklist).
## Tables 2, 5, 7 and 12 are text only.
##   Rscript R/make_tables.R [output folder]          (default: output/tables)
## Writes tableNN.csv, one cell per printed cell as plain text, and tableNN.md for the companion page.
## Reads docs/data/ (or the folder in DML_OUT); Table 10 also reads the engine's settings from R/estimators.R.
REPO <- getwd()
while (!file.exists(file.path(REPO, "R", "estimators.R")) && dirname(REPO) != REPO) REPO <- dirname(REPO)
DATA <- Sys.getenv("DML_OUT", file.path(REPO, "docs", "data"))
args <- commandArgs(trailingOnly = TRUE)
TAB  <- if (length(args)) args[1] else file.path(REPO, "output", "tables")
dir.create(TAB, showWarnings = FALSE, recursive = TRUE)
rd <- function(f) read.csv(file.path(DATA, f), stringsAsFactors = FALSE, check.names = FALSE)

## ---- formatting as printed ------------------------------------------------------
f2  <- function(x) sprintf("%.2f", x)
f3  <- function(x) sprintf("%.3f", x)
usd <- function(x, plus = FALSE) paste0(ifelse(plus & round(x) > 0, "+", ""),
                                        formatC(round(x), format = "d", big.mark = ","))
ci2   <- function(lo, hi) sprintf("(%s, %s)", f2(lo), f2(hi))
ciusd <- function(lo, hi) sprintf("(%s, %s)", usd(lo), usd(hi))
le    <- function(x) sprintf("≤ %.2f", ceiling(x * 100) / 100)      # an upper bound at two decimals

write_table <- function(num, df, title) {
  write.csv(df, file.path(TAB, sprintf("table%02d.csv", num)), row.names = FALSE)
  esc <- function(x) gsub("|", "\\|", x, fixed = TRUE)
  md <- c(sprintf("**Table %d.** %s", num, title), "",
          paste0("| ", paste(esc(names(df)), collapse = " | "), " |"),
          paste0("|", paste(rep("---", ncol(df)), collapse = "|"), "|"),
          apply(df, 1, function(r) paste0("| ", paste(esc(r), collapse = " | "), " |")))
  writeLines(md, file.path(TAB, sprintf("table%02d.md", num)))
}

## ---- inputs ----------------------------------------------------------------------
ds  <- rd("data_summary.csv"); nh <- ds[grepl("^NHEFS", ds$dataset), ]
ov  <- rd("nhefs_overlap_data.csv")
p   <- rd("nhefs_panel_table.csv")
gf  <- function(e) p[p$estimator == e & is.na(p$r), ]
gs  <- function(e, r) p[p$estimator == e & !is.na(p$r) & p$r == r, ]
SPL <- c("DML (cross-fit AIPW)", "CV-TMLE", "DC-AIPW", "DC-TMLE")
sw  <- rd("nhefs_software.csv"); swr <- function(pat) sw[grepl(pat, sw$source), ][1, ]
r3  <- rd("results_all3.csv"); la <- r3[grepl("psid", r3$dataset, ignore.case = TRUE), ]
pt  <- rd("positivity_trimming_lalonde.csv")
nr  <- rd("nhefs_reps_raw.csv")
ta  <- rd("nhefs_truncation_audit.csv")            # the analysis fits (Table 3 uses the diagnostic fit)

## the engine's settings, read from R/estimators.R without running a fit
eng <- new.env()
for (x in parse(file.path(REPO, "R", "estimators.R")))
  if (is.call(x) && identical(x[[1]], as.name("<-")) &&
      as.character(x[[2]]) %in% c("default_SL", "nuisance_full", "nuisance_crossfit", "nuisance_double"))
    eval(x, eng)
gbound <- eval(formals(eng$nuisance_full)$gbound)
stopifnot(identical(gbound, eval(formals(eng$nuisance_crossfit)$gbound)),
          identical(gbound, eval(formals(eng$nuisance_double)$gbound)))

## ---- Table 1: the NHEFS analytic sample -------------------------------------------
write_table(1, data.frame(Group = c("Quit smoking (A=1)", "Continued smoking (A=0)", "Total analytic sample"),
                          n = c(nh$n_treated, nh$n - nh$n_treated, nh$n), check.names = FALSE),
            "NHEFS analytic sample.")

## ---- Table 3: NHEFS overlap diagnostics --------------------------------------------
capped <- sum(ov$pi < gbound[1] | ov$pi > gbound[2])
write_table(3, data.frame(
  Diagnostic = c("Sample size n", "Quit smoking (A=1)", "Estimated propensity range", "Truncation window",
                 "Observations with capped propensity"),
  Value = c(nh$n, sprintf("%d (%.1f%%)", nh$n_treated, 100 * nh$n_treated / nh$n),
            sprintf("%s–%s", f3(min(ov$pi)), f3(max(ov$pi))), sprintf("[%s, %s]", f3(gbound[1]), f3(gbound[2])),
            sprintf("%d (%.1f%%)", capped, 100 * capped / nrow(ov))), check.names = FALSE),
  "NHEFS overlap diagnostics for the fitted propensity model.")

## ---- Table 4: what moves the NHEFS estimate ----------------------------------------
label_gap <- max(abs(gf("AIPW (one-step)")$estimate - gf("TMLE")$estimate),
                 abs(gs("DML (cross-fit AIPW)", 100)$estimate - gs("CV-TMLE", 100)$estimate),
                 abs(gs("DC-AIPW", 100)$estimate - gs("DC-TMLE", 100)$estimate))
scheme_move <- max(diff(range(c(gf("AIPW (one-step)")$estimate, gs("DML (cross-fit AIPW)", 100)$estimate,
                                gs("DC-AIPW", 100)$estimate))),
                   diff(range(c(gf("TMLE")$estimate, gs("CV-TMLE", 100)$estimate, gs("DC-TMLE", 100)$estimate))))
dr6 <- c(gf("AIPW (one-step)")$estimate, gf("TMLE")$estimate, sapply(SPL, function(e) gs(e, 100)$estimate))
swv <- c(unlist(swr("^hand")[c("AIPW", "TMLE", "DML_crossfit")]), swr("^tmle \\(")$TMLE, swr("Zhong")$AIPW,
         swr("DoubleML")$DML_crossfit, swr("tmle3")$TMLE)
write_table(4, data.frame(
  `Knob varied` = c("Estimator label (AIPW vs. TMLE)", "Cross-fitting scheme",
                    "Label and scheme together (six estimates)", "Whole software implementation (single runs)"),
  `Held fixed` = c("library, folds, truncation", "label, library, truncation", "library, truncation", "estimand only"),
  `Spread in estimate (kg)` = c(le(label_gap), f2(scheme_move), f2(diff(range(dr6))), f2(diff(range(swv)))),
  check.names = FALSE),
  "What moves the NHEFS estimate and what does not.")

## ---- Table 6: singly-robust baselines and the doubly-robust repair ------------------
t6 <- list(c("g-computation", "g-computation", "outcome model Q"), c("IPTW", "IPTW (stabilized)", "propensity model g"),
           c("AIPW (one-step)", "AIPW (EIF one-step)", "Q or g"), c("TMLE", "TMLE (EIF targeting)", "Q or g"))
write_table(6, do.call(rbind, lapply(t6, function(x) { r <- gf(x[1]); data.frame(
  Estimator = x[2], `Must be correct (for consistency)` = x[3], `Estimate (kg)` = f2(r$estimate),
  `95% CI` = if (is.na(r$ci_lo)) "—" else ci2(r$ci_lo, r$ci_hi), check.names = FALSE) })),
  "Singly-robust baselines and the doubly-robust repair on NHEFS.")

## ---- Table 8: LaLonde/NSW-PSID, every route ------------------------------------------
t8 <- list(c("g-computation", "g-computation", "full-sample"), c("IPTW", "IPTW", "full-sample"),
           c("AIPW", "AIPW", "full-sample"), c("TMLE", "TMLE", "full-sample"),
           c("AIPW = DML", "AIPW = DML", "cross-fit"), c("TMLE (CV-TMLE)", "CV-TMLE", "cross-fit"),
           c("AIPW (DC-AIPW)", "DC-AIPW", "double-cross-fit"), c("TMLE (DC-TMLE)", "DC-TMLE", "double-cross-fit"))
write_table(8, do.call(rbind, lapply(t8, function(x) { r <- la[la$estimator == x[1], ]; stopifnot(nrow(r) == 1)
  data.frame(Estimator = x[2], `Nuisance fitting` = x[3], `Estimate (US$)` = usd(r$estimate), check.names = FALSE) })),
  "Poor overlap on LaLonde/NSW–PSID affects every route.")

## ---- Table 9: LaLonde/NSW-PSID trimming ------------------------------------------------
write_table(9, data.frame(
  `Trimming rule` = c("None (all data)", "ĝ ∈ [0.05, 0.95]", "ĝ ∈ [0.10, 0.90]"),
  `n kept` = usd(pt$n_kept), `% kept` = as.character(pt$pct_kept), `Estimate (US$)` = usd(pt$AIPW_DML, plus = TRUE),
  `95% CI (US$)` = ciusd(pt$ci_lo, pt$ci_hi), check.names = FALSE),
  "Trimming to the estimated common support on LaLonde/NSW–PSID.")

## ---- Table 10: the reporting checklist's numeric entries and settings -------------------------
sdsplit <- tapply(nr$estimate, nr$method, sd)
write_table(10, data.frame(
  Item = c("Estimand", "Estimand", "Super Learner library", "Cross-fitting", "Cross-fitting", "Cross-fitting",
           "Cross-fitting", "Cross-fitting", "Propensity truncation", "Propensity truncation",
           "Trimmed fraction and sensitivity", "Trimmed fraction and sensitivity", "Trimmed fraction and sensitivity"),
  Quantity = c("n", "quitters", "library (default_SL)", "single cross-fit folds V", "double cross-fit folds V",
               "repetitions r", "SD of single-split estimates, cross-fit (kg)",
               "SD of single-split estimates, double cross-fit (kg)", "bounds", "bound on the clever covariate |H|",
               "observations the bounds bind on, over the audited analysis fits", "span of the six doubly-robust estimates at r = 100 (kg)",
               "share of the sample discarded by trimming (Step 6)"),
  Value = c(nh$n, nh$n_treated, paste(eng$default_SL, collapse = ", "),
            formals(eng$nuisance_crossfit)$V, formals(eng$nuisance_double)$V, max(p$r, na.rm = TRUE),
            f2(mean(sdsplit[c("DML (cross-fit AIPW)", "CV-TMLE")])), f2(mean(sdsplit[c("DC-AIPW", "DC-TMLE")])),
            sprintf("[%s, %s]", f3(gbound[1]), f3(gbound[2])), round(1 / gbound[1]), sum(ta$total_clipped),
            f2(diff(range(dr6))), sprintf("%.0f–%.0f%%", 100 * (1 - pt$n_kept[2] / pt$n_kept[1]),
                                          100 * (1 - pt$n_kept[3] / pt$n_kept[1]))),
  check.names = FALSE),
  "Reporting checklist for the NHEFS analysis: numeric entries and settings.")

## ---- Table 11: one estimand, five implementations ----------------------------------------------
h <- swr("^hand")
write_table(11, data.frame(
  Implementation = c("Hand-coded (dr3)", "Hand-coded (dr3)", "Hand-coded (dr3)", "tmle", "AIPW", "DoubleML", "tmle3"),
  `Combination rule` = c("One-step", "Targeting", "One-step", "Targeting", "One-step", "One-step", "Targeting"),
  `Nuisance-fitting scheme` = c("Full sample", "Full sample", "Five-fold cross-fit",
                                "Out-of-fold initial Q̂; full-sample ĝ", "Five-fold cross-fit", "Five-fold cross-fit",
                                "Cross-fit initial fits (CV-TMLE)"),
  `Estimate (kg)` = f2(c(h$AIPW, h$TMLE, h$DML_crossfit, swr("^tmle \\(")$TMLE, swr("Zhong")$AIPW,
                         swr("DoubleML")$DML_crossfit, swr("tmle3")$TMLE)), check.names = FALSE),
  "The same NHEFS estimand across five independent implementations, one row per reported estimate.")

cat("wrote Tables 1, 3, 4, 6, 8, 9, 10 and 11 to", TAB, "\n")
