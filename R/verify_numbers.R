## =============================================================================
## verify_numbers.R -- deterministic numeric audit. Re-derives every data-bearing number reported in the
## manuscript (arXiv:2609.26142) from the raw result CSVs and compares it with the printed value:
## PASS = matches at the printed precision; FAIL = flag for review. It reads the raw CSVs, not the table
## generators, so it cannot inherit their errors: docs/data/ (the frozen copies), or the folder named by
## DML_OUT (DML_OUT=output audits freshly generated files). One check (Step 6's clipped propensities)
## reads the corrected double-cross-fit run in results/dcfix1_213b1138cc/ instead. Settings the text states (the propensity
## truncation, the packages' bound, the trimming windows) are read from the parsed scripts, which are not
## run. Exhibit numbers in the labels follow the manuscript prepared for arXiv version 2 (v2 forthcoming); the
## README's exhibit map gives the v1 numbers.
##   Rscript R/verify_numbers.R      (from anywhere in the repository; exit status 1 if any check FAILs)
REPO <- getwd()
while (!file.exists(file.path(REPO, "R", "estimators.R")) && dirname(REPO) != REPO) REPO <- dirname(REPO)
OUT <- Sys.getenv("DML_OUT", file.path(REPO, "docs", "data"))
rd  <- function(f) read.csv(file.path(OUT, f), stringsAsFactors = FALSE, check.names = FALSE)

.rows <- list()
chk <- function(label, printed, computed, digits = 2) {
  ok <- length(computed) == 1 && !is.na(computed) &&
        round(as.numeric(computed), digits) == round(as.numeric(printed), digits)
  .rows[[length(.rows) + 1]] <<- data.frame(
    check = label, printed = printed,
    computed = if (length(computed) == 1 && !is.na(computed)) round(as.numeric(computed), digits) else NA,
    status = if (ok) "PASS" else "FAIL", stringsAsFactors = FALSE)
}

## ---- NHEFS estimator panel: Table 6, Figure 5, Table 4, abstract (nhefs_panel_table.csv) ----
p  <- rd("nhefs_panel_table.csv")
pf <- function(e) p[p$estimator == e & is.na(p$r), ]
ps <- function(e, r) p[p$estimator == e & !is.na(p$r) & p$r == r, ]
gf <- function(e) pf(e)$estimate
gs <- function(e, r) ps(e, r)$estimate
SPL <- c("DML (cross-fit AIPW)", "CV-TMLE", "DC-AIPW", "DC-TMLE")
chk("NHEFS g-computation full (Table 6 = 3.27)", 3.27, gf("g-computation"))
chk("NHEFS IPTW full (Table 6 = 3.19)",          3.19, gf("IPTW"))
chk("NHEFS AIPW full (Table 6 = 3.32)",          3.32, gf("AIPW (one-step)"))
chk("NHEFS TMLE full (Table 6 = 3.33)",          3.33, gf("TMLE"))
chk("NHEFS AIPW full CI lower (2.50)", 2.50, pf("AIPW (one-step)")$ci_lo); chk("NHEFS AIPW full CI upper (4.14)", 4.14, pf("AIPW (one-step)")$ci_hi)
chk("NHEFS TMLE full CI lower (2.50)", 2.50, pf("TMLE")$ci_lo);            chk("NHEFS TMLE full CI upper (4.15)", 4.15, pf("TMLE")$ci_hi)
chk("NHEFS DML r=100 (Fig 5 = 3.38)",     3.38, gs("DML (cross-fit AIPW)", 100))
chk("NHEFS CV-TMLE r=100 (Fig 5 = 3.39)", 3.39, gs("CV-TMLE", 100))
chk("NHEFS DC-AIPW r=100 (Fig 5 = 3.37)", 3.37, gs("DC-AIPW", 100))
chk("NHEFS DC-TMLE r=100 (Fig 5 = 3.36)", 3.36, gs("DC-TMLE", 100))
chk("NHEFS DML r=1 (Fig 5 = 3.45)",       3.45, gs("DML (cross-fit AIPW)", 1))
chk("NHEFS CV-TMLE r=1 (Fig 5 = 3.45)",   3.45, gs("CV-TMLE", 1))
chk("NHEFS DC-AIPW r=1 (Fig 5 = 3.42)",   3.42, gs("DC-AIPW", 1))
chk("NHEFS DC-TMLE r=1 (Fig 5 = 3.42)",   3.42, gs("DC-TMLE", 1))
chk("NHEFS DML r=100 CI lower (2.41)", 2.41, ps("DML (cross-fit AIPW)", 100)$ci_lo)
chk("NHEFS DML r=100 CI upper (4.36)", 4.36, ps("DML (cross-fit AIPW)", 100)$ci_hi)
full2 <- c(gf("AIPW (one-step)"), gf("TMLE"))
dr100 <- c(full2, sapply(SPL, gs, r = 100)); dr1 <- c(full2, sapply(SPL, gs, r = 1))
chk("NHEFS six-DR span r=100 min (3.32)", 3.32, min(dr100)); chk("NHEFS six-DR span r=100 max (3.39)", 3.39, max(dr100))
chk("NHEFS six-DR spread r=100 (0.07)", 0.07, diff(range(dr100)))
chk("NHEFS six-DR spread incl. r=1 rows (0.13)", 0.13, diff(range(dr1)))
chk("NHEFS split-based r=100 min (Step 7 = 3.36)", 3.36, min(sapply(SPL, gs, r = 100)))
chk("NHEFS split-based r=100 max (Step 7 = 3.39)", 3.39, max(sapply(SPL, gs, r = 100)))
chk("NHEFS label swap, largest (Table 4 <= 0.01)", 0.01,
    max(abs(gf("AIPW (one-step)") - gf("TMLE")), abs(gs("DML (cross-fit AIPW)", 100) - gs("CV-TMLE", 100)),
        abs(gs("DC-AIPW", 100) - gs("DC-TMLE", 100))))
chk("NHEFS fold scheme, largest move (Table 4 = 0.07)", 0.07,
    max(diff(range(c(gf("AIPW (one-step)"), gs("DML (cross-fit AIPW)", 100), gs("DC-AIPW", 100)))),
        diff(range(c(gf("TMLE"), gs("CV-TMLE", 100), gs("DC-TMLE", 100))))))
sh <- sapply(SPL, function(e) gs(e, 1) - gs(e, 100))
chk("NHEFS r=1 -> r=100 shift, smallest (0.05)", 0.05, min(sh)); chk("NHEFS r=1 -> r=100 shift, largest (0.06)", 0.06, max(sh))
w <- c(pf("AIPW (one-step)")$ci_width, pf("TMLE")$ci_width, sapply(SPL, function(e) ps(e, 100)$ci_width))
chk("NHEFS CI width, narrowest (1.6)", 1.6, min(w), 1); chk("NHEFS CI width, widest (2.0)", 2.0, max(w), 1)
chk("NHEFS full-sample CI width (1.6)", 1.6, pf("AIPW (one-step)")$ci_width, 1)
chk("NHEFS cross-fit CI width r=1 (1.9)", 1.9, ps("DML (cross-fit AIPW)", 1)$ci_width, 1)
chk("NHEFS double-cross-fit CI width r=1 (2.0)", 2.0, ps("DC-AIPW", 1)$ci_width, 1)
## the text gives one width for each pair of routes; the partners of the three widths above
chk("NHEFS full-sample CI width, TMLE (1.6)", 1.6, pf("TMLE")$ci_width, 1)
chk("NHEFS cross-fit CI width r=1, CV-TMLE (1.9)", 1.9, ps("CV-TMLE", 1)$ci_width, 1)
chk("NHEFS double-cross-fit CI width r=1, DC-TMLE (2.0)", 2.0, ps("DC-TMLE", 1)$ci_width, 1)

## ---- NHEFS repeated splits: SDs across the 100 splits and single-split spans (nhefs_reps_raw.csv) ----
nr <- rd("nhefs_reps_raw.csv"); sdn <- tapply(nr$estimate, nr$method, sd)
chk("NHEFS SD across splits, DML (0.08)", 0.08, sdn[["DML (cross-fit AIPW)"]]); chk("NHEFS SD across splits, CV-TMLE (0.08)", 0.08, sdn[["CV-TMLE"]])
chk("NHEFS SD across splits, DC-AIPW (0.14)", 0.14, sdn[["DC-AIPW"]]);       chk("NHEFS SD across splits, DC-TMLE (0.14)", 0.14, sdn[["DC-TMLE"]])
spans <- tapply(nr$estimate, nr$rep, function(e) diff(range(c(e, full2))))
chk("NHEFS single-split six-DR span, median (0.15)", 0.15, median(spans))
chk("NHEFS single-split span, lower quartile (0.09)", 0.09, unname(quantile(spans, 0.25)))
chk("NHEFS single-split span, upper quartile (0.22)", 0.22, unname(quantile(spans, 0.75)))

## ---- NHEFS software comparison: Table 11 (nhefs_software.csv) ----
sw <- rd("nhefs_software.csv"); swr <- function(pat) sw[grepl(pat, sw$source, ignore.case = TRUE), ][1, ]
h <- swr("hand"); chk("NHEFS hand AIPW (Table 11 = 3.32)", 3.32, h$AIPW)
chk("NHEFS hand TMLE (Table 11 = 3.33)", 3.33, h$TMLE); chk("NHEFS hand DML (Table 11 = 3.41)", 3.41, h$DML_crossfit)
chk("NHEFS tmle-pkg (Table 11 = 3.38)", 3.38, swr("^tmle \\(")$TMLE)
chk("NHEFS AIPW-pkg/Zhong (Table 11 = 3.49)", 3.49, swr("Zhong")$AIPW)
chk("NHEFS DoubleML (Table 11 = 3.46)", 3.46, swr("DoubleML")$DML_crossfit)
chk("NHEFS tmle3 (Table 11 = 3.40)", 3.40, swr("tmle3")$TMLE)
swv <- c(h$AIPW, h$TMLE, h$DML_crossfit, swr("^tmle \\(")$TMLE, swr("Zhong")$AIPW, swr("DoubleML")$DML_crossfit, swr("tmle3")$TMLE)
## range from the unrounded estimates (0.1752 kg), then rounded; the printed endpoints stay 3.32-3.49
chk("NHEFS software spread (Tables 4 and 11 = 0.18)", 0.18, diff(range(swv)))
chk("NHEFS software range, lower end (3.32)", 3.32, min(swv)); chk("NHEFS software range, upper end (3.49)", 3.49, max(swv))

## ---- Software reconciliation: 20 seeds each (nhefs_matched_repeated.csv) ----
mr <- rd("nhefs_matched_repeated.csv"); med <- sapply(mr[c("hand_DML", "AIPW_pkg", "tmle_pkg")], median, na.rm = TRUE)
chk("NHEFS 20-seed median, hand DML (3.39)", 3.39, med[["hand_DML"]])
chk("NHEFS 20-seed median, AIPW pkg (3.39)", 3.39, med[["AIPW_pkg"]])
chk("NHEFS 20-seed median, tmle pkg (3.39)", 3.39, med[["tmle_pkg"]])
chk("NHEFS 20-seed medians, largest gap (0.007)", 0.007, diff(range(med)), 3)
chk("NHEFS SD across seeds, single cross-fit (0.08)", 0.08, sd(mr$hand_DML, na.rm = TRUE))
chk("NHEFS SD across seeds, tmle pkg (0.04)", 0.04, sd(mr$tmle_pkg, na.rm = TRUE))
chk("NHEFS AIPW pkg over seeds, min (3.27)", 3.27, min(mr$AIPW_pkg, na.rm = TRUE))
chk("NHEFS AIPW pkg over seeds, max (3.56)", 3.56, max(mr$AIPW_pkg, na.rm = TRUE))
## Step 7 gaps from the unrounded estimates (the printed endpoints would give 0.10 and 0.08)
chk("NHEFS AIPW pkg first seed above its median (0.11)", 0.11, mr$AIPW_pkg[mr$rep == 1] - med[["AIPW_pkg"]])
chk("NHEFS tmle pkg vs DoubleML single runs, gap (0.075)", 0.075,
    swr("DoubleML")$DML_crossfit - swr("^tmle \\(")$TMLE, 3)
## the full-sample fit sits 0.068-0.075 kg below the three medians (Step 7: "about 0.07")
chk("NHEFS full-sample AIPW below the 20-seed medians, mean gap (0.07)", 0.07, mean(med - h$AIPW))
## the packages' propensity bound in the 20-seed runs, read from the parsed script (our DML there uses the
## engine default, checked below)
code_of <- function(f) paste(unlist(lapply(parse(file.path(REPO, f), keep.source = FALSE), deparse, width.cutoff = 500L)),
                             collapse = "\n")
mrc <- code_of("R/nhefs_matched_repeated.R")
pkb <- as.numeric(sub(".*= ", "", regmatches(mrc, gregexpr("g\\.?bound = [0-9.]+", mrc))[[1]]))
chk("NHEFS 20-seed runs, package calls setting a bound (2)", 2, length(pkb), 0)
chk("NHEFS 20-seed runs, the packages' bound (0.025)", 0.025, unique(pkb), 3)

## ---- Overlap and truncation: Step 1, Table 3, Box 3, Step 7 (nhefs_overlap_data.csv, nhefs_truncation_audit.csv) ----
ov <- rd("nhefs_overlap_data.csv")
chk("NHEFS diagnostic propensity min (0.044)", 0.044, min(ov$pi), 3); chk("NHEFS diagnostic propensity max (0.759)", 0.759, max(ov$pi), 3)
ta <- rd("nhefs_truncation_audit.csv")
chk("NHEFS audited analysis fits: units clipped (0)", 0, sum(ta$total_clipped), 0)
chk("NHEFS audited double-cross-fit runs (101)", 101, ta$fits[grepl("double", ta$fit)], 0)
## the truncation [0.025, 0.975] stated throughout: the default of every engine function that truncates, read
## from the parsed engine without running it, and no engine call in an analysis script that changes it
eng_fn <- Filter(function(e) is.call(e) && identical(e[[1]], as.name("<-")) && is.call(e[[3]]) &&
                   identical(e[[3]][[1]], as.name("function")) && "gbound" %in% names(e[[3]][[2]]) &&
                   !identical(e[[3]][[2]][["gbound"]], quote(expr = )),
                 as.list(parse(file.path(REPO, "R", "estimators.R"), keep.source = FALSE)))
gb <- sapply(eng_fn, function(e) eval(e[[3]][[2]]$gbound))
chk("Engine functions with a truncation bound (3)", 3, length(eng_fn), 0)
chk("Engine propensity truncation, lower bound (0.025)", 0.025, unique(gb[1, ]), 3)
chk("Engine propensity truncation, upper bound (0.975)", 0.975, unique(gb[2, ]), 3)
eng_names <- vapply(eng_fn, function(e) as.character(e[[2]]), "")
n_gb <- function(e) {
  if (!is.call(e)) return(0L)
  f <- if (is.name(e[[1]])) as.character(e[[1]]) else
    if (is.call(e[[1]]) && identical(e[[1]][[1]], as.name("$"))) as.character(e[[1]][[3]]) else ""
  as.integer(f %in% eng_names && "gbound" %in% names(e)) +
    sum(vapply(seq_along(e)[-1], function(i) if (identical(e[[i]], quote(expr = ))) 0L else n_gb(e[[i]]), 0L))
}
scr <- c(list.files(file.path(REPO, "R"), "[.]R$", recursive = TRUE, full.names = TRUE),
         list.files(file.path(REPO, "extras"), "[.]R$", full.names = TRUE))
scr <- scr[!grepl("/tests/", scr) & basename(scr) != "estimators.R"]
chk("Engine calls in analysis scripts that change the bound (0)", 0,
    sum(vapply(scr, function(f) sum(vapply(as.list(parse(f, keep.source = FALSE)), n_gb, 0L)), 0L)), 0)

## ---- IHDP calibration: Figure 7, Step 5, Box 13, Discussion (grade_table.csv) ----
g <- rd("grade_table.csv"); gv <- function(e, c) g[[c]][g$estimator == e]
chk("IHDP AIPW coverage (0.90)", 0.90, gv("AIPW", "coverage"));            chk("IHDP TMLE coverage (0.93)", 0.93, gv("TMLE", "coverage"))
chk("IHDP DML coverage (0.99)", 0.99, gv("AIPW = DML", "coverage"));       chk("IHDP CV-TMLE coverage (0.99)", 0.99, gv("TMLE (CV-TMLE)", "coverage"))
chk("IHDP DC-AIPW coverage (0.96)", 0.96, gv("AIPW (DC-AIPW)", "coverage")); chk("IHDP DC-TMLE coverage (0.97)", 0.97, gv("TMLE (DC-TMLE)", "coverage"))
chk("IHDP AIPW bias (-0.18)", -0.18, gv("AIPW", "mean_bias"));             chk("IHDP TMLE bias (-0.16)", -0.16, gv("TMLE", "mean_bias"))
chk("IHDP DML bias (-0.15)", -0.15, gv("AIPW = DML", "mean_bias"));       chk("IHDP CV-TMLE bias (-0.15)", -0.15, gv("TMLE (CV-TMLE)", "mean_bias"))
chk("IHDP DC-AIPW bias (-0.23)", -0.23, gv("AIPW (DC-AIPW)", "mean_bias")); chk("IHDP DC-TMLE bias (-0.23)", -0.23, gv("TMLE (DC-TMLE)", "mean_bias"))
chk("IHDP g-computation bias (-0.22)", -0.22, gv("g-computation", "mean_bias")); chk("IHDP g-computation RMSE (0.67)", 0.67, gv("g-computation", "rmse"))
fc <- c("AIPW", "TMLE", "AIPW = DML", "TMLE (CV-TMLE)")
chk("IHDP RMSE full/cross-fit, lowest (0.52)", 0.52, min(sapply(fc, gv, c = "rmse")))
chk("IHDP RMSE full/cross-fit, highest (0.56)", 0.56, max(sapply(fc, gv, c = "rmse")))
chk("IHDP DC-AIPW RMSE (0.74)", 0.74, gv("AIPW (DC-AIPW)", "rmse"));      chk("IHDP DC-TMLE RMSE (0.72)", 0.72, gv("TMLE (DC-TMLE)", "rmse"))
chk("IHDP full-sample CI width (1.04)", 1.04, gv("AIPW", "mean_CI_width"))
chk("IHDP DML CI width (1.34)", 1.34, gv("AIPW = DML", "mean_CI_width")); chk("IHDP CV-TMLE CI width (1.33)", 1.33, gv("TMLE (CV-TMLE)", "mean_CI_width"))
chk("IHDP DC-AIPW CI width (1.45)", 1.45, gv("AIPW (DC-AIPW)", "mean_CI_width")); chk("IHDP DC-TMLE CI width (1.44)", 1.44, gv("TMLE (DC-TMLE)", "mean_CI_width"))
chk("IHDP cross-fit vs full width, % wider (30)", 30, 100 * (gv("AIPW = DML", "mean_CI_width") / gv("AIPW", "mean_CI_width") - 1), 0)
chk("IHDP AIPW coverage MCSE (0.035)", 0.035, gv("AIPW", "cov_mcse"), 3);  chk("IHDP TMLE coverage MCSE (0.032)", 0.032, gv("TMLE", "cov_mcse"), 3)
chk("IHDP cross-fit coverage MCSE (0.005)", 0.005, gv("AIPW = DML", "cov_mcse"), 3)
chk("IHDP DC-AIPW coverage MCSE (0.015)", 0.015, gv("AIPW (DC-AIPW)", "cov_mcse"), 3); chk("IHDP DC-TMLE coverage MCSE (0.012)", 0.012, gv("TMLE (DC-TMLE)", "cov_mcse"), 3)
chk("IHDP AIPW coverage below 0.95, in MCSEs (1.5)", 1.5, (0.95 - gv("AIPW", "coverage")) / gv("AIPW", "cov_mcse"), 1)
chk("IHDP DML bias MCSE (0.065)", 0.065, gv("AIPW = DML", "bias_mcse"), 3); chk("IHDP AIPW bias MCSE (0.075)", 0.075, gv("AIPW", "bias_mcse"), 3)
bm <- sapply(c(fc, "AIPW (DC-AIPW)", "TMLE (DC-TMLE)"), gv, c = "bias_mcse")
chk("IHDP bias MCSE, smallest (0.06)", 0.06, min(bm)); chk("IHDP bias MCSE, largest (0.08)", 0.08, max(bm))

## ---- IHDP paired contrasts and the r = 20 pilot (grade_raw.csv, ihdp_contrasts.csv, ihdp_aggregated_pilot.csv) ----
gr <- rd("grade_raw.csv"); tr <- tapply(gr$truth, gr$realization, unique)
chk("IHDP sample-average target, min (3.16)", 3.16, min(tr)); chk("IHDP sample-average target, max (12.65)", 12.65, max(tr))
chk("IHDP realizations (50)", 50, length(tr), 0); chk("IHDP fold seeds per realization (20)", 20, length(unique(gr$seed)), 0)
ct <- rd("ihdp_contrasts.csv"); cv <- function(a, b, q, col = "difference") ct[[col]][ct$from == a & ct$to == b & ct$quantity == q]
## Step 5's small-pilot example: realizations 1-10 alone
p10 <- function(e) { s <- gr[gr$estimator == e & gr$realization <= 10, ]; tapply(as.numeric(s$covers), s$realization, mean) }
d10 <- p10("AIPW = DML") - p10("AIPW")
chk("IHDP realizations 1-10: AIPW->DML coverage gap (0.275)", 0.275, mean(d10), 3)
chk("IHDP realizations 1-10: paired MCSE (0.123)", 0.123, sd(d10) / sqrt(length(d10)), 3)
chk("IHDP AIPW->DML coverage gap (0.09)", 0.09, cv("AIPW", "AIPW = DML", "coverage"));   chk("IHDP AIPW->DML paired t (2.8)", 2.8, cv("AIPW", "AIPW = DML", "coverage", "t"), 1)
chk("IHDP TMLE->CV-TMLE coverage gap (0.06)", 0.06, cv("TMLE", "TMLE (CV-TMLE)", "coverage")); chk("IHDP TMLE->CV-TMLE paired t (2.0)", 2.0, cv("TMLE", "TMLE (CV-TMLE)", "coverage", "t"), 1)
chk("IHDP DC-AIPW vs DML coverage (-0.035)", -0.035, cv("AIPW = DML", "AIPW (DC-AIPW)", "coverage"), 3)
chk("IHDP DC-AIPW vs DML paired t (-2.9)", -2.9, cv("AIPW = DML", "AIPW (DC-AIPW)", "coverage", "t"), 1)
chk("IHDP DC-TMLE vs CV-TMLE coverage (-0.026)", -0.026, cv("TMLE (CV-TMLE)", "TMLE (DC-TMLE)", "coverage"), 3)
chk("IHDP DC-TMLE vs CV-TMLE paired t (-2.5)", -2.5, cv("TMLE (CV-TMLE)", "TMLE (DC-TMLE)", "coverage", "t"), 1)
chk("IHDP DC-AIPW vs AIPW coverage t (1.7)", 1.7, cv("AIPW", "AIPW (DC-AIPW)", "coverage", "t"), 1)
chk("IHDP DC-TMLE vs TMLE coverage t (1.2)", 1.2, cv("TMLE", "TMLE (DC-TMLE)", "coverage", "t"), 1)
chk("IHDP DC MSE vs AIPW, paired t (2.0)", 2.0, cv("AIPW", "AIPW (DC-AIPW)", "mse", "t"), 1)
chk("IHDP DC MSE vs TMLE, paired t (2.1)", 2.1, cv("TMLE", "TMLE (DC-TMLE)", "mse", "t"), 1)
chk("IHDP DC MSE vs DML, paired t (1.8)", 1.8, cv("AIPW = DML", "AIPW (DC-AIPW)", "mse", "t"), 1)
chk("IHDP DC MSE vs CV-TMLE, paired t (1.6)", 1.6, cv("TMLE (CV-TMLE)", "TMLE (DC-TMLE)", "mse", "t"), 1)
chk("IHDP DML vs AIPW bias, paired t (1.2)", 1.2, cv("AIPW", "AIPW = DML", "bias", "t"), 1)
## all 15 pairwise bias contrasts among the six doubly-robust routes, from the raw grading
dr6 <- c("AIPW", "TMLE", "AIPW = DML", "TMLE (CV-TMLE)", "AIPW (DC-AIPW)", "TMLE (DC-TMLE)")
pb  <- sapply(dr6, function(e) tapply(gr$bias[gr$estimator == e], gr$realization[gr$estimator == e], mean))
bt  <- do.call(rbind, lapply(1:5, function(i) do.call(rbind, lapply((i + 1):6, function(j) {
  d <- pb[, j] - pb[, i]; data.frame(a = dr6[i], b = dr6[j], diff = mean(d), t = mean(d) / (sd(d) / sqrt(length(d)))) }))))
bres <- abs(bt$t) >= 2
chk("IHDP resolved bias pairs, count (2)", 2, sum(bres), 0)
chk("IHDP TMLE less biased than AIPW, by (0.023)", 0.023, bt$diff[bt$a == "AIPW" & bt$b == "TMLE"], 3)
chk("IHDP TMLE vs AIPW bias, paired t (3.2)", 3.2, bt$t[bt$a == "AIPW" & bt$b == "TMLE"], 1)
chk("IHDP CV-TMLE less biased than DML, by (0.006)", 0.006, bt$diff[bt$a == "AIPW = DML" & bt$b == "TMLE (CV-TMLE)"], 3)
chk("IHDP CV-TMLE vs DML bias, paired t (2.8)", 2.8, bt$t[bt$a == "AIPW = DML" & bt$b == "TMLE (CV-TMLE)"], 1)
chk("IHDP other bias pairs, largest |t| below 1.9 (1.8)", 1.8, max(abs(bt$t[!bres])), 1)
chk("IHDP mean bias, least negative DR (-0.15)", -0.15, max(colMeans(pb))); chk("IHDP mean bias, most negative DR (-0.23)", -0.23, min(colMeans(pb)))
pil <- rd("ihdp_aggregated_pilot.csv")
for (e in pil$estimator) chk(sprintf("IHDP r=20 aggregated, %s included (50)", e), 50, pil$included[pil$estimator == e], 0)
chk("IHDP 50/50 Clopper-Pearson lower limit (0.93)", 0.93, qbeta(0.025, 50, 1))

## ---- RHC panel: Figure 6, Step 4 (real_application.csv, rhc_reps_table.csv, rhc_reps_raw.csv) ----
ra <- rd("real_application.csv"); rf <- function(e) ra[grepl("^RHC", ra$dataset) & ra$estimator == e, ]
chk("RHC g-computation (0.026)", 0.026, rf("g-computation")$estimate, 3); chk("RHC IPTW (0.063)", 0.063, rf("IPTW")$estimate, 3)
chk("RHC AIPW full (0.033)", 0.033, rf("AIPW")$estimate, 3);              chk("RHC TMLE full (0.053)", 0.053, rf("TMLE")$estimate, 3)
chk("RHC seed-1 DML (Step 4 = 0.044)", 0.044, rf("AIPW = DML")$estimate, 3); chk("RHC seed-1 CV-TMLE (0.047)", 0.047, rf("TMLE (CV-TMLE)")$estimate, 3)
rr <- rd("rhc_reps_table.csv"); rm_ <- function(m, r) rr[rr$method == m & rr$r == r, ]
RSP <- c("DML (cross-fit)", "CV-TMLE (cross-fit)", "DC-AIPW", "DC-TMLE")
chk("RHC DML r=100 (0.043)", 0.043, rm_("DML (cross-fit)", 100)$estimate, 3); chk("RHC CV-TMLE r=100 (0.047)", 0.047, rm_("CV-TMLE (cross-fit)", 100)$estimate, 3)
chk("RHC DC-AIPW r=100 (0.045)", 0.045, rm_("DC-AIPW", 100)$estimate, 3); chk("RHC DC-TMLE r=100 (0.049)", 0.049, rm_("DC-TMLE", 100)$estimate, 3)
e100 <- sapply(RSP, function(m) rm_(m, 100)$estimate); e1 <- sapply(RSP, function(m) rm_(m, 1)$estimate)
chk("RHC split-based spread r=100 (0.006)", 0.006, diff(range(e100)), 3); chk("RHC split-based spread r=1 (0.007)", 0.007, diff(range(e1)), 3)
chk("RHC split-based r=1, max (0.050)", 0.050, max(e1), 3)
chk("RHC DC vs single cross-fit at r=100, largest gap (0.002)", 0.002,
    max(abs(e100[["DC-AIPW"]] - e100[["DML (cross-fit)"]]), abs(e100[["DC-TMLE"]] - e100[["CV-TMLE (cross-fit)"]])), 3)
se100 <- sapply(RSP, function(m) rm_(m, 100)$se)
chk("RHC SE r=100, smallest (0.012)", 0.012, min(se100), 3); chk("RHC SE r=100, largest (0.013)", 0.013, max(se100), 3)
lo <- c(sapply(RSP, function(m) rm_(m, 1)$ci_lo), sapply(RSP, function(m) rm_(m, 100)$ci_lo), rf("AIPW")$ci_lo, rf("TMLE")$ci_lo)
chk("RHC every interval excludes zero (all lower limits > 0)", 1, as.numeric(all(lo > 0)), 0)
rw <- rd("rhc_reps_raw.csv"); sdr <- tapply(rw$estimate, rw$method, sd)
chk("RHC SD across splits, DML (0.0017)", 0.0017, sdr[["DML (cross-fit)"]], 4); chk("RHC SD across splits, CV-TMLE (0.0021)", 0.0021, sdr[["CV-TMLE (cross-fit)"]], 4)
chk("RHC SD across splits, DC-AIPW (0.0030)", 0.0030, sdr[["DC-AIPW"]], 4);     chk("RHC SD across splits, DC-TMLE (0.0037)", 0.0037, sdr[["DC-TMLE"]], 4)
d <- rw$estimate[rw$method == "DML (cross-fit)"]
chk("RHC DML single split, min (0.038)", 0.038, min(d), 3); chk("RHC DML single split, max (0.047)", 0.047, max(d), 3)
chk("RHC DML single-split swing (0.0095)", 0.0095, diff(range(d)), 4)
## running medians: first r after which every later value stays within 0.001 of the r = 100 value
settle <- function(m, tol = 0.001) { s <- rr[rr$method == m, ]; s <- s[order(s$r), ]; ok <- abs(s$estimate - s$estimate[s$r == 100]) < tol
  for (k in seq_along(ok)) if (all(ok[k:length(ok)])) return(s$r[k]); NA }
chk("RHC DML settles within 0.001 from r = (7)", 7, settle("DML (cross-fit)"), 0);  chk("RHC CV-TMLE settles within 0.001 from r = (12)", 12, settle("CV-TMLE (cross-fit)"), 0)
chk("RHC DC-AIPW settles within 0.001 from r = (54)", 54, settle("DC-AIPW"), 0);                 chk("RHC DC-TMLE settles within 0.001 from r = (41)", 41, settle("DC-TMLE"), 0)
exc <- function(m) { s <- rr[rr$method == m, ]; max(abs(s$estimate - s$estimate[s$r == 100])) }
chk("RHC DML largest distance from r=100 value (0.0014)", 0.0014, exc("DML (cross-fit)"), 4)
chk("RHC CV-TMLE largest distance (0.0016)", 0.0016, exc("CV-TMLE (cross-fit)"), 4)
chk("RHC DC-AIPW largest distance (0.0029)", 0.0029, exc("DC-AIPW"), 4); chk("RHC DC-TMLE largest distance (0.0041)", 0.0041, exc("DC-TMLE"), 4)
chk("RHC net change r=1 -> r=100, largest (0.0009)", 0.0009, max(abs(e100 - e1)), 4)
## at r = 8 both double-cross-fit medians sat below their r = 100 values (Step 4)
dc8 <- sapply(c("DC-AIPW", "DC-TMLE"), function(m) rm_(m, 100)$estimate - rm_(m, 8)$estimate)
chk("RHC DC medians at r=8, both below their r=100 values", 1, as.numeric(all(dc8 > 0)), 0)
chk("RHC DC medians at r=8, largest shortfall (0.004)", 0.004, max(dc8), 3)
gw <- sapply(RSP, function(m) 100 * (rm_(m, 100)$ci_width / rm_(m, 1)$ci_width - 1))
chk("RHC CI width growth, DML % (1.5)", 1.5, gw[["DML (cross-fit)"]], 1);  chk("RHC CI width growth, CV-TMLE % (1.7)", 1.7, gw[["CV-TMLE (cross-fit)"]], 1)
chk("RHC CI width growth, DC-AIPW % (4.6)", 4.6, gw[["DC-AIPW"]], 1);       chk("RHC CI width growth, DC-TMLE % (6.8)", 6.8, gw[["DC-TMLE"]], 1)
## estimator-label gap relative to CI width: reported estimates (Figures 5-6) and single RHC splits
lab_ratio <- function(tab, key, a, b, r) { x <- tab[tab[[key]] == a & tab$r == r, ]; y <- tab[tab[[key]] == b & tab$r == r, ]
  abs(x$estimate - y$estimate) / ((x$ci_width + y$ci_width) / 2) }
rep_rat <- c(sapply(c(1, 100), function(r) c(lab_ratio(rr, "method", "DML (cross-fit)", "CV-TMLE (cross-fit)", r), lab_ratio(rr, "method", "DC-AIPW", "DC-TMLE", r))),
             sapply(c(1, 100), function(r) c(lab_ratio(p, "estimator", "DML (cross-fit AIPW)", "CV-TMLE", r), lab_ratio(p, "estimator", "DC-AIPW", "DC-TMLE", r))))
chk("Label gap / CI width, reported estimates, max (0.1)", 0.1, max(rep_rat), 1)
wr_ <- reshape(rw[, c("rep", "method", "estimate", "se")], idvar = "rep", timevar = "method", direction = "wide")
chk("Label gap / CI width, single RHC split, max (0.15)", 0.15,
    max(abs(wr_[["estimate.DC-TMLE"]] - wr_[["estimate.DC-AIPW"]]) / (2 * 1.96 * (wr_[["se.DC-AIPW"]] + wr_[["se.DC-TMLE"]]) / 2)))

## ---- RHC residual diagnostic: Step 4 (rhc_residual_diagnostic.csv) ----
rd_ <- rd("rhc_residual_diagnostic.csv"); dv <- function(s, c) rd_[[c]][grepl(s, rd_$scheme)]
chk("RHC residual RMS, in-sample (0.30)", 0.30, dv("full", "resid_rms"));      chk("RHC residual RMS, out of fold (0.42)", 0.42, dv("cross", "resid_rms"))
chk("RHC AIPW correction, in-sample (0.007)", 0.007, dv("full", "correction"), 3)
chk("RHC AIPW correction, out of fold (0.014)", 0.014, dv("cross", "correction"), 3)
chk("RHC g-computation plug-in, in-sample (0.026)", 0.026, dv("full", "gcomp"), 3)
chk("RHC diagnostic AIPW = Figure 6 full-sample AIPW", 1, as.numeric(abs(dv("full", "aipw") - rf("AIPW")$estimate) < 1e-8), 0)
chk("RHC diagnostic cross-fit AIPW = seed-1 DML", 1,
    as.numeric(abs(dv("cross", "aipw") - ra$estimate[grepl("^RHC", ra$dataset) & ra$estimator == "AIPW = DML"]) < 1e-8), 0)

## ---- Poor overlap: Table 8 and Step 6 (results_all3.csv, lalonde_seed_sweep.csv) ----
r3 <- rd("results_all3.csv"); lw <- function(e) r3[grepl("psid", r3$dataset, ignore.case = TRUE) & r3$estimator == e, ]
lv <- function(e) lw(e)$estimate
chk("LaLonde g-computation (Table 8 = -819)", -819, lv("g-computation"), 0); chk("LaLonde IPTW (Table 8 = -14629)", -14629, lv("IPTW"), 0)
chk("LaLonde AIPW full (Table 8 = -833)", -833, lv("AIPW"), 0);              chk("LaLonde TMLE full (Table 8 = -1254)", -1254, lv("TMLE"), 0)
chk("LaLonde DML (Table 8 = -1284)", -1284, lv("AIPW = DML"), 0);            chk("LaLonde CV-TMLE (Table 8 = -5186)", -5186, lv("TMLE (CV-TMLE)"), 0)
chk("LaLonde DC-AIPW (Table 8 = -881)", -881, lv("AIPW (DC-AIPW)"), 0);      chk("LaLonde DC-TMLE (Table 8 = -5644)", -5644, lv("TMLE (DC-TMLE)"), 0)
chk("LaLonde DML CI lower (-2278)", -2278, lw("AIPW = DML")$ci_lo, 0);       chk("LaLonde DML CI upper (-291)", -291, lw("AIPW = DML")$ci_hi, 0)
sw8 <- rd("lalonde_seed_sweep.csv"); sv <- function(e) sw8$estimate[sw8$estimator == e]
chk("LaLonde CV-TMLE 8-seed median (-5077)", -5077, median(sv("TMLE (CV-TMLE)")), 0)
chk("LaLonde CV-TMLE 8-seed range, nearest zero (-3641)", -3641, max(sv("TMLE (CV-TMLE)")), 0)
chk("LaLonde CV-TMLE 8-seed range, farthest (-5384)", -5384, min(sv("TMLE (CV-TMLE)")), 0)
chk("LaLonde DC-TMLE 8-seed median (-4830)", -4830, median(sv("TMLE (DC-TMLE)")), 0)
chk("LaLonde DC-TMLE 8-seed range, nearest zero (-4125)", -4125, max(sv("TMLE (DC-TMLE)")), 0)
chk("LaLonde DC-TMLE 8-seed range, farthest (-5644)", -5644, min(sv("TMLE (DC-TMLE)")), 0)
chk("LaLonde every DR estimate negative, all seeds", 1,
    as.numeric(all(sw8$estimate[!sw8$estimator %in% c("g-computation", "IPTW")] < 0)), 0)
## Step 6: the propensity bound in the reported double-cross-fit run. The seed-1 DC rows of the corrected run
## (results/dcfix1_213b1138cc/, read directly), which must reproduce Table 8's DC estimates and agree on the count.
dcl <- read.csv(file.path(REPO, "results", "dcfix1_213b1138cc", "lalonde_r1.csv"), stringsAsFactors = FALSE)
dcl <- dcl[dcl$study == "lalonde_r1" & dcl$seed == 1 & dcl$method %in% c("DC-AIPW", "DC-TMLE"), ]
dcl_ok <- nrow(dcl) == 2 && setequal(round(dcl$estimate), c(-881, -5644)) && length(unique(dcl$n_clipped)) == 1
chk("LaLonde DC run, seed 1: propensities clipped (Step 6 = 2093)", 2093, if (dcl_ok) dcl$n_clipped[1] else NA, 0)

## ---- Trimming: Table 9 (positivity_trimming_lalonde.csv) ----
pt <- rd("positivity_trimming_lalonde.csv")
chk("LaLonde none: DML (-1284)", -1284, pt$AIPW_DML[1], 0)
chk("LaLonde none: CI lower (Table 9 = -2278)", -2278, pt$ci_lo[1], 0); chk("LaLonde none: CI upper (Table 9 = -291)", -291, pt$ci_hi[1], 0)
chk("LaLonde trim[.05,.95] est (+484)", 484, pt$AIPW_DML[2], 0); chk("LaLonde trim[.05,.95] % kept (16)", 16, pt$pct_kept[2], 0)
chk("LaLonde trim[.10,.90] est (-27)", -27, pt$AIPW_DML[3], 0);  chk("LaLonde trim[.10,.90] % kept (11)", 11, pt$pct_kept[3], 0)
chk("LaLonde trim[.05,.95] CI lower (-1229)", -1229, pt$ci_lo[2], 0); chk("LaLonde trim[.05,.95] CI upper (2196)", 2196, pt$ci_hi[2], 0)
chk("LaLonde trim[.10,.90] CI lower (-2214)", -2214, pt$ci_lo[3], 0); chk("LaLonde trim[.10,.90] CI upper (2161)", 2161, pt$ci_hi[3], 0)
## the windows of Table 9, read from the result labels; each label names the window its rule applies
## (rules in R/positivity_diagnostics.R, parsed without running)
wl <- lapply(regmatches(pt$rule[2:3], gregexpr("[0-9.]+", pt$rule[2:3])), as.numeric)
chk("LaLonde trimming window 1, lower (Table 9 = 0.05)", 0.05, wl[[1]][1]); chk("LaLonde trimming window 1, upper (Table 9 = 0.95)", 0.95, wl[[1]][2])
chk("LaLonde trimming window 2, lower (Table 9 = 0.10)", 0.10, wl[[2]][1]); chk("LaLonde trimming window 2, upper (Table 9 = 0.90)", 0.90, wl[[2]][2])
pdc <- paste(unlist(lapply(parse(file.path(REPO, "R", "positivity_diagnostics.R"), keep.source = FALSE), deparse,
                           width.cutoff = 500L)), collapse = "\n")
tw <- lapply(regmatches(pdc, gregexpr("`trim to \\[[0-9.]+, [0-9.]+\\]` = pi >= [0-9.]+ & pi <= [0-9.]+", pdc))[[1]],
             function(x) as.numeric(regmatches(x, gregexpr("[0-9.]+", x))[[1]]))
chk("LaLonde trimming rules, label = window applied (2 of 2)", 2, sum(vapply(tw, function(x) all(x[1:2] == x[3:4]), NA)), 0)

## ---- Sample sizes and design facts from the raw data: Tables 1 and 3, Steps 4 and 6, Figure 7 ----
## (data_summary.csv, written by R/data_summary.R from the raw data)
ds <- rd("data_summary.csv"); dsr <- function(p) ds[grepl(p, ds$dataset), ]
chk("NHEFS analytic sample (Tables 1 and 3 = 1566)", 1566, dsr("^NHEFS")$n, 0)
chk("NHEFS quit smoking (Table 1 = 403)", 403, dsr("^NHEFS")$n_treated, 0)
chk("NHEFS continued smoking (Table 1 = 1163)", 1163, dsr("^NHEFS")$n - dsr("^NHEFS")$n_treated, 0)
chk("NHEFS quit smoking, % (Table 3 = 25.7)", 25.7, 100 * dsr("^NHEFS")$n_treated / dsr("^NHEFS")$n, 1)
chk("NHEFS diagnostic scores capped (Table 3 = 0)", 0, sum(ov$pi < 0.025 | ov$pi > 0.975), 0)
chk("RHC sample size (Figure 6, Step 4 = 5735)", 5735, dsr("^RHC")$n, 0)
chk("RHC 30-day mortality, % (Step 4 = 33)", 33, 100 * dsr("^RHC")$outcome_mean, 0)
chk("RHC one-hot covariate columns within 60-70", 1,
    as.numeric(dsr("^RHC")$n_covariate_columns >= 60 & dsr("^RHC")$n_covariate_columns <= 70), 0)
chk("LaLonde/NSW-PSID units (Table 9 = 2675)", 2675, dsr("^LaLonde")$n, 0)
chk("LaLonde/NSW treated units (185)", 185, dsr("^LaLonde")$n_treated, 0)
chk("LaLonde/PSID comparison units (2490)", 2490, dsr("^LaLonde")$n - dsr("^LaLonde")$n_treated, 0)
chk("IHDP sample size (Figure 7, Step 5 = 747)", 747, dsr("^IHDP")$n, 0)
chk("IHDP covariates and treatment fixed over 50 realizations", 1, as.numeric(isTRUE(as.logical(dsr("^IHDP")$design_fixed_all_50))), 0)

## ---- LaLonde reference, trimming sizes and half-widths: Tables 8 and 9, Step 6 ----
chk("LaLonde experimental benchmark (Table 8 = 1794)", 1794, unique(lw("AIPW")$truth), 0)
chk("LaLonde none: n kept (Table 9 = 2675)", 2675, pt$n_kept[1], 0)
chk("LaLonde trim[.05,.95] n kept (Table 9 = 423)", 423, pt$n_kept[2], 0)
chk("LaLonde trim[.10,.90] n kept (Table 9 = 300)", 300, pt$n_kept[3], 0)
chk("LaLonde trimming discards, fewest % (84)", 84, 100 * (1 - pt$n_kept[2] / pt$n_kept[1]), 0)
chk("LaLonde trimming discards, most % (89)", 89, 100 * (1 - pt$n_kept[3] / pt$n_kept[1]), 0)
hwt <- (pt$ci_hi[2:3] - pt$ci_lo[2:3]) / 2
chk("LaLonde trimmed CI half-width, smaller (1,700)", 1700, round(min(hwt), -2), 0)
chk("LaLonde trimmed CI half-width, larger (2,200)", 2200, round(max(hwt), -2), 0)
chk("LaLonde trimmed CI half-widths, to the nearest 1,000 (roughly 2,000)", 2000, unique(round(hwt, -3)), 0)

## ---- NHEFS: the confidence half-width and the software spread: Table 4, Steps 6 and 7 ----
hw6 <- c(with(pf("AIPW (one-step)"), (ci_hi - ci_lo) / 2), with(pf("TMLE"), (ci_hi - ci_lo) / 2),
         sapply(SPL, function(e) with(ps(e, 100), (ci_hi - ci_lo) / 2)))
chk("NHEFS doubly-robust CI half-width, mean of six (0.9)", 0.9, mean(hw6), 1)
chk("NHEFS software spread / that half-width (one fifth)", 0.2, diff(range(swv)) / mean(hw6), 1)
## the cited Hernan-Robins benchmark of about 3.4-3.5 kg: every doubly-robust interval of Figure 5 covers it, and
## every software estimate lies within that half-width of it (Step 7)
dri <- rbind(pf("AIPW (one-step)"), pf("TMLE"), p[p$estimator %in% SPL & !is.na(p$r) & p$r %in% c(1, 100), ])
chk("NHEFS every doubly-robust interval covers 3.4-3.5 kg (Figure 5)", 1, as.numeric(all(dri$ci_lo <= 3.4 & dri$ci_hi >= 3.5)), 0)
chk("NHEFS software estimates within the half-width of 3.4-3.5 kg", 1,
    as.numeric(all(swv >= 3.4 - mean(hw6) & swv <= 3.5 + mean(hw6))), 0)
## Table 11's hand-coded DML, tmle and AIPW entries are the first of the 20 seeded runs (Figure 10)
r1 <- mr[mr$rep == 1, ]
chk("Table 11 hand-coded DML = seed 1 of Figure 10", 0, abs(h$DML_crossfit - r1$hand_DML), 6)
chk("Table 11 tmle = seed 1 of Figure 10", 0, abs(swr("^tmle \\(")$TMLE - r1$tmle_pkg), 6)
chk("Table 11 AIPW = seed 1 of Figure 10", 0, abs(swr("Zhong")$AIPW - r1$AIPW_pkg), 6)

## ---- IHDP: Figure 7's mean CI widths and design, Box 13's three-decimal values ----
chk("IHDP panels per estimator (Figure 7 = 1000)", 1000, unique(g$n_panels), 0)
for (w in list(c("AIPW", 1.04), c("TMLE", 1.04), c("AIPW = DML", 1.34), c("TMLE (CV-TMLE)", 1.33),
               c("AIPW (DC-AIPW)", 1.45), c("TMLE (DC-TMLE)", 1.44)))
  chk(sprintf("IHDP mean CI width, %s (Figure 7 = %s)", w[1], w[2]), as.numeric(w[2]), gv(w[1], "mean_CI_width"))
chk("IHDP CV-TMLE coverage MCSE (Box 13 = 0.004)", 0.004, gv("TMLE (CV-TMLE)", "cov_mcse"), 3)
chk("IHDP DML coverage (Box 13 = 0.991)", 0.991, gv("AIPW = DML", "coverage"), 3)
chk("IHDP CV-TMLE coverage (Box 13 = 0.994)", 0.994, gv("TMLE (CV-TMLE)", "coverage"), 3)

## ---- Figure 9's density inputs ----
## md5 after removing carriage returns, so the check holds with either line ending
md5_lf <- function(f) { b <- readBin(f, "raw", file.size(f)); tf <- tempfile(); writeBin(b[b != as.raw(13)], tf)
  unname(tools::md5sum(tf)) }
chk("Figure 9 input positivity_overlap.csv (md5 cad194b0)", 1,
    as.numeric(md5_lf(file.path(OUT, "positivity_overlap.csv")) == "cad194b055af08cb8f5b995cd69f44c5"), 0)
chk("Figure 9 input nhefs_overlap_data.csv (md5 04e6240c)", 1,
    as.numeric(md5_lf(file.path(OUT, "nhefs_overlap_data.csv")) == "04e6240c0303b2a5d17c027fac104596"), 0)
po <- rd("positivity_overlap.csv"); pl <- po[grepl("^LaLonde", po$dataset), ]
chk("Figure 9 LaLonde treated scores (185)", 185, sum(pl$group == "treated"), 0)
chk("Figure 9 LaLonde comparison scores (2490)", 2490, sum(pl$group == "control"), 0)
chk("Figure 9 NHEFS quit scores (403)", 403, sum(ov$group == "quit smoking"), 0)
chk("Figure 9 NHEFS did-not-quit scores (1163)", 1163, sum(ov$group == "did not quit"), 0)
## the clipped peak: R's default bandwidth on a fine grid (the default 512 points give about 123)
chk("Figure 9 LaLonde comparison peak density (154)", 154,
    max(density(pl$pi[pl$group == "control"], n = 2^14)$y), 0)

## ---- report ----
res <- do.call(rbind, .rows)
cat(sprintf("%-58s %9s %9s  %s\n", "check", "printed", "computed", "status"))
for (i in seq_len(nrow(res)))
  cat(sprintf("%-58s %9s %9s  %s\n", substr(res$check[i], 1, 58), res$printed[i],
              ifelse(is.na(res$computed[i]), "NA", res$computed[i]), res$status[i]))
nfail <- sum(res$status == "FAIL")
cat(sprintf("\n%d checks: %d PASS, %d FAIL\n", nrow(res), sum(res$status == "PASS"), nfail))
if (nfail) { cat("FAILING:\n"); print(res[res$status == "FAIL", c("check", "printed", "computed")], row.names = FALSE)
  if (!interactive() && sys.nframe() == 0L) quit(status = 1) }
