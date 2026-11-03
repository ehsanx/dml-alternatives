## =============================================================================
## comparison_all3.R  --  the estimator panel on IHDP, LaLonde/NSW-PSID (observational) and Twins.
## One Super Learner library for every estimator, with folds matched within each fitting scheme. IHDP and Twins
## estimates are graded against their known sample-average effects (bias; does the 95% CI cover it).
## LaLonde/NSW-PSID is not graded: the loader's truth, $1,794, is the experimental estimate for a different
## population (the NSW experimental sample) and is kept only as a reference.
## Writes: output/results_all3.csv  and  output/results_all3.md  (LaLonde/NSW-PSID rows: Table 8, Figure 9)
source("R/utils_data.R")
source("R/estimators.R")
## reproducibility: single-thread the BLAS/OpenMP backends and set a seed before each
## run_panel so nuisance_full's internal SuperLearner CV folds are pinned too (the cross-fit
## and double-cross-fit folds are already seeded via make_folds). Mirrors real_application.R /
## nhefs_analysis.R; without this the full-sample estimates drift run-to-run.
Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1")
SEED <- as.integer(Sys.getenv("DML_SEED", "1"))

OUT <- "output"; dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
lib <- if (nzchar(Sys.getenv("SL_FAST"))) c("SL.mean", "SL.glm") else default_SL

datasets <- list(
  load_ihdp(realization = 1L),
  load_lalonde(setting = "psid"),                 # observational, poor overlap; truth column = $1,794, an
                                                  # experimental estimate for a different population (reference only)
  load_twins(n_sub = 3000L)                        # subsampled for runtime
)

all_rows <- list()
for (dat in datasets) {
  cat("\n==== ", dat$name, " (n =", length(dat$Y),
      if (grepl("^LaLonde", dat$name)) ", experimental reference (not graded) =" else ", known effect =",
      round(dat$truth, 4), ",", dat$outcome_type, ") ====\n")
  t0 <- Sys.time()
  set.seed(SEED)
  r <- run_panel(dat$Y, dat$A, dat$X, family_y = dat$family_y,
                 SL.library = lib, truth = dat$truth, seed = 1L)
  r$dataset <- dat$name
  r$truth <- dat$truth
  r$covers <- ifelse(is.na(r$ci_lo), NA, r$truth >= r$ci_lo & r$truth <= r$ci_hi)
  cat("   elapsed:", round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1), "s\n")
  all_rows[[length(all_rows) + 1]] <- r
}
res <- do.call(rbind, all_rows)
res <- res[, c("dataset", "truth", "estimator", "nuisance",
               "estimate", "se", "ci_lo", "ci_hi", "bias", "covers")]

## The CSV keeps bias and covers for every row, as in the frozen copy; for LaLonde/NSW-PSID they are relative
## to the experimental reference and are not used.
write.csv(res, file.path(OUT, "results_all3.csv"), row.names = FALSE)

## ---- pretty markdown table --------------------------------------------------
fmt <- function(x, d = 3) formatC(x, format = "f", digits = d)
graded <- !grepl("^LaLonde", res$dataset)
md <- c("# DML vs AIPW vs TMLE — three-benchmark comparison",
        "",
        "One Super Learner library for every estimator, with folds matched within each fitting scheme.",
        "IHDP and Twins are graded against their known sample-average effects (`covers`: does the 95% CI contain it).",
        "LaLonde/NSW–PSID is not graded: the $1,794 experimental estimate concerns a different population.",
        "")
for (dn in unique(res$dataset)) {
  sub <- res[res$dataset == dn, ]
  g <- !grepl("^LaLonde", dn)
  md <- c(md, sprintf(if (g) "## %s  (known sample-average effect = %s)" else
                        "## %s  (experimental reference for a different population = %s; not graded)",
                      dn, fmt(sub$truth[1], 3)), "",
          "| estimator | nuisance | estimate | SE | 95% CI | bias | covers |",
          "|---|---|---:|---:|:--:|---:|:--:|")
  for (i in seq_len(nrow(sub))) with(sub[i, ],
    md <<- c(md, sprintf("| %s | %s | %s | %s | (%s, %s) | %s | %s |",
      estimator, nuisance, fmt(estimate),
      ifelse(is.na(se), "—", fmt(se)),
      ifelse(is.na(ci_lo), "—", fmt(ci_lo)), ifelse(is.na(ci_hi), "—", fmt(ci_hi)),
      if (g) fmt(bias) else "—", if (!g || is.na(covers)) "—" else if (covers) "yes" else "**no**")))
  md <- c(md, "")
}
writeLines(md, file.path(OUT, "results_all3.md"))

cat("\n================ RESULTS (all three datasets) ================\n")
shown <- res
shown$bias[!graded] <- NA; shown$covers[!graded] <- NA
print(transform(shown, truth = round(truth, 3), estimate = round(estimate, 3),
                se = round(se, 3), bias = round(bias, 3),
                ci_lo = round(ci_lo, 3), ci_hi = round(ci_hi, 3)), row.names = FALSE)
cat("LaLonde/NSW-PSID: truth is the $1,794 experimental reference for a different population; not graded.\n")
cat("\nWrote:", file.path(OUT, "results_all3.csv"), "and results_all3.md\n")
