## =============================================================================
## real_application.R -- the panel on REAL observational epi data with NO ground
## truth (RHC, NHEFS). Here we cannot grade bias/coverage; instead we report the
## estimates + CIs and how well the estimators AGREE (a practical reassurance check
## on real data). Writes output/real_application.{csv,md} (RHC rows: Figure 6, Step 4).
source("R/utils_data.R")
source("R/estimators.R")
## reproducibility: single-thread the BLAS/OpenMP backends (SL.ranger is already
## num.threads=1) and set a global seed before each panel so nuisance_full's internal
## SuperLearner CV folds are pinned too (the cross-fit/double folds are seeded via
## make_folds). Without this the full-sample estimates drift run-to-run.
Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1")
SEED <- as.integer(Sys.getenv("DML_SEED", "1"))
OUT <- "output"; dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
lib <- if (nzchar(Sys.getenv("SL_FAST"))) c("SL.mean", "SL.glm") else default_SL

apps <- list(load_rhc(outcome = "dth30"), load_nhefs())
all <- list()
for (d in apps) {
  cat("\n====", d$name, " (n =", length(d$Y), ",", d$outcome_type, ") ====\n")
  t0 <- Sys.time()
  set.seed(SEED)
  r <- run_panel(d$Y, d$A, d$X, family_y = d$family_y, SL.library = lib, truth = NA, seed = 1L)
  r$dataset <- d$name
  dr <- r[grepl("AIPW|TMLE", r$estimator), ]
  cat(sprintf("   elapsed %.0fs | DR estimates span [%.4f, %.4f] (agreement spread = %.4f)\n",
              as.numeric(difftime(Sys.time(), t0, units = "secs")),
              min(dr$estimate), max(dr$estimate), diff(range(dr$estimate))))
  all[[length(all) + 1]] <- r
}
res <- do.call(rbind, all)
res <- res[, c("dataset", "estimator", "nuisance", "estimate", "se", "ci_lo", "ci_hi")]
write.csv(res, file.path(OUT, "real_application.csv"), row.names = FALSE)

fmt <- function(x, dg = 3) ifelse(is.na(x), "—", formatC(x, format = "f", digits = dg))
md <- c("# Real-application datasets (RHC, NHEFS) — NO ground truth",
        "", "These are observational: we report estimates, SEs and 95% CIs and check that the estimators",
        "AGREE (they cannot be graded for bias/coverage). RHC outcome = 30-day mortality (risk difference);",
        "NHEFS outcome = weight change (kg). Identical Super Learner library and folds throughout.", "")
for (dn in unique(res$dataset)) {
  sub <- res[res$dataset == dn, ]
  md <- c(md, sprintf("## %s", dn), "",
          "| estimator | nuisance | estimate | SE | 95% CI |",
          "|---|---|---:|---:|:--:|")
  for (i in seq_len(nrow(sub))) md <- c(md, sprintf(
    "| %s | %s | %s | %s | (%s, %s) |", sub$estimator[i], sub$nuisance[i], fmt(sub$estimate[i]),
    fmt(sub$se[i]), fmt(sub$ci_lo[i]), fmt(sub$ci_hi[i])))
  md <- c(md, "")
}
writeLines(md, file.path(OUT, "real_application.md"))
cat("\n"); print(res, row.names = FALSE, digits = 4); cat("\nDONE real_application\n")
