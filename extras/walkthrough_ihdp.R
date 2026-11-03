## walkthrough_ihdp.R -- supplementary analysis (not a manuscript exhibit): a step-by-step walkthrough of
## the estimator panel on one dataset (IHDP realization 1). Run from the repository root:
##       Rscript extras/walkthrough_ihdp.R          (full Super Learner library)
##       SL_FAST=1 Rscript extras/walkthrough_ihdp.R (fast glm-only, for debugging)
## =============================================================================
source("R/utils_data.R")
source("R/estimators.R")

dat <- load_ihdp(1)
cat("Dataset :", dat$name, "\n")
cat("n =", length(dat$Y), " | treated =", sum(dat$A), " | control =", sum(dat$A == 0), "\n")
cat("True (sample) ATE        :", round(dat$truth, 4), "\n")
cat("Naive difference in means:", round(mean(dat$Y[dat$A == 1]) - mean(dat$Y[dat$A == 0]), 4),
    "  <- confounded\n\n")

lib <- if (nzchar(Sys.getenv("SL_FAST"))) c("SL.mean", "SL.glm") else default_SL
cat("Super Learner library:", paste(lib, collapse = ", "), "\n\n")

res <- run_panel(dat$Y, dat$A, dat$X, family_y = dat$family_y,
                 SL.library = lib, truth = dat$truth, seed = 1L)

res_print <- transform(res,
  estimate = round(estimate, 3), se = round(se, 3),
  ci_lo = round(ci_lo, 3), ci_hi = round(ci_hi, 3), bias = round(bias, 3))
print(res_print, row.names = FALSE)

## ---- cross-check hand-coded TMLE vs the tmle package (same SL library) ------
if (requireNamespace("tmle", quietly = TRUE)) {
  suppressMessages(library(tmle))
  set.seed(1)
  tm <- tmle(Y = dat$Y, A = dat$A, W = as.data.frame(dat$X),
             Q.SL.library = lib, g.SL.library = lib)
  hand <- res$estimate[res$estimator == "TMLE" & res$nuisance == "full"]
  cat(sprintf("\nCross-check  hand TMLE (full) = %.4f   vs   tmle package = %.4f\n",
              hand, tm$estimates$ATE$psi))
}
