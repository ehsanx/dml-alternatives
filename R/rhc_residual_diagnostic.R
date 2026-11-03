## =============================================================================
## rhc_residual_diagnostic.R -- the RHC residual diagnostic of Step 4: with in-sample nuisances
## the outcome model's residuals are smaller than out of fold, and so is AIPW's mean weighted
## correction. Refits the full-sample and five-fold cross-fit nuisances exactly as
## real_application.R does (set.seed before the full-sample fit; cross-fit folds seeded with 1)
## and reports, for each fitting scheme:
##   resid_rms   root-mean-square of Y - Q(A, X)
##   gcomp       mean of Q(1, X) - Q(0, X)               (the plug-in)
##   correction  mean of H(A, X) {Y - Q(A, X)}            (AIPW's one-step correction)
##   aipw        gcomp + correction                      (equals the AIPW / DML estimate)
## Writes output/rhc_residual_diagnostic.csv. Needs the RHC data (DML_DATA_DIR, see utils_data.R).
## Run from the repository root:  Rscript R/rhc_residual_diagnostic.R
source("R/utils_data.R")
source("R/estimators.R")
Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1")
SEED <- as.integer(Sys.getenv("DML_SEED", "1"))
OUT <- "output"; dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
lib <- default_SL

d <- load_rhc(outcome = "dth30")
set.seed(SEED)
nu_F <- nuisance_full(d$Y, d$A, d$X, lib, d$family_y)
nu_X <- nuisance_crossfit(d$Y, d$A, d$X, 5L, lib, d$family_y, seed = 1L)

diag_one <- function(nu, scheme) {
  H <- d$A / nu$pi - (1 - d$A) / (1 - nu$pi)
  gc <- mean(nu$Q1 - nu$Q0); co <- mean(H * (d$Y - nu$QA))
  data.frame(scheme = scheme, resid_rms = sqrt(mean((d$Y - nu$QA)^2)),
             gcomp = gc, correction = co, aipw = gc + co)
}
res <- rbind(diag_one(nu_F, "full sample (in-sample)"), diag_one(nu_X, "cross-fit V=5 (out of fold)"))
write.csv(res, file.path(OUT, "rhc_residual_diagnostic.csv"), row.names = FALSE)
print(res, row.names = FALSE, digits = 4)
