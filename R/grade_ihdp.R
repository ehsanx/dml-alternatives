## =============================================================================
## grade_ihdp.R -- SEED-AVERAGED IHDP grading (Figure 7, Step 5, Box 13). Single-split coverage depends on
## the arbitrary fold partition, so we run the FULL 8-estimator panel over each of N_REAL known-truth
## realizations x N_SEED fold seeds, and average bias / RMSE / coverage / CI width over BOTH realizations and
## seeds, with a Monte Carlo SE that reflects the realization spread. This reports the honest single-split
## coverage AVERAGED over partitions (not median-aggregation, which over-covers).
##
## Run from the repository root (sources R/utils_data.R and R/estimators.R), with the IHDP files
## under $DML_DATA_DIR/ihdp/. utils_data.R imports `haven` at source time.
## Parallelised over (realization x seed) PAIRS. Resume-safe (part per pair).
## Config env vars: DML_ROOT DML_DATA_DIR DML_OUT N_REAL N_SEED NCORES FRESH.
## Writes output/grade_raw.csv (per realization x seed x estimator) and output/grade_table.csv
## (per estimator: mean bias, RMSE, coverage, CI width, MCSEs); frozen copies in docs/data/.
suppressMessages({
  library(SuperLearner); library(ranger); library(glmnet); library(earth); library(parallel)
})
Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1")
say <- function(...) { cat(...); cat("\n"); flush.console() }

ROOT   <- Sys.getenv("DML_ROOT", getwd())          # run from the repository root
LOADER <- file.path(ROOT, "R", "utils_data.R"); ESTIM <- file.path(ROOT, "R", "estimators.R")
DDIR   <- Sys.getenv("DML_DATA_DIR", file.path(ROOT, "data")); ODIR <- file.path(ROOT, "output")
if (!file.exists(LOADER) || !file.exists(ESTIM))
  stop(sprintf("R/utils_data.R and R/estimators.R not found under DML_ROOT=%s; run from the repository root", ROOT))
source(LOADER); source(ESTIM)
DATA_DIR <- Sys.getenv("DML_DATA_DIR", DDIR)
OUT      <- Sys.getenv("DML_OUT", ODIR)
PARTS    <- file.path(OUT, "grade_parts")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
dir.create(PARTS, showWarnings = FALSE, recursive = TRUE)
if (nzchar(Sys.getenv("FRESH"))) invisible(file.remove(Sys.glob(file.path(PARTS, "pair_*.csv"))))

lib     <- default_SL                                  # full 5-learner library from estimators.R
navail  <- length(Sys.glob(file.path(DATA_DIR, "ihdp", "ihdp_npci_*.csv")))
if (navail == 0L) stop("no IHDP files under ", file.path(DATA_DIR, "ihdp"))
N_REAL  <- min(as.integer(Sys.getenv("N_REAL", "50")), navail)
N_SEED  <- as.integer(Sys.getenv("N_SEED", "20"))
NCORE   <- as.integer(Sys.getenv("NCORES", as.character(max(1L, detectCores() - 1L))))

pairs <- expand.grid(i = seq_len(N_REAL), s = seq_len(N_SEED))
NCORE <- max(1L, min(NCORE, nrow(pairs)))
say(sprintf("seed-averaged IHDP grade: %d realizations x %d seeds = %d panels, lib={%s}, cores=%d",
            N_REAL, N_SEED, nrow(pairs), paste(lib, collapse = ", "), NCORE))
say(sprintf("  DML_ROOT=%s ; engine=%s", ROOT, ESTIM))
say(sprintf("  DATA_DIR=%s (IHDP files: %d) ; OUT=%s", DATA_DIR, navail, OUT))

one_pair <- function(k) {
  i <- pairs$i[k]; s <- pairs$s[k]
  pf <- file.path(PARTS, sprintf("pair_%03d_%02d.csv", i, s))
  d <- load_ihdp(i)
  p <- run_panel(d$Y, d$A, d$X, family_y = d$family_y, SL.library = lib,
                 seed = 100L + s, truth = d$truth)
  p$realization <- i; p$seed <- s; p$truth <- d$truth
  p$covers <- p$truth >= p$ci_lo & p$truth <= p$ci_hi          # NA when ci is NA (g-comp/IPTW)
  p$ci_width <- p$ci_hi - p$ci_lo
  write.csv(p, pf, row.names = FALSE)
  cat(sprintf("  [ok] realization %d seed %d\n", i, s)); flush.console()
  p
}

## resume: skip pairs whose part file exists
have <- basename(Sys.glob(file.path(PARTS, "pair_*.csv")))
todo <- which(!sprintf("pair_%03d_%02d.csv", pairs$i, pairs$s) %in% have)
if (length(todo) < nrow(pairs)) say(sprintf("resume: %d/%d panels already done; %d remaining",
                                             nrow(pairs) - length(todo), nrow(pairs), length(todo)))

t0 <- Sys.time()
if (length(todo)) {
  NCORE <- max(1L, min(NCORE, length(todo)))
  cl <- makeCluster(NCORE, outfile = "")
  on.exit(try(stopCluster(cl), silent = TRUE), add = TRUE)
  clusterExport(cl, c("ROOT", "LOADER", "ESTIM", "DDIR", "DATA_DIR", "PARTS", "lib", "pairs", "one_pair"))
  ok <- clusterEvalQ(cl, {
    suppressMessages({library(SuperLearner); library(ranger); library(glmnet); library(earth)})
    Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1")
    source(LOADER); source(ESTIM)
    DATA_DIR <- Sys.getenv("DML_DATA_DIR", DDIR)
    exists("load_ihdp") && exists("run_panel")
  })
  if (!all(unlist(ok))) stop("workers failed to load engine (is `haven` installed? install.packages('haven'))")
  say(sprintf("cluster of %d workers; running %d panels ...", NCORE, length(todo)))
  invisible(parLapply(cl, todo, function(k) tryCatch(one_pair(k),
    error = function(e) { cat(sprintf("  [FAIL] pair %d: %s\n", k, conditionMessage(e))); flush.console(); NULL })))
  stopCluster(cl)
}

parts <- Sys.glob(file.path(PARTS, "pair_*.csv"))
if (!length(parts)) stop("no panels produced; see [FAIL] messages above")
raw <- do.call(rbind, lapply(parts, read.csv, stringsAsFactors = FALSE))
write.csv(raw, file.path(OUT, "grade_raw.csv"), row.names = FALSE)

## aggregate per estimator; MCSE from the per-realization (seed-averaged) values
est_order <- c("g-computation","IPTW","AIPW","TMLE","AIPW = DML","TMLE (CV-TMLE)","AIPW (DC-AIPW)","TMLE (DC-TMLE)")
agg <- do.call(rbind, lapply(est_order, function(e) {
  s <- raw[raw$estimator == e, ]; if (!nrow(s)) return(NULL)
  per_real_cov  <- tapply(s$covers, s$realization, function(x) mean(x, na.rm = TRUE))
  per_real_bias <- tapply(s$bias,   s$realization, mean)
  data.frame(estimator = e, nuisance = s$nuisance[1],
             mean_bias = mean(s$bias), bias_mcse = stats::sd(per_real_bias) / sqrt(length(per_real_bias)),
             rmse = sqrt(mean(s$bias^2)),
             coverage = mean(s$covers, na.rm = TRUE),
             cov_mcse = if (all(is.na(s$covers))) NA_real_ else stats::sd(per_real_cov) / sqrt(sum(!is.na(per_real_cov))),
             mean_CI_width = mean(s$ci_width, na.rm = TRUE),
             n_panels = nrow(s), row.names = NULL)
}))
write.csv(agg, file.path(OUT, "grade_table.csv"), row.names = FALSE)

say(sprintf("\nDONE in %.1f min; %d panels (%d realizations x %d seeds)",
            as.numeric(difftime(Sys.time(), t0, units = "mins")), nrow(raw) / length(est_order), N_REAL, N_SEED))
say("seed-averaged IHDP grade (coverage MCSE over realizations):")
print(agg[, c("estimator","mean_bias","rmse","coverage","cov_mcse","mean_CI_width")], row.names = FALSE, digits = 3)
