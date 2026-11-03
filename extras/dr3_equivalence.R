## dr3_equivalence.R -- supplementary check (not a manuscript exhibit): does the dr3 package compute what this
## repository's engine, R/estimators.R, computes? Both run side by side on NHEFS with the same data, Super Learner
## library, folds, propensity truncation and random-number initialization, and their estimates are compared:
##   full-sample AIPW and TMLE; seed-1 cross-fit DML and CV-TMLE; seed-1 DC-AIPW and DC-TMLE.
##   Rscript extras/dr3_equivalence.R              36 Super Learner fits, a few minutes
##   Rscript extras/dr3_equivalence.R --seeds 20   also the cross-fit DML for seeds 2 to 20 (Figure 10's seeds),
##                                                 about half an hour
## Prints each pair of estimates, their absolute difference and PASS or FAIL against TOL, and exits with status 1
## if any comparison fails. NHEFS comes from causaldata, so no raw data are needed. RHC is not used: the engine
## reads Vanderbilt's CSV export of that dataset and dr3 its R-format export, which differ slightly, so a
## difference there would not be an implementation difference.
## Needs dr3 (remotes::install_github("ehsanx/dr3")). Run from the repository root.
## =============================================================================
TOL <- 1e-10
args <- commandArgs(trailingOnly = TRUE)
nseeds <- if (length(i <- which(args == "--seeds"))) suppressWarnings(as.integer(args[i[1] + 1])) else 1L
if (is.na(nseeds) || nseeds < 1L) stop("--seeds needs a positive integer")

## as in the analysis scripts: one BLAS/OpenMP thread
Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1")
source("R/utils_data.R")
source("R/estimators.R")
suppressMessages(library(dr3))
cat("dr3", as.character(packageVersion("dr3")), "| SuperLearner", as.character(packageVersion("SuperLearner")),
    "|", R.version.string, "\n")

## ---- the same inputs: data, library, folds, truncation -----------------------
d   <- load_nhefs()
d3  <- dr3::dr_nhefs()
lib <- default_SL
gb  <- c(0.025, 0.975)
n   <- length(d$Y)
inputs <- c(
  data    = identical(d$Y, d3$Y) && identical(d$A, d3$A) && identical(d$X, d3$X),
  library = identical(lib, dr3::default_dr_library),
  folds   = identical(make_folds(n, 5L, 1L), dr3::dr_folds(n, 5L, 1L)) &&
            identical(make_folds(n, 3L, 1L), dr3::dr_folds(n, 3L, 1L)))
cat(sprintf("identical %s: %s\n", names(inputs), ifelse(inputs, "yes", "NO")), sep = "")
if (!all(inputs)) { cat("FAIL: the inputs differ, so the estimates cannot be compared\n"); quit(status = 1) }

## ---- fits: the same random-number initialization for both ----------------------
## The full-sample scheme takes no seed, so set.seed() comes first, as in R/nhefs_analysis.R; the cross-fit
## schemes seed their own fold construction, which also fixes the random-number state of the fits after it.
fit_pair <- function(scheme, seed) switch(scheme,
  full = list(
    engine = { set.seed(seed); nuisance_full(d$Y, d$A, d$X, lib, d$family_y, gb) },
    dr3    = { set.seed(seed); dr3::dr_nuisance(d$Y, d$A, d$X, "full", SL.library = lib,
                                                family_y = d$family_y, gbound = gb) }),
  crossfit = list(
    engine = nuisance_crossfit(d$Y, d$A, d$X, 5L, lib, d$family_y, gb, seed = seed),
    dr3    = dr3::dr_nuisance(d$Y, d$A, d$X, "crossfit", V = 5L, SL.library = lib,
                              family_y = d$family_y, gbound = gb, seed = seed)),
  double = list(
    engine = nuisance_double(d$Y, d$A, d$X, 3L, lib, d$family_y, gb, seed = seed),
    dr3    = dr3::dr_nuisance(d$Y, d$A, d$X, "double", Vd = 3L, SL.library = lib,
                              family_y = d$family_y, gbound = gb, seed = seed)))

yr <- range(d$Y)
rows <- list(); nuis <- 0
add <- function(estimator, engine, dr3)
  rows[[length(rows) + 1]] <<- data.frame(estimator, engine, dr3, stringsAsFactors = FALSE)
nuisance_gap <- function(f) max(abs(unlist(f$engine[c("pi", "Q1", "Q0")]) - unlist(f$dr3[c("pi", "Q1", "Q0")])))
run_scheme <- function(scheme, seed, labels, tmle = TRUE) {
  t0 <- Sys.time(); f <- fit_pair(scheme, seed)
  nuis <<- max(nuis, nuisance_gap(f))
  add(labels[1], est_aipw(d$Y, d$A, f$engine)$estimate, dr3::dr_aipw(d$Y, d$A, f$dr3)$estimate)
  if (tmle) add(labels[2], est_tmle(d$Y, d$A, f$engine, yr)$estimate, dr3::dr_tmle(d$Y, d$A, f$dr3, yr)$estimate)
  cat(sprintf("  %-8s seed %2d: %.1f min\n", scheme, seed, as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}
cat("fitting (engine, then dr3, for each scheme):\n")
run_scheme("full", 1L, c("AIPW, full sample", "TMLE, full sample"))
run_scheme("crossfit", 1L, c("DML, seed 1", "CV-TMLE, seed 1"))
run_scheme("double", 1L, c("DC-AIPW, seed 1", "DC-TMLE, seed 1"))
for (s in seq_len(nseeds)[-1]) run_scheme("crossfit", s, sprintf("DML, seed %d", s), tmle = FALSE)

## ---- the comparison --------------------------------------------------------------
res <- do.call(rbind, rows)
res$abs_diff <- abs(res$engine - res$dr3)
res$result <- ifelse(res$abs_diff <= TOL, "PASS", "FAIL")
out <- transform(res, engine = sprintf("%.12f", engine), dr3 = sprintf("%.12f", dr3),
                 abs_diff = sprintf("%.1e", abs_diff))
cat("\nNHEFS estimates (kg), engine vs dr3, tolerance", format(TOL), "\n")
print(out, row.names = FALSE, right = FALSE)
cat(sprintf("\nlargest absolute difference: estimates %.1e, nuisance predictions (pi, Q1, Q0) %.1e\n",
            max(res$abs_diff), nuis))

## for information, not PASS/FAIL: the engine against the archived values the manuscript reports. Exact agreement
## needs the package versions of renv.lock on the platform that produced them (see the README).
sw <- read.csv("docs/data/nhefs_software.csv", stringsAsFactors = FALSE); h <- sw[grepl("^hand", sw$source), ][1, ]
dc <- read.csv("results/dcfix1_213b1138cc/nhefs_r1.csv", stringsAsFactors = FALSE)
mr <- read.csv("docs/data/nhefs_matched_repeated.csv", stringsAsFactors = FALSE)
arch <- c("AIPW, full sample" = h$AIPW, "TMLE, full sample" = h$TMLE, "DML, seed 1" = h$DML_crossfit,
          "DC-AIPW, seed 1" = dc$estimate[dc$method == "DC-AIPW"], "DC-TMLE, seed 1" = dc$estimate[dc$method == "DC-TMLE"],
          setNames(mr$hand_DML, sprintf("DML, seed %d", mr$rep))[-1])
hit <- res$estimator %in% names(arch)
cat(sprintf("for information, engine vs the archived values (%d estimates): largest absolute difference %.1e\n",
            sum(hit), max(abs(res$engine[hit] - arch[res$estimator[hit]]))))

if (any(res$result == "FAIL")) {
  cat("FAIL:", sum(res$result == "FAIL"), "of", nrow(res), "comparisons exceed the tolerance\n")
  quit(status = 1)
}
cat("PASS: dr3 reproduces the engine in all", nrow(res), "comparisons\n")
