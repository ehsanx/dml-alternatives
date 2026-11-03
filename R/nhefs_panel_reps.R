## =============================================================================
## nhefs_panel_reps.R -- the FULL NHEFS estimator panel (Table 6, Figure 5) with repeated cross-fitting,
## and the single-split estimates behind the split-to-split SDs (Step 4, Table 10).
##
## NHEFS is ONE real dataset (n=1566, weight change in kg, NO ground truth). This
## script produces, on ONE shared 5-learner Super Learner library and one truncation
## rule, EVERY estimator in the panel:
##   * FULL-SAMPLE (r-independent, computed once, pinned seed):
##       g-computation, IPTW, AIPW (one-step), TMLE
##   * SPLIT-BASED, at every r = 1..R_MAX (median-aggregated a la Karim & Mondol 2025):
##       DML (cross-fit AIPW), CV-TMLE, DC-AIPW, DC-TMLE
##   -> so r=1 (single arbitrary split) AND r=100 (repeated-cross-fit median) are
##      both available for the split rows, as in the RHC panel (Figure 6).
##
## REPRODUCIBILITY: single-thread BLAS/OpenMP (SL.ranger is already num.threads=1)
## and a fixed seed (DML_SEED) before the full-sample fit pin nuisance_full's internal SL CV
## folds; the cross-fit/double folds of repetition j are seeded via make_folds(seed=100+j).
##
## Parallelised OVER THE REPETITIONS (independent). Resume-safe (part per rep).
## Run from the repository root (sources R/utils_data.R and R/estimators.R).
## Config via env vars: DML_ROOT DML_DATA_DIR DML_OUT R_MAX R_GRID NCORES DML_SEED FRESH.
## Writes  output/nhefs_reps_raw.csv       (per rep x split-method: estimate, se)
##         output/nhefs_full_panel.csv      (the 4 full-sample estimators)
##         output/nhefs_panel_table.csv     (COMBINED panel: full-sample + split at each r)
## The frozen nhefs_panel_table.csv and nhefs_reps_raw.csv are in docs/data/.
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
PARTS    <- file.path(OUT, "nhefs_reps_parts")
dir.create(OUT,   showWarnings = FALSE, recursive = TRUE)
dir.create(PARTS, showWarnings = FALSE, recursive = TRUE)
if (nzchar(Sys.getenv("FRESH"))) invisible(file.remove(Sys.glob(file.path(PARTS, "rep_*.csv"))))

lib   <- c("SL.mean", "SL.glm", "SL.glmnet", "SL.ranger", "SL.earth")   # FULL library
R_MAX <- as.integer(Sys.getenv("R_MAX", "100"))
SEED  <- as.integer(Sys.getenv("DML_SEED", "1"))
rg    <- Sys.getenv("R_GRID", "")
GRID  <- if (nzchar(rg)) as.integer(strsplit(rg, ",")[[1]]) else seq_len(R_MAX)
GRID  <- sort(unique(pmin(GRID, R_MAX)))
gridstr <- if (identical(as.integer(GRID), seq_len(R_MAX))) sprintf("1..%d (every r)", R_MAX) else paste(GRID, collapse = ",")
NCORE <- as.integer(Sys.getenv("NCORES", as.character(max(1L, detectCores() - 1L))))
NCORE <- max(1L, min(NCORE, R_MAX))

d0 <- tryCatch(load_nhefs(), error = function(e)
  stop(sprintf("could not load NHEFS : %s", conditionMessage(e))))
yr0 <- range(d0$Y)
say(sprintf("NHEFS panel: n=%d (%d quit smoking), R_MAX=%d, report r={%s}, lib={%s}, cores=%d, seed=%d",
            length(d0$Y), sum(d0$A), R_MAX, gridstr, paste(lib, collapse = ", "), NCORE, SEED))
say(sprintf("  DML_ROOT = %s", ROOT))
say(sprintf("  engine   = %s , %s", LOADER, ESTIM))
say(sprintf("  OUT      = %s   (parts -> %s)", OUT, PARTS))

## ---- 1. FULL-SAMPLE panel (r-independent; computed once with a pinned seed) --
say("computing full-sample panel (g-computation, IPTW, AIPW, TMLE) ...")
set.seed(SEED)
nuF <- nuisance_full(d0$Y, d0$A, d0$X, lib, d0$family_y)
mkfull <- function(res, label) data.frame(
  estimator = label, nuisance = "full sample", r = NA_integer_,
  estimate = res$estimate, se = res$se, ci_lo = res$ci[1], ci_hi = res$ci[2],
  ci_width = res$ci[2] - res$ci[1], partition_sd = NA_real_, row.names = NULL)
full_panel <- rbind(
  mkfull(est_gcomp(d0$Y, d0$A, nuF),         "g-computation"),
  mkfull(est_iptw (d0$Y, d0$A, nuF),         "IPTW"),
  mkfull(est_aipw (d0$Y, d0$A, nuF),         "AIPW (one-step)"),
  mkfull(est_tmle (d0$Y, d0$A, nuF, yr0),    "TMLE"))
write.csv(full_panel, file.path(OUT, "nhefs_full_panel.csv"), row.names = FALSE)
say("  full-sample panel done.")

## ---- 2. SPLIT-BASED estimators: repeat the split r=1..R_MAX -----------------
M <- list(`DML (cross-fit AIPW)` = c("crossfit", "aipw"),
          `CV-TMLE`              = c("crossfit", "tmle"),
          `DC-AIPW`              = c("double",   "aipw"),
          `DC-TMLE`              = c("double",   "tmle"))
NUIS <- c(`DML (cross-fit AIPW)` = "cross-fit", `CV-TMLE` = "cross-fit",
          `DC-AIPW` = "double cross-fit", `DC-TMLE` = "double cross-fit")

combine <- function(th, vr, rr) {                # median-aggregate first rr reps
  th <- th[seq_len(rr)]; vr <- vr[seq_len(rr)]
  p <- stats::median(th); v <- stats::median(vr + (th - p)^2)
  c(est = p, se = sqrt(v))
}

one_rep <- function(j) {
  pf <- file.path(PARTS, sprintf("rep_%03d.csv", j))
  if (file.exists(pf)) return(read.csv(pf))       # resume
  t0 <- Sys.time()
  d  <- load_nhefs(); yr <- range(d$Y)
  nuX <- nuisance_crossfit(d$Y, d$A, d$X, 5L, lib, d$family_y, seed = 100 + j)
  nuD <- nuisance_double  (d$Y, d$A, d$X, 3L, lib, d$family_y, seed = 100 + j)
  rows <- list()
  for (nm in names(M)) {
    nu <- if (M[[nm]][1] == "crossfit") nuX else nuD
    e  <- if (M[[nm]][2] == "aipw") est_aipw(d$Y, d$A, nu) else est_tmle(d$Y, d$A, nu, yr)
    rows[[length(rows) + 1]] <- data.frame(rep = j, method = nm,
                                           estimate = e$estimate, se = e$se, row.names = NULL)
  }
  df <- do.call(rbind, rows)
  write.csv(df, pf, row.names = FALSE)
  cat(sprintf("  [ok]   rep %d/%d  (%.0fs)\n", j, R_MAX,
              as.numeric(difftime(Sys.time(), t0, units = "secs")))); flush.console()
  df
}

t0 <- Sys.time()
cl <- makeCluster(NCORE, outfile = "")
on.exit(try(stopCluster(cl), silent = TRUE), add = TRUE)
clusterExport(cl, c("ROOT", "LOADER", "ESTIM", "DDIR", "DATA_DIR", "PARTS", "lib", "R_MAX", "M", "one_rep"))
initok <- clusterEvalQ(cl, {
  suppressMessages({library(SuperLearner); library(ranger); library(glmnet); library(earth)})
  Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1")
  source(LOADER); source(ESTIM)
  DATA_DIR <- Sys.getenv("DML_DATA_DIR", DDIR)
  exists("load_nhefs") && exists("nuisance_crossfit") && exists("est_tmle")
})
if (!all(unlist(initok))) stop("workers failed to load the engine -- check package availability on the node")
say(sprintf("cluster of %d workers initialised; running %d repetitions on NHEFS ...", NCORE, R_MAX))

res_list <- parLapply(cl, seq_len(R_MAX), function(j) {
  tryCatch(one_rep(j),
           error = function(e) {
             cat(sprintf("  [FAIL] rep %d: %s\n", j, conditionMessage(e))); flush.console()
             data.frame(rep = j, error = conditionMessage(e))
           })
})
stopCluster(cl)

parts  <- Sys.glob(file.path(PARTS, "rep_*.csv"))
failed <- Filter(function(x) is.data.frame(x) && "error" %in% names(x), res_list)
if (length(parts) == 0L) {
  say("\nNO repetition produced output. Failures:")
  for (f in failed) say(sprintf("  rep %d: %s", f$rep[1], f$error[1]))
  stop("all repetitions failed -- see messages above; no CSV written")
}
raw <- do.call(rbind, lapply(parts, read.csv))
write.csv(raw, file.path(OUT, "nhefs_reps_raw.csv"), row.names = FALSE)

## ---- 3. aggregate split methods per r, then COMBINE with the full panel -----
split_rows <- list()
for (nm in names(M)) {
  s  <- raw[raw$method == nm, ]; s <- s[order(s$rep), ]
  th <- s$estimate; vr <- s$se^2
  for (rr in GRID) {
    cc <- combine(th, vr, rr); est <- unname(cc["est"]); se <- unname(cc["se"])
    split_rows[[length(split_rows) + 1]] <- data.frame(
      estimator = nm, nuisance = unname(NUIS[nm]), r = rr, estimate = est, se = se,
      ci_lo = est - 1.96 * se, ci_hi = est + 1.96 * se, ci_width = 2 * 1.96 * se,
      partition_sd = stats::sd(th), row.names = NULL)
  }
}
panel <- rbind(full_panel, do.call(rbind, split_rows))
write.csv(panel, file.path(OUT, "nhefs_panel_table.csv"), row.names = FALSE)

## ---- 4. console summary -----------------------------------------------------
say(sprintf("\nDONE in %.1f min; reps ok: %d/%d%s",
            as.numeric(difftime(Sys.time(), t0, units = "mins")), length(parts), R_MAX,
            if (length(failed)) sprintf("  (FAILED: %s)",
              paste(vapply(failed, function(f) as.character(f$rep[1]), ""), collapse = ",")) else ""))
say(sprintf("wrote %s , %s , %s",
            file.path(OUT, "nhefs_full_panel.csv"), file.path(OUT, "nhefs_reps_raw.csv"),
            file.path(OUT, "nhefs_panel_table.csv")))
say("\nNHEFS estimator panel (weight change, kg):")
say("  full sample:")
print(full_panel[, c("estimator", "estimate", "se", "ci_lo", "ci_hi")], row.names = FALSE, digits = 4)
kr <- intersect(c(1, 25, 50, 100), GRID); if (!length(kr)) kr <- sort(unique(c(min(GRID), max(GRID))))
say(sprintf("  split-based at r in {%s}  [partition_sd = spread of single-split estimates]:", paste(kr, collapse = ",")))
sp <- do.call(rbind, split_rows)
print(sp[sp$r %in% kr, c("estimator", "r", "estimate", "se", "ci_lo", "ci_hi", "partition_sd")],
      row.names = FALSE, digits = 4)
