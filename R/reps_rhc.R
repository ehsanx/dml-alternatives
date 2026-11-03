## =============================================================================
## reps_rhc.R -- REPETITION study on the RHC dataset (Karim & Mondol 2025; Figures 6 and 8). RHC is ONE
## real dataset (n=5735, 30-day mortality, NO ground truth), so there is no simulation and no coverage/bias
## grading -- we repeat the random split r=1..R_MAX times on the SAME data and median-aggregate, and report
## how each split-based estimator's point estimate and CI STABILISE as r grows (removing dependence on any
## single partition). Full 5-learner Super Learner library; the four split-based estimators only (DML,
## CV-TMLE, DC-AIPW, DC-TMLE) -- full-sample AIPW/TMLE and g-comp/IPTW have no split to repeat.
##
## Parallelised OVER THE REPETITIONS (they are independent) over NCORES workers; resumes from
## per-repetition part files. Run from the repository root (sources R/utils_data.R and R/estimators.R),
## with the RHC data under $DML_DATA_DIR/rhc/.
## Config via env vars: DML_ROOT DML_DATA_DIR DML_OUT R_MAX R_GRID NCORES  (see README).
## Writes output/rhc_reps_raw.csv (per rep x method) and output/rhc_reps_table.csv
## (per method x r: median-aggregated estimate, se, CI, CI width, partition SD); frozen copies in docs/data/.
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
PARTS    <- file.path(OUT, "rhc_reps_parts")
dir.create(OUT,   showWarnings = FALSE, recursive = TRUE)
dir.create(PARTS, showWarnings = FALSE, recursive = TRUE)
if (nzchar(Sys.getenv("FRESH"))) invisible(file.remove(Sys.glob(file.path(PARTS, "rep_*.csv"))))  # FRESH=1 forces a clean restart

lib   <- c("SL.mean", "SL.glm", "SL.glmnet", "SL.ranger", "SL.earth")   # FULL library
R_MAX <- as.integer(Sys.getenv("R_MAX", "100"))
rg    <- Sys.getenv("R_GRID", "")               # empty -> report every r from 1..R_MAX
GRID  <- if (nzchar(rg)) as.integer(strsplit(rg, ",")[[1]]) else seq_len(R_MAX)
GRID  <- sort(unique(pmin(GRID, R_MAX)))
gridstr <- if (identical(as.integer(GRID), seq_len(R_MAX))) sprintf("1..%d (every r)", R_MAX) else paste(GRID, collapse = ",")
NCORE <- as.integer(Sys.getenv("NCORES", as.character(max(1L, detectCores() - 1L))))
NCORE <- max(1L, min(NCORE, R_MAX))

## load once to report n and fail early if the data are missing
d0 <- tryCatch(load_rhc(), error = function(e)
  stop(sprintf("could not load RHC from %s : %s", file.path(DATA_DIR, "rhc", "rhc.csv"), conditionMessage(e))))
say(sprintf("RHC repetition study: n=%d (%d treated, %d deaths), R_MAX=%d, report r={%s}, lib={%s}, cores=%d",
            length(d0$Y), sum(d0$A), sum(d0$Y), R_MAX, gridstr, paste(lib, collapse = ", "), NCORE))
say(sprintf("  DML_ROOT = %s", ROOT))
say(sprintf("  engine   = %s , %s", LOADER, ESTIM))
say(sprintf("  DATA_DIR = %s   (rhc.csv %s)", DATA_DIR,
            if (file.exists(file.path(DATA_DIR, "rhc", "rhc.csv"))) "found" else "MISSING"))
say(sprintf("  OUT      = %s   (parts -> %s)", OUT, PARTS))

M <- list(`DML (cross-fit)`     = c("crossfit", "aipw"),
          `CV-TMLE (cross-fit)` = c("crossfit", "tmle"),
          `DC-AIPW`             = c("double",   "aipw"),
          `DC-TMLE`             = c("double",   "tmle"))

combine <- function(th, vr, rr) {                # median-aggregate first rr repetitions
  th <- th[seq_len(rr)]; vr <- vr[seq_len(rr)]
  p <- stats::median(th); v <- stats::median(vr + (th - p)^2)
  c(est = p, se = sqrt(v))
}

one_rep <- function(j) {                         # one random split of the RHC data
  t0 <- Sys.time()
  d  <- load_rhc(); yr <- range(d$Y)
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
  write.csv(df, file.path(PARTS, sprintf("rep_%03d.csv", j)), row.names = FALSE)
  cat(sprintf("  [ok]   rep %d/%d  (%.0fs)\n", j, R_MAX,
              as.numeric(difftime(Sys.time(), t0, units = "secs")))); flush.console()
  df
}

## resume: skip repetitions whose part file already exists, so a rerun picks up where it stopped
done <- as.integer(sub(".*rep_(\\d+)\\.csv$", "\\1", Sys.glob(file.path(PARTS, "rep_*.csv"))))
done <- done[!is.na(done) & done >= 1 & done <= R_MAX]
todo <- setdiff(seq_len(R_MAX), done)
if (length(done)) say(sprintf("resume: %d/%d repetitions already have a part file; %d remaining", length(done), R_MAX, length(todo)))

t0 <- Sys.time(); res_list <- list()
if (length(todo)) {
  NCORE <- max(1L, min(NCORE, length(todo)))
  cl <- makeCluster(NCORE, outfile = "")
  on.exit(try(stopCluster(cl), silent = TRUE), add = TRUE)
  clusterExport(cl, c("ROOT", "LOADER", "ESTIM", "DDIR", "DATA_DIR", "PARTS", "lib", "R_MAX", "M", "one_rep"))
  initok <- clusterEvalQ(cl, {
    suppressMessages({library(SuperLearner); library(ranger); library(glmnet); library(earth)})
    Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1")
    source(LOADER); source(ESTIM)
    DATA_DIR <- Sys.getenv("DML_DATA_DIR", DDIR)
    exists("load_rhc") && exists("nuisance_crossfit") && exists("est_tmle")
  })
  if (!all(unlist(initok))) stop("workers failed to load the engine -- check package availability on the node")
  say(sprintf("cluster of %d workers initialised; running %d repetitions on the RHC dataset ...", NCORE, length(todo)))
  res_list <- parLapply(cl, todo, function(j) {
    tryCatch(one_rep(j),
             error = function(e) {
               cat(sprintf("  [FAIL] rep %d: %s\n", j, conditionMessage(e))); flush.console()
               data.frame(rep = j, error = conditionMessage(e))
             })
  })
  stopCluster(cl)
} else say("all repetitions already have part files; aggregating existing parts.")

parts  <- Sys.glob(file.path(PARTS, "rep_*.csv"))
failed <- Filter(function(x) is.data.frame(x) && "error" %in% names(x), res_list)
if (length(parts) == 0L) {
  say("\nNO repetition produced output. Failures:")
  for (f in failed) say(sprintf("  rep %d: %s", f$rep[1], f$error[1]))
  stop("all repetitions failed -- see messages above; no CSV written")
}
raw <- do.call(rbind, lapply(parts, read.csv))
write.csv(raw, file.path(OUT, "rhc_reps_raw.csv"), row.names = FALSE)

## aggregate per method x r: median of first r reps; also the partition SD (spread
## of the single-split estimates across all reps -- "how much the split matters")
agg <- list()
for (nm in names(M)) {
  s  <- raw[raw$method == nm, ]; s <- s[order(s$rep), ]
  th <- s$estimate; vr <- s$se^2
  for (rr in GRID) {
    cc <- combine(th, vr, rr); est <- unname(cc["est"]); se <- unname(cc["se"])
    agg[[length(agg) + 1]] <- data.frame(method = nm, r = rr, estimate = est, se = se,
      ci_lo = est - 1.96 * se, ci_hi = est + 1.96 * se, ci_width = 2 * 1.96 * se,
      partition_sd = stats::sd(th), row.names = NULL)
  }
}
agg <- do.call(rbind, agg)
write.csv(agg, file.path(OUT, "rhc_reps_table.csv"), row.names = FALSE)

say(sprintf("\nDONE in %.1f min; reps ok: %d/%d%s",
            as.numeric(difftime(Sys.time(), t0, units = "mins")), length(parts), R_MAX,
            if (length(failed)) sprintf("  (FAILED: %s)",
              paste(vapply(failed, function(f) as.character(f$rep[1]), ""), collapse = ",")) else ""))
say(sprintf("wrote %s and %s", file.path(OUT, "rhc_reps_raw.csv"), file.path(OUT, "rhc_reps_table.csv")))
kr <- intersect(c(1, 25, 50, 100), GRID); if (!length(kr)) kr <- sort(unique(c(min(GRID), max(GRID))))
say("\nRHC estimate (risk difference) stabilising with r  [partition_sd = spread of single-split estimates]:")
print(agg[agg$r %in% kr, ], row.names = FALSE, digits = 4)
