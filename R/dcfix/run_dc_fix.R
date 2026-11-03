## run_dc_fix.R -- rerun every double-cross-fit (DC-AIPW / DC-TMLE) result with the splitting rule of
## R/estimators.R.
##
## Only the DC rows are rerun; the other estimators do not use that splitting rule. Each study reuses the
## fold seeds of the script that produced it (in brackets), so the DC estimators see that script's partitions:
##   nhefs_r1      NHEFS, seed 1                      (nhefs_analysis.R, real_application.R)
##   nhefs_reps    NHEFS, reps j = 1..100, seed 100+j (R/nhefs_panel_reps.R)
##   rhc_r1        RHC (30-day death), seed 1         (real_application.R)
##   rhc_reps      RHC, reps j = 1..100, seed 100+j   (R/reps_rhc.R)
##   ihdp_grade    IHDP realizations 1..50 x seeds s = 1..20, seed 100+s (R/grade_ihdp.R)
##   ihdp_r1       IHDP realization 1, seed 1         (comparison_all3.R)
##   lalonde_r1    LaLonde/NSW-PSID, seed 1           (comparison_all3.R)
##   lalonde_sweep LaLonde/NSW-PSID, seeds 1..8       (lalonde_seed_sweep.R)
##   nhefs_trunc   NHEFS truncation audit: pre-clipping propensity range and count clipped
##                 for the full-sample and seed-1 cross-fit fits (DC fits record theirs too)
## (Twins, in comparison_all3.R, is not used in the manuscript and is not rerun.)
##
## Output goes to a versioned folder results/<RUN_ID>/, where RUN_ID carries the md5 of the engine file
## (computed on LF line endings, so it does not depend on the checkout). Parts are resumed only inside that
## folder, so results from a different engine are never reused. The folder's MANIFEST.txt records the
## engine hash, git commit, R and package versions, seeds, settings and timing. Run from the repo root:
##   Rscript R/dcfix/run_dc_fix.R            # full run (NCORES workers, default 10)
##   DCFIX_SMOKE=1 Rscript R/dcfix/run_dc_fix.R   # 1 unit per study, tiny library
## Benchmark data: set DML_DATA_DIR (see R/utils_data.R). NHEFS needs no download.
## =============================================================================
suppressMessages({ library(SuperLearner); library(parallel) })
Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1")
ROOT   <- normalizePath(Sys.getenv("DML_ROOT", getwd()), winslash = "/")
ESTIM  <- file.path(ROOT, "R", "estimators.R")
LOADER <- file.path(ROOT, "R", "utils_data.R")
DDIR   <- Sys.getenv("DML_DATA_DIR", file.path(ROOT, "data"))
SMOKE  <- nzchar(Sys.getenv("DCFIX_SMOKE"))
NCORES <- as.integer(Sys.getenv("NCORES", "10"))
lib    <- if (SMOKE) c("SL.mean", "SL.glm") else
          c("SL.mean", "SL.glm", "SL.glmnet", "SL.ranger", "SL.earth")   # the full library of the analysis scripts

.lf <- tempfile(); local({ con <- file(.lf, "wb"); writeLines(readLines(ESTIM, warn = FALSE), con, sep = "\n"); close(con) })
engine_md5 <- unname(tools::md5sum(.lf))   # LF-normalized, so a CRLF checkout gives the same id
RUN_ID <- sprintf("dcfix1_%s%s", substr(engine_md5, 1, 10), if (SMOKE) "_smoke" else "")
OUT    <- file.path(ROOT, "results", RUN_ID)
PARTS  <- file.path(OUT, "parts")
dir.create(PARTS, recursive = TRUE, showWarnings = FALSE)
say <- function(...) { cat(sprintf(...), "\n"); flush.console() }

## ---- task list (most expensive first, for load balancing) ---------------------
tasks <- c(
  lapply(seq_len(100), function(j) list(study = "rhc_reps",   unit = sprintf("rep%03d", j), seed = 100L + j)),
  list(list(study = "rhc_r1", unit = "r1", seed = 1L)),
  lapply(seq_len(8),   function(s) list(study = "lalonde_sweep", unit = sprintf("seed%02d", s), seed = s)),
  list(list(study = "lalonde_r1", unit = "r1", seed = 1L)),
  lapply(seq_len(100), function(j) list(study = "nhefs_reps", unit = sprintf("rep%03d", j), seed = 100L + j)),
  list(list(study = "nhefs_r1", unit = "r1", seed = 1L), list(study = "nhefs_trunc", unit = "audit", seed = 1L)),
  unlist(lapply(seq_len(50), function(i) lapply(seq_len(20), function(s)
    list(study = "ihdp_grade", unit = sprintf("real%02d_seed%02d", i, s), seed = 100L + s, realization = i))),
    recursive = FALSE),
  list(list(study = "ihdp_r1", unit = "r1", seed = 1L, realization = 1L)))
if (SMOKE) tasks <- tasks[!duplicated(vapply(tasks, `[[`, "", "study"))]

## ---- one task: the DC nuisances, then DC-AIPW and DC-TMLE -----------------------
one_task <- function(t) {
  pf <- file.path(PARTS, t$study, paste0(t$unit, ".csv"))
  if (file.exists(pf)) return(invisible(pf))                    # resume within THIS run only
  dir.create(dirname(pf), recursive = TRUE, showWarnings = FALSE)
  t0 <- Sys.time()
  d <- switch(sub("_.*", "", t$study),
              nhefs = load_nhefs(), rhc = load_rhc(outcome = "dth30"),
              ihdp = load_ihdp(t$realization), lalonde = load_lalonde("psid"))
  truth <- if (is.null(d$truth)) NA_real_ else d$truth
  yr <- range(d$Y)
  if (t$study == "nhefs_trunc") {
    set.seed(t$seed)
    nuF <- nuisance_full(d$Y, d$A, d$X, lib, d$family_y)
    nuX <- nuisance_crossfit(d$Y, d$A, d$X, 5L, lib, d$family_y, seed = t$seed)
    out <- data.frame(study = t$study, unit = c("full sample (seed 1)", "cross-fit V=5 (seed 1)"),
                      clip_raw_min = c(nuF$clip["raw_min"], nuX$clip["raw_min"]),
                      clip_raw_max = c(nuF$clip["raw_max"], nuX$clip["raw_max"]),
                      n_clipped = c(nuF$clip["n_clipped"], nuX$clip["n_clipped"]), row.names = NULL)
  } else {
    nu <- nuisance_double(d$Y, d$A, d$X, 3L, lib, d$family_y, seed = t$seed)
    ea <- est_aipw(d$Y, d$A, nu); et <- est_tmle(d$Y, d$A, nu, yr)
    out <- data.frame(study = t$study, unit = t$unit, seed = t$seed,
                      realization = if (is.null(t$realization)) NA_integer_ else t$realization,
                      method = c("DC-AIPW", "DC-TMLE"),
                      estimate = c(ea$estimate, et$estimate), se = c(ea$se, et$se),
                      ci_lo = c(ea$ci[1], et$ci[1]), ci_hi = c(ea$ci[2], et$ci[2]),
                      truth = truth,
                      tmle_converged = c(NA, et$converged), tmle_score_resid = c(NA, et$score_resid),
                      clip_raw_min = nu$clip["raw_min"], clip_raw_max = nu$clip["raw_max"],
                      n_clipped = nu$clip["n_clipped"], row.names = NULL)
    out$covers <- out$truth >= out$ci_lo & out$truth <= out$ci_hi
  }
  out$seconds <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  write.csv(out, pf, row.names = FALSE)
  invisible(pf)
}

## ---- manifest ------------------------------------------------------------------
git <- function(...) tryCatch(system2("git", c("-C", shQuote(ROOT), ...), stdout = TRUE, stderr = FALSE),
                              error = function(e) NA_character_)
pkgs <- c("SuperLearner", "ranger", "glmnet", "earth", "haven")
man <- c(sprintf("run_id: %s", RUN_ID), sprintf("started: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
         sprintf("engine: R/estimators.R md5 %s (LF line endings; rule: score fold v with g from fold v+1, Q from fold v+2, mod 3)", engine_md5),
         sprintf("git HEAD: %s; engine file uncommitted changes: %s", paste(git("rev-parse", "--short", "HEAD"), collapse = ""),
                 if (length(git("status", "--porcelain", "--", "R/estimators.R"))) "yes" else "no"),
         sprintf("R: %s", R.version.string),
         sprintf("packages: %s", paste(sprintf("%s %s", pkgs, vapply(pkgs, function(p) as.character(packageVersion(p)), "")), collapse = "; ")),
         sprintf("library: %s", paste(lib, collapse = ", ")),
         "settings: V = 3 folds for DC; propensity truncation [0.025, 0.975]; TMLE outcome rescaled by its observed range; SL internal CV V = 5",
         sprintf("workers: %d; smoke run: %s", NCORES, SMOKE),
         "studies and seeds: see the header of R/dcfix/run_dc_fix.R",
         sprintf("tasks: %d", length(tasks)))
mf <- file.path(OUT, "MANIFEST.txt")
if (file.exists(mf)) {                      # same engine hash = same run: log the relaunch
  cat(sprintf("relaunched: %s (%d parts already on disk are reused)", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
              length(Sys.glob(file.path(PARTS, "*", "*.csv")))), file = mf, sep = "\n", append = TRUE)
} else writeLines(man, mf)
say("run %s: %d tasks -> %s", RUN_ID, length(tasks), OUT)

## ---- run ------------------------------------------------------------------------
t0 <- Sys.time()
cl <- makeCluster(NCORES, outfile = file.path(OUT, "workers.log"))
clusterExport(cl, c("ESTIM", "LOADER", "DDIR", "PARTS", "lib", "one_task"))
clusterEvalQ(cl, {
  suppressMessages({ library(SuperLearner); library(ranger); library(glmnet); library(earth) })
  Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1", DML_DATA_DIR = DDIR)
  source(LOADER); source(ESTIM); TRUE
})
## chunk.size = 1: hand out one task at a time, so the slow RHC reps spread over all workers
res <- parLapplyLB(cl, tasks, function(t) tryCatch({ one_task(t); "ok" },
                                                   error = function(e) paste("FAIL", t$study, t$unit, conditionMessage(e))),
                   chunk.size = 1L)
stopCluster(cl)
fails <- unlist(res[res != "ok"])
mins <- as.numeric(difftime(Sys.time(), t0, units = "mins"))

## ---- collect --------------------------------------------------------------------
for (st in unique(vapply(tasks, `[[`, "", "study"))) {
  f <- Sys.glob(file.path(PARTS, st, "*.csv"))
  if (length(f)) write.csv(do.call(rbind, lapply(f, read.csv)), file.path(OUT, paste0(st, ".csv")), row.names = FALSE)
}
cat(c(sprintf("finished: %s (%.1f min); failures: %d", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), mins, length(fails)),
      if (length(fails)) paste("  ", fails)), file = file.path(OUT, "MANIFEST.txt"), sep = "\n", append = TRUE)
say("done in %.1f min; failures: %d", mins, length(fails))
if (length(fails)) print(head(fails, 20))
