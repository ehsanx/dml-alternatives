## =============================================================================
## lalonde_seed_sweep.R -- seed sensitivity of the LaLonde/NSW-PSID estimates (poor overlap).
## On one split the cross-fit and double-cross-fit TMLE estimates lie far below the full-sample
## estimates. This varies the fold/SL seed to show that the exact values move but the sign, and
## the position far below the full-sample estimates, do not (Table 8, Step 6). The $1,794 experimental
## estimate concerns a different population and is not a benchmark for these estimates.
## Writes output/lalonde_seed_sweep.{csv,md}.  Run: Rscript R/lalonde_seed_sweep.R
## (env SEEDS_N, default 8).  Thread-pin: OMP_NUM_THREADS=1 ... Rscript ...
source("R/utils_data.R")
source("R/estimators.R")
Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1")
OUT <- "output"; dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
lib <- if (nzchar(Sys.getenv("SL_FAST"))) c("SL.mean", "SL.glm") else default_SL

d <- load_lalonde("psid")
N <- as.integer(Sys.getenv("SEEDS_N", "8"))
cat(sprintf("LaLonde/PSID seed sweep: n=%d, truth=%.0f, %d seeds\n", length(d$Y), d$truth, N))

rows <- list()
for (s in seq_len(N)) {
  t0 <- Sys.time(); set.seed(s)
  r <- run_panel(d$Y, d$A, d$X, family_y = d$family_y, SL.library = lib, truth = d$truth, seed = s)
  r$seed <- s
  rows[[s]] <- r[, c("seed", "estimator", "nuisance", "estimate")]
  cat(sprintf("  seed %d/%d  %.0fs\n", s, N, as.numeric(difftime(Sys.time(), t0, units = "secs"))))
}
res <- do.call(rbind, rows)
write.csv(res, file.path(OUT, "lalonde_seed_sweep.csv"), row.names = FALSE)

## per-estimator median + range across seeds
smry <- do.call(rbind, lapply(split(res, res$estimator), function(x) data.frame(
  estimator = x$estimator[1], nuisance = x$nuisance[1],
  median = round(median(x$estimate)), min = round(min(x$estimate)), max = round(max(x$estimate)))))
smry <- smry[order(smry$median), ]
md <- c("# LaLonde/NSW-PSID seed sweep (poor overlap; no known truth for this population)", "",
        sprintf("%d seeds; median [min, max] of each estimator across seeds.", N), "",
        "| estimator | nuisance | median | min | max |", "|---|---|--:|--:|--:|")
for (i in seq_len(nrow(smry))) md <- c(md, sprintf("| %s | %s | %d | %d | %d |",
  smry$estimator[i], smry$nuisance[i], smry$median[i], smry$min[i], smry$max[i]))
writeLines(md, file.path(OUT, "lalonde_seed_sweep.md"))
cat("\n"); print(smry, row.names = FALSE)
cat("\nwrote output/lalonde_seed_sweep.{csv,md}\n")
