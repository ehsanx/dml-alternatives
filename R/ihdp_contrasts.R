## =============================================================================
## ihdp_contrasts.R -- realization-level analyses of the IHDP grading (Step 5).
## Reads output/grade_raw.csv (one row per realization x fold seed x estimator; frozen copy
## in docs/data/) and writes, without refitting any model:
##   output/ihdp_contrasts.csv         paired contrasts between routes (coverage, bias,
##                                     squared error, CI width), with the realization as the
##                                     Monte Carlo unit: difference, paired MCSE, t
##   output/ihdp_aggregated_pilot.csv  exploratory check of the recommended aggregated
##                                     estimator: each realization's 20 single splits combined
##                                     by the median rule (r = 20); one pass/fail per realization
## Run from the repository root:
##   Rscript R/ihdp_contrasts.R
OUT <- "output"
gr  <- read.csv(file.path(OUT, "grade_raw.csv"), stringsAsFactors = FALSE)

## the repeated-split rule (Chernozhukov et al. 2018, Sec. 3.4): median of the estimates, and
## median of (variance + squared distance from that median)
combine <- function(th, vr) {
  p <- stats::median(th); v <- stats::median(vr + (th - p)^2); c(est = p, se = sqrt(v))
}

## ---- paired contrasts: average over the 20 splits within each realization, then pair ----
per_real <- function(e, f) tapply(f(gr[gr$estimator == e, ]), gr$realization[gr$estimator == e], mean)
pairs <- list(c("AIPW", "AIPW = DML"), c("TMLE", "TMLE (CV-TMLE)"), c("AIPW", "AIPW (DC-AIPW)"),
              c("TMLE", "TMLE (DC-TMLE)"), c("AIPW = DML", "AIPW (DC-AIPW)"), c("TMLE (CV-TMLE)", "TMLE (DC-TMLE)"))
fs <- list(coverage = function(s) as.numeric(s$covers), bias = function(s) s$bias,
           mse = function(s) s$bias^2, ci_width = function(s) s$ci_width)
ct <- do.call(rbind, lapply(pairs, function(pq) do.call(rbind, lapply(names(fs), function(nm) {
  d <- per_real(pq[2], fs[[nm]]) - per_real(pq[1], fs[[nm]])
  data.frame(from = pq[1], to = pq[2], quantity = nm, difference = mean(d),
             paired_mcse = sd(d) / sqrt(length(d)), t = mean(d) / (sd(d) / sqrt(length(d))), row.names = NULL)
}))))
write.csv(ct, file.path(OUT, "ihdp_contrasts.csv"), row.names = FALSE)

## ---- exploratory: the aggregated estimator over the 20 stored splits (r = 20) ----------
split_e <- c("AIPW = DML", "TMLE (CV-TMLE)", "AIPW (DC-AIPW)", "TMLE (DC-TMLE)")
pil <- do.call(rbind, lapply(split_e, function(e) {
  s <- gr[gr$estimator == e, ]
  inc <- vapply(split(s, s$realization), function(x) {
    x <- x[order(x$seed), ]; cc <- combine(x$estimate, x$se^2)
    x$truth[1] >= cc[["est"]] - 1.96 * cc[["se"]] && x$truth[1] <= cc[["est"]] + 1.96 * cc[["se"]] }, TRUE)
  data.frame(estimator = e, realizations = length(inc), included = sum(inc), inclusion = mean(inc), row.names = NULL)
}))
write.csv(pil, file.path(OUT, "ihdp_aggregated_pilot.csv"), row.names = FALSE)

print(ct, row.names = FALSE, digits = 3)
print(pil, row.names = FALSE)
