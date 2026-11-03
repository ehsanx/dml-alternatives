## tmle3_crosscheck.R -- supplementary analysis (not a manuscript exhibit): a tlverse tmle3 arm for the
## package cross-check. tmle3 is cross-validated by construction (CV-TMLE), so it is compared mainly with
## the hand-coded CV-TMLE (cross-fit TMLE), using an sl3 Super Learner stack that mirrors the SuperLearner
## library. Writes output/tmle3_crosscheck.md.
## =============================================================================
source("R/utils_data.R")
source("R/estimators.R")
suppressMessages({ library(sl3); library(tmle3); library(data.table) })
options(sl3.verbose = FALSE)
OUT <- "output"
lib <- c("SL.mean", "SL.glm", "SL.glmnet", "SL.ranger", "SL.earth")

sl3_sl <- function() {
  stack <- Stack$new(Lrnr_mean$new(), Lrnr_glm$new(), Lrnr_glmnet$new(),
                     Lrnr_ranger$new(), Lrnr_earth$new())
  Lrnr_sl$new(learners = stack, metalearner = Lrnr_nnls$new())
}

run_one <- function(d) {
  yr  <- if (all(d$Y %in% 0:1)) c(0, 1) else range(d$Y)
  nuF <- nuisance_full(d$Y, d$A, d$X, lib, d$family_y)
  nuX <- nuisance_crossfit(d$Y, d$A, d$X, 5L, lib, d$family_y, seed = 1L)
  hand_full <- est_tmle(d$Y, d$A, nuF, yr)$estimate
  hand_cf   <- est_tmle(d$Y, d$A, nuX, yr)$estimate
  Xdf <- as.data.frame(d$X)
  dt  <- data.table(A = d$A, Y = d$Y, Xdf)
  nodes <- list(W = names(Xdf), A = "A", Y = "Y")
  sl <- sl3_sl()
  fit <- tmle3(tmle_ATE(1, 0), dt, nodes, list(A = sl, Y = sl))
  t3 <- as.numeric(fit$summary$tmle_est)
  data.frame(dataset = d$name, hand_TMLE_full = hand_full, hand_CV_TMLE = hand_cf,
             tmle3_CV_TMLE = t3, absdiff_vs_handCV = abs(t3 - hand_cf))
}

rows <- list()
for (d in list(load_ihdp(1), load_twins(n_sub = 2000L, seed = 1L))) {
  r <- tryCatch(run_one(d), error = function(e) {
    cat("tmle3 failed on", d$name, ":", conditionMessage(e), "\n"); NULL })
  if (!is.null(r)) rows[[length(rows) + 1]] <- r
}
res <- do.call(rbind, rows)
print(res, row.names = FALSE, digits = 5)

if (!is.null(res)) {
  fmt <- function(x) formatC(x, format = "g", digits = 5)
  md <- c("# tmle3 (tlverse) cross-check", "",
          "`tmle3` is cross-validated by construction (CV-TMLE). We compare it to our hand-coded CV-TMLE",
          "(cross-fit TMLE), with an `sl3` Super Learner stack mirroring our `SuperLearner` library. Small",
          "`absdiff_vs_handCV` confirms the hand-coded cross-fit TMLE reproduces an independent tlverse impl.", "",
          "| dataset | hand TMLE (full) | hand CV-TMLE | tmle3 (CV-TMLE) | abs diff vs hand CV-TMLE |",
          "|---|---:|---:|---:|---:|")
  for (i in seq_len(nrow(res))) md <- c(md, sprintf("| %s | %s | %s | %s | %s |",
      res$dataset[i], fmt(res$hand_TMLE_full[i]), fmt(res$hand_CV_TMLE[i]),
      fmt(res$tmle3_CV_TMLE[i]), fmt(res$absdiff_vs_handCV[i])))
  writeLines(md, file.path(OUT, "tmle3_crosscheck.md"))
  cat("\nwrote", file.path(OUT, "tmle3_crosscheck.md"), "\n")
}
