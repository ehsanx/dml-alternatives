## package_crosscheck.R -- supplementary analysis (not a manuscript exhibit): check that the hand-coded
## estimators reproduce established packages:
##   hand TMLE (full)        vs  tmle::tmle
##   hand AIPW (cross-fit=DML) vs AIPW::AIPW (k_split cross-fitting)
##   hand DML  (cross-fit AIPW) vs DoubleML::DoubleMLIRM (ATE score)
## tmle and AIPW use the same Super Learner library as the hand code; DoubleML uses mlr3 ranger learners.
## Writes output/crosscheck.md and output/crosscheck.csv.
## =============================================================================
source("R/utils_data.R")
source("R/estimators.R")
OUT <- "output"; dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
lib <- c("SL.mean", "SL.glm", "SL.glmnet", "SL.ranger", "SL.earth")

rows <- list()
add <- function(dataset, method, hand, pkg)
  rows[[length(rows) + 1]] <<- data.frame(dataset, method, hand = hand, package = pkg,
                                          abs_diff = abs(hand - pkg))

run_dataset <- function(d, do_dml = TRUE) {
  cat("\n====", d$name, "====\n")
  y_range <- if (all(d$Y %in% 0:1)) c(0, 1) else range(d$Y)
  fam <- d$family_y
  nu_F <- nuisance_full(d$Y, d$A, d$X, lib, fam)
  nu_X <- nuisance_crossfit(d$Y, d$A, d$X, 5L, lib, fam, seed = 1L)
  h_tmle <- est_tmle(d$Y, d$A, nu_F, y_range)$estimate
  h_dml  <- est_aipw(d$Y, d$A, nu_X)$estimate

  ## --- tmle::tmle (TMLE, full-data targeting w/ SL) ---
  tryCatch({
    suppressMessages(library(tmle)); set.seed(1)
    tm <- tmle(Y = d$Y, A = d$A, W = as.data.frame(d$X),
               Q.SL.library = lib, g.SL.library = lib,
               family = if (all(d$Y %in% 0:1)) "binomial" else "gaussian")
    add(d$name, "TMLE  (hand-full vs tmle pkg)", h_tmle, tm$estimates$ATE$psi)
  }, error = function(e) cat("  tmle pkg failed:", conditionMessage(e), "\n"))

  ## --- AIPW::AIPW (cross-fit AIPW = DML) ---
  tryCatch({
    suppressMessages(library(AIPW)); set.seed(1)
    ap <- AIPW$new(Y = d$Y, A = d$A, W = as.data.frame(d$X),
                   Q.SL.library = lib, g.SL.library = lib,
                   k_split = 5, verbose = FALSE)$fit()$summary()
    rd <- ap$estimates$RD["Estimate"]
    add(d$name, "AIPW/DML (hand-crossfit vs AIPW pkg)", h_dml, as.numeric(rd))
  }, error = function(e) cat("  AIPW pkg failed:", conditionMessage(e), "\n"))

  ## --- DoubleML::DoubleMLIRM (DML, ATE) --- (continuous outcome only, mlr3 ranger)
  if (do_dml) tryCatch({
    suppressMessages({library(DoubleML); library(mlr3); library(mlr3learners)})
    suppressMessages(lgr::get_logger("mlr3")$set_threshold("warn"))
    df <- data.frame(Y = d$Y, A = d$A, as.data.frame(d$X))
    xcols <- setdiff(names(df), c("Y", "A"))
    dml_data <- DoubleMLData$new(df, y_col = "Y", d_cols = "A", x_cols = xcols)
    set.seed(1)
    irm <- DoubleMLIRM$new(dml_data,
             ml_g = lrn("regr.ranger"), ml_m = lrn("classif.ranger", predict_type = "prob"),
             n_folds = 5, score = "ATE")
    irm$fit()
    add(d$name, "DML (hand-crossfit vs DoubleML pkg)", h_dml, as.numeric(irm$coef))
  }, error = function(e) cat("  DoubleML pkg failed:", conditionMessage(e), "\n"))
}

run_dataset(load_nhefs(), do_dml = TRUE)                       # continuous -- NHEFS (the software-bridge dataset; matched library)
run_dataset(load_ihdp(1), do_dml = TRUE)                       # continuous
run_dataset(load_twins(n_sub = 2000L, seed = 1L), do_dml = FALSE)  # binary

res <- do.call(rbind, rows)
fmt <- function(x) formatC(x, format = "g", digits = 5)
md <- c("# Package cross-check — hand-coded vs established packages", "",
        "Small `abs_diff` = the explicit step-by-step code reproduces the package. `tmle` and `AIPW`",
        "use the identical Super Learner library; `DoubleML` uses mlr3 ranger (independent learner).", "",
        "| dataset | comparison | hand | package | abs diff |",
        "|---|---|---:|---:|---:|")
for (i in seq_len(nrow(res)))
  md <- c(md, sprintf("| %s | %s | %s | %s | %s |", res$dataset[i], res$method[i],
                      fmt(res$hand[i]), fmt(res$package[i]), fmt(res$abs_diff[i])))
writeLines(md, file.path(OUT, "crosscheck.md"))
write.csv(res, file.path(OUT, "crosscheck.csv"), row.names = FALSE)
cat("\n================ CROSS-CHECK ================\n"); print(res, row.names = FALSE, digits = 5)
cat("\nWrote", file.path(OUT, "crosscheck.md"), "\n")
