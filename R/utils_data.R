## =============================================================================
## utils_data.R -- loaders for the benchmark datasets (IHDP, LaLonde/NSW, Twins) and the real datasets
## (NHEFS, RHC). Each loader returns a common structure:
##   list(name, Y, A, X (data.frame of covariates), family_y, outcome_type, truth)
## where `truth` is the benchmark ATE used to grade the estimators (NA for the real datasets).

suppressMessages({ library(haven) })

## Set DATA_DIR to your local benchmark-data folder (with ihdp/, lalonde_nsw/, twins/, rhc/ subfolders).
## NHEFS needs no download (causaldata::nhefs). Override with the DML_DATA_DIR environment variable.
DATA_DIR <- Sys.getenv("DML_DATA_DIR", "data")

.impute_median <- function(df) {
  for (j in seq_along(df)) {
    x <- df[[j]]; if (anyNA(x)) df[[j]][is.na(x)] <- median(x, na.rm = TRUE)
  }
  df
}

## ---- IHDP (semi-synthetic; TRUE ATE known from mu0, mu1) --------------------
## columns (no header): treatment, y_factual, y_cfactual, mu0, mu1, x1..x25
load_ihdp <- function(realization = 1L, dir = file.path(DATA_DIR, "ihdp")) {
  f <- file.path(dir, sprintf("ihdp_npci_%d.csv", realization))
  d <- read.csv(f, header = FALSE)
  names(d) <- c("treatment", "y_factual", "y_cfactual", "mu0", "mu1",
                paste0("x", 1:25))
  list(name = sprintf("IHDP (realization %d)", realization),
       Y = d$y_factual, A = d$treatment,
       X = d[, paste0("x", 1:25)],
       family_y = gaussian(), outcome_type = "continuous",
       truth = mean(d$mu1 - d$mu0))          # true sample ATE
}

## ---- LaLonde / NSW ---------------------------------------------------------
## setting = "experimental" : NSW randomized sample (nsw_dw)  -> truth = randomized diff
##           "cps" / "psid" : NSW treated + observational controls -> recover the experimental truth
load_lalonde <- function(setting = c("experimental", "cps", "psid"),
                         dir = file.path(DATA_DIR, "lalonde_nsw")) {
  setting <- match.arg(setting)
  covs <- c("age", "education", "black", "hispanic", "married", "nodegree", "re74", "re75")
  dw <- as.data.frame(read_dta(file.path(dir, "nsw_dw.dta")))
  ## experimental benchmark ATE (randomized -> unadjusted difference is unbiased)
  truth <- mean(dw$re78[dw$treat == 1]) - mean(dw$re78[dw$treat == 0])
  if (setting == "experimental") {
    d <- dw
  } else {
    ctrl_file <- if (setting == "cps") "cps_controls.dta" else "psid_controls.dta"
    ctrl <- as.data.frame(read_dta(file.path(dir, ctrl_file)))
    treated <- dw[dw$treat == 1, ]
    keep <- c("treat", covs, "re78")
    d <- rbind(treated[, keep], ctrl[, keep])
  }
  list(name = paste0("LaLonde/NSW (", setting, ")"),
       Y = d$re78, A = d$treat, X = .impute_median(d[, covs]),
       family_y = gaussian(), outcome_type = "continuous", truth = truth)
}

## ---- Twins (both potential outcomes stored; TRUE ATE known) ----------------
## The file gives 30 covariates + outcome(t=0), outcome(t=1) (1-yr mortality of the
## lighter/heavier twin). There is no treatment column: as in the CEVAE/GANITE
## benchmark we SIMULATE a confounded assignment T | X, reveal only the factual
## outcome, and hide the counterfactual.
## twins_sample(): the fixed semi-synthetic sample BEFORE treatment assignment.
## Returns y0, y1 (both potential 1-yr-mortality outcomes), X, and truth = mean(y1-y0).
twins_sample <- function(dir = file.path(DATA_DIR, "twins"),
                         n_sub = NULL, sample_seed = 1L) {
  d <- read.csv(gzfile(file.path(dir, "Twin_data.csv.gz")),
                header = TRUE, check.names = FALSE, quote = "")
  ## header fields are wrapped in literal quote characters (the two outcome columns
  ## even use a smart closing quote); strip all quote-like chars from the names.
  names(d) <- gsub("[\"'‘’“”]", "", names(d))
  ## Outcome = 1-year mortality. The raw columns store day-of-death with 9999 = survived;
  ## GANITE's loader defines the binary outcome as (raw < 9999)  =>  died = 1.
  y0 <- as.integer(d[["outcome(t=0)"]] < 9999)
  y1 <- as.integer(d[["outcome(t=1)"]] < 9999)
  Xall <- d[, setdiff(names(d), c("outcome(t=0)", "outcome(t=1)")), drop = FALSE]
  ok <- !(is.na(y0) | is.na(y1)); y0 <- y0[ok]; y1 <- y1[ok]; Xall <- Xall[ok, , drop = FALSE]
  Xall <- .impute_median(as.data.frame(lapply(Xall, as.numeric)))
  if (!is.null(n_sub) && n_sub < length(y0)) {
    set.seed(sample_seed); idx <- sample(length(y0), n_sub)
    y0 <- y0[idx]; y1 <- y1[idx]; Xall <- Xall[idx, , drop = FALSE]
  }
  list(y0 = y0, y1 = y1, X = Xall, truth = mean(y1 - y0))
}

## twins_assign(): SIMULATE a confounded treatment on a fixed sample, reveal only the
## factual outcome. Varying `seed` (with the sample held fixed) gives Monte-Carlo replicates
## whose estimand `truth` is constant -> valid coverage/bias/RMSE study.
twins_assign <- function(s, seed = 1L, conf_strength = 0.75) {
  set.seed(seed)
  Xs  <- scale(as.matrix(s$X)); w <- rnorm(ncol(Xs))
  lin <- as.numeric(Xs %*% w) / sqrt(ncol(Xs))
  A   <- rbinom(length(s$y0), 1, plogis(conf_strength * lin))
  Y   <- A * s$y1 + (1 - A) * s$y0
  list(name = "Twins (simulated confounding)", Y = Y, A = A, X = s$X,
       family_y = binomial(), outcome_type = "binary (risk difference)", truth = s$truth)
}

## fixed sample + one treatment draw (used by comparison_all3.R). It calls set.seed() and draws, so it sets
## the random-number state of the fits that follow; R/positivity_diagnostics.R calls it for that reason only.
load_twins <- function(dir = file.path(DATA_DIR, "twins"),
                       seed = 1L, conf_strength = 0.75, n_sub = NULL) {
  twins_assign(twins_sample(dir, n_sub, sample_seed = seed),
               seed = seed, conf_strength = conf_strength)
}

## ---- NHEFS (real observational epi; NO ground truth) -----------------------
## The canonical Hernan & Robins "What If" example: quitting smoking (qsmk) ->
## 1971-1982 weight change (wt82_71). truth is unknown (observational); the
## benchmark point estimate in the literature is ~3.4-3.5 kg. Use it to show the
## *workflow* and estimator agreement, not bias/coverage grading.
load_nhefs <- function() {
  suppressMessages(library(causaldata))
  d <- causaldata::nhefs
  covs <- c("sex", "age", "race", "education", "smokeintensity",
            "smokeyrs", "exercise", "active", "wt71")
  keep <- stats::complete.cases(d[, c("qsmk", "wt82_71", covs)])
  d <- d[keep, ]
  list(name = "NHEFS (smoking cessation -> weight change)",
       Y = d$wt82_71, A = as.integer(d$qsmk),
       X = .impute_median(as.data.frame(lapply(d[, covs], as.numeric))),
       family_y = gaussian(), outcome_type = "continuous (kg)",
       truth = NA_real_)          # observational: no known ATE
}

## ---- RHC / Right Heart Catheterization (real observational epi; NO truth) ---
## Connors et al. 1996. Treatment swang1 (RHC within 24h), outcome = 30-day mortality
## (dth30) by default, or overall death. ~49 baseline covariates after dropping IDs,
## dates, post-treatment vars, and the 3 high-missingness columns (cat2/adld3p/urin1).
## The canonical dataset behind Mondol & Karim (2025); use for the real-data workflow.
load_rhc <- function(dir = file.path(DATA_DIR, "rhc"), outcome = c("dth30", "death")) {
  outcome <- match.arg(outcome)
  d <- read.csv(file.path(dir, "rhc.csv"), stringsAsFactors = TRUE)
  d[[1]] <- NULL                                   # drop the unnamed row-index column
  A <- as.integer(d$swang1 == "RHC")
  Y <- as.integer(d[[outcome]] == "Yes")
  drop <- c("swang1", "death", "dth30", "sadmdte", "dschdte", "dthdte", "lstctdte",
            "t3d30", "ptid", "cat2", "adld3p", "urin1")  # treatment, outcomes, dates/ids, high-NA
  covs <- droplevels(d[, setdiff(names(d), drop), drop = FALSE])
  X <- as.data.frame(model.matrix(~ . - 1, data = covs))  # one-hot the factors
  colnames(X) <- make.names(colnames(X), unique = TRUE)   # ranger rejects spaces/&//$ in names
  list(name = paste0("RHC (right heart cath -> ", outcome, ")"),
       Y = Y, A = A, X = X,
       family_y = binomial(), outcome_type = "binary (mortality risk difference)",
       truth = NA_real_)          # observational: no known ATE
}
