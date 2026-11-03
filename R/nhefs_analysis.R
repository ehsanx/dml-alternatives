## =============================================================================
## nhefs_analysis.R -- NHEFS-specific artifacts for the tutorial:
##  (1) the six estimators on one shared SL library, one run (seed 1);
##  (2) a SOFTWARE-COMPARISON table (hand vs tmle / AIPW / DoubleML / tmle3; Table 11, Figure 10);
##  (3) the propensity-overlap diagnostic (Step 1 "diagnose overlap first"; Figures 2 and 9).
## Writes output/nhefs_panel.csv, nhefs_software.{csv,md}, nhefs_overlap_data.csv and nhefs_overlap.png.
source("R/utils_data.R")
source("R/estimators.R")
suppressMessages({ library(ggplot2) })
OUT <- "output"; dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
lib <- c("SL.mean", "SL.glm", "SL.glmnet", "SL.ranger", "SL.earth")
## reproducibility: pin BLAS/OpenMP threads (also avoids the sl3/tmle3 load segfault) and set a
## seed before each full-sample fit -- nuisance_full() takes no seed, so its SuperLearner CV folds
## are otherwise unpinned and g-comp/IPTW/AIPW/TMLE drift run-to-run (matches real_application.R).
Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1")

d <- load_nhefs()
cat("NHEFS n =", length(d$Y), " treated(qsmk=1) =", sum(d$A), "\n")

## ---- (1) the six-estimator panel, one run ---------------------------------
set.seed(1L)
panel <- run_panel(d$Y, d$A, d$X, family_y = d$family_y, SL.library = lib, truth = NA, seed = 1L)
panel <- panel[, c("estimator", "nuisance", "estimate", "se", "ci_lo", "ci_hi")]
write.csv(panel, file.path(OUT, "nhefs_panel.csv"), row.names = FALSE)
dr <- panel[grepl("AIPW|TMLE", panel$estimator), ]
cat(sprintf("DR estimates span [%.3f, %.3f] kg (spread %.3f)\n",
            min(dr$estimate), max(dr$estimate), diff(range(dr$estimate))))

## ---- (2) software comparison: same estimand, four packages -----------------
set.seed(1L)
nuF <- nuisance_full(d$Y, d$A, d$X, lib, d$family_y)
nuX <- nuisance_crossfit(d$Y, d$A, d$X, 5L, lib, d$family_y, seed = 1L)
hand_aipw <- est_aipw(d$Y, d$A, nuF)$estimate
hand_dml  <- est_aipw(d$Y, d$A, nuX)$estimate
hand_tmle <- est_tmle(d$Y, d$A, nuF, range(d$Y))$estimate
soft <- data.frame(source = "hand-coded (this paper)", AIPW = hand_aipw, TMLE = hand_tmle,
                   DML_crossfit = hand_dml, stringsAsFactors = FALSE)
add_soft <- function(src, aipw = NA, tmle = NA, dml = NA)
  soft[nrow(soft) + 1, ] <<- list(src, aipw, tmle, dml)

tryCatch({ suppressMessages(library(tmle)); set.seed(1)
  tm <- tmle(Y = d$Y, A = d$A, W = as.data.frame(d$X), Q.SL.library = lib, g.SL.library = lib, family = "gaussian")
  add_soft("tmle (Gruber & van der Laan)", tmle = tm$estimates$ATE$psi)
}, error = function(e) cat("tmle failed:", conditionMessage(e), "\n"))

tryCatch({ suppressMessages(library(AIPW)); set.seed(1)
  ap <- AIPW$new(Y = d$Y, A = d$A, W = as.data.frame(d$X), Q.SL.library = lib, g.SL.library = lib,
                 k_split = 5, verbose = FALSE)$fit()$summary()
  add_soft("AIPW (Zhong et al.)", aipw = as.numeric(ap$estimates$RD["Estimate"]))
}, error = function(e) cat("AIPW failed:", conditionMessage(e), "\n"))

tryCatch({ suppressMessages({library(DoubleML); library(mlr3); library(mlr3learners)})
  suppressMessages(lgr::get_logger("mlr3")$set_threshold("warn"))
  df <- data.frame(Y = d$Y, A = d$A, as.data.frame(d$X)); xcols <- setdiff(names(df), c("Y", "A"))
  dd <- DoubleMLData$new(df, y_col = "Y", d_cols = "A", x_cols = xcols); set.seed(1)
  irm <- DoubleMLIRM$new(dd, ml_g = lrn("regr.ranger"), ml_m = lrn("classif.ranger", predict_type = "prob"),
                         n_folds = 5, score = "ATE"); irm$fit()
  add_soft("DoubleML (Bach et al.)", dml = as.numeric(irm$coef))
}, error = function(e) cat("DoubleML failed:", conditionMessage(e), "\n"))

tryCatch({ suppressMessages({library(sl3); library(tmle3); library(data.table)}); options(sl3.verbose = FALSE)
  stack <- Stack$new(Lrnr_mean$new(), Lrnr_glm$new(), Lrnr_glmnet$new(), Lrnr_ranger$new(), Lrnr_earth$new())
  sl <- Lrnr_sl$new(learners = stack, metalearner = Lrnr_nnls$new())
  Xdf <- as.data.frame(d$X); dt <- data.table(A = d$A, Y = d$Y, Xdf)
  fit <- tmle3(tmle_ATE(1, 0), dt, list(W = names(Xdf), A = "A", Y = "Y"), list(A = sl, Y = sl))
  add_soft("tmle3 (tlverse)", tmle = as.numeric(fit$summary$tmle_est))
}, error = function(e) cat("tmle3 failed:", conditionMessage(e), "\n"))

write.csv(soft, file.path(OUT, "nhefs_software.csv"), row.names = FALSE)
fmt <- function(x) ifelse(is.na(x), "—", formatC(x, format = "f", digits = 3))
md <- c("# NHEFS software comparison (effect of smoking cessation on weight change, kg)", "",
        "Same estimand, same data; one run of each implementation's doubly-robust estimate. The runs differ in",
        "nuisance-fitting scheme (in-sample or out-of-fold) and in their realized fits. Their spread is compatible",
        "with fitting-scheme differences and single-split noise; this comparison does not isolate the effect of an",
        "individual setting. R/nhefs_matched_repeated.R repeats three implementations over 20 seeds.", "",
        "| Implementation | AIPW | TMLE | DML (cross-fit) |", "|---|---:|---:|---:|")
for (i in seq_len(nrow(soft))) md <- c(md, sprintf("| %s | %s | %s | %s |",
    soft$source[i], fmt(soft$AIPW[i]), fmt(soft$TMLE[i]), fmt(soft$DML_crossfit[i])))
writeLines(md, file.path(OUT, "nhefs_software.md"))

## ---- (3) overlap diagnostic ------------------------------------------------
g <- .fit_sl(d$A, as.data.frame(d$X), binomial(), c("SL.glm", "SL.glmnet", "SL.ranger"))
pp <- data.frame(pi = g(as.data.frame(d$X)),
                 group = factor(d$A, 0:1, c("did not quit", "quit smoking")))
write.csv(pp, file.path(OUT, "nhefs_overlap_data.csv"), row.names = FALSE)
p <- ggplot(pp, aes(pi, fill = group)) + geom_density(alpha = 0.5, colour = NA) +
  scale_fill_manual(values = c("#D55E00", "#0072B2")) +
  labs(title = "NHEFS: estimated propensity overlap (good common support)",
       x = "estimated P(quit smoking | covariates)", y = "density", fill = NULL) +
  theme_bw(base_size = 11) + theme(legend.position = "top")
ggsave(file.path(OUT, "nhefs_overlap.png"), p, width = 7, height = 3.6, dpi = 130)

cat("\n== panel ==\n"); print(panel, row.names = FALSE, digits = 4)
cat("\n== software ==\n"); print(soft, row.names = FALSE, digits = 4)
cat("\nDONE nhefs_analysis\n")
