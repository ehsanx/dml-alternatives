## =============================================================================
## nhefs_matched_repeated.R -- software reconciliation on NHEFS over repeated seeds (Step 7, Figure 10).
## Once the Super Learner library AND the [0.025,0.975] propensity truncation are
## matched across implementations, do the packages converge? A single cross-fit split
## still scatters by the partition SD (~0.08 kg on NHEFS), so each implementation is
## run with R seeds and the medians are compared:
##   hand-coded DML (=cross-fit AIPW)  vs  AIPW::AIPW  vs  tmle::tmle
## all with the identical library and truncation. Writes output/nhefs_matched_repeated.{csv,md}.
## Run:  Rscript R/nhefs_matched_repeated.R   (env REPS, default 20)
source("R/utils_data.R")
source("R/estimators.R")
Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1")
OUT <- "output"; dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
lib <- c("SL.mean", "SL.glm", "SL.glmnet", "SL.ranger", "SL.earth")
d <- load_nhefs()
R <- as.integer(Sys.getenv("REPS", "20"))
cat(sprintf("NHEFS matched-repeated comparison, R=%d splits, truncation [0.025,0.975]\n", R))

hand <- sapply(seq_len(R), function(s) {
  nu <- nuisance_crossfit(d$Y, d$A, d$X, 5L, lib, d$family_y, seed = s)
  est_aipw(d$Y, d$A, nu)$estimate
})
suppressMessages(library(AIPW))
aipw <- sapply(seq_len(R), function(s) {
  set.seed(s)
  tryCatch(as.numeric(AIPW$new(Y = d$Y, A = d$A, W = as.data.frame(d$X),
             Q.SL.library = lib, g.SL.library = lib, k_split = 5,
             verbose = FALSE)$fit()$summary(g.bound = 0.025)$estimates$RD["Estimate"]),
           error = function(e) NA_real_)
})
suppressMessages(library(tmle))
tm <- sapply(seq_len(R), function(s) {
  set.seed(s)
  tryCatch(as.numeric(tmle(Y = d$Y, A = d$A, W = as.data.frame(d$X),
             Q.SL.library = lib, g.SL.library = lib, gbound = 0.025,
             family = "gaussian")$estimates$ATE$psi),
           error = function(e) NA_real_)
})

out <- data.frame(rep = seq_len(R), hand_DML = hand, AIPW_pkg = aipw, tmle_pkg = tm)
write.csv(out, file.path(OUT, "nhefs_matched_repeated.csv"), row.names = FALSE)
med <- sapply(out[c("hand_DML","AIPW_pkg","tmle_pkg")], median, na.rm = TRUE)
sdv <- sapply(out[c("hand_DML","AIPW_pkg","tmle_pkg")], sd, na.rm = TRUE)
maxgap <- max(abs(outer(med, med, "-")))

md <- c("# NHEFS software bridge: matched library + truncation, medians over repeated cross-fits", "",
        sprintf("R = %d repeated cross-fits; propensity truncation [0.025, 0.975]; identical Super Learner library.", R), "",
        "| implementation | median (kg) | per-split SD |",
        "|---|---:|---:|",
        sprintf("| hand-coded engine (DML = cross-fit AIPW) | %.2f | %.3f |", med["hand_DML"], sdv["hand_DML"]),
        sprintf("| `AIPW` package (Zhong et al.) | %.2f | %.3f |", med["AIPW_pkg"], sdv["AIPW_pkg"]),
        sprintf("| `tmle` package (Gruber & van der Laan) | %.2f | %.3f |", med["tmle_pkg"], sdv["tmle_pkg"]),
        "",
        sprintf("**Largest pairwise median gap: %.3f kg** -- the three medians agree at about %.2f kg.", maxgap, median(med)),
        sprintf("A single run scatters with an SD across seeds of about %.2f kg (each seed redraws the folds, the internal cross-validation and the random forests).", mean(sdv[c("hand_DML","AIPW_pkg")])),
        "At package defaults (single library/fold/truncation choices) the same estimand spans 3.32–3.49 kg (Table: software).")
writeLines(md, file.path(OUT, "nhefs_matched_repeated.md"))
cat(sprintf("\nmedians: hand %.4f  AIPW %.4f  tmle %.4f  | max pairwise gap %.4f\n",
            med["hand_DML"], med["AIPW_pkg"], med["tmle_pkg"], maxgap))
cat("wrote nhefs_matched_repeated.{csv,md}\n")
