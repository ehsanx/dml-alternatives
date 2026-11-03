## =============================================================================
## positivity_diagnostics.R -- show WHERE the doubly-robust estimators break.
## (1) Estimated-propensity overlap by treatment group for IHDP / LaLonde-PSID.
## (2) Trimming demonstration on LaLonde-PSID: re-estimate the ATE after restricting
##     to the common-support region, showing how sensitive the answer is to overlap.
## Writes output/positivity_overlap.png and output/positivity_trimming_lalonde.{csv,md}
## =============================================================================
source("R/utils_data.R")
source("R/estimators.R")
suppressMessages(library(ggplot2))
OUT <- "output"; dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
lib_ps <- c("SL.glm", "SL.glmnet", "SL.ranger")

## ---- (1) propensity overlap plot -------------------------------------------
prop_of <- function(d) {
  g <- .fit_sl(d$A, as.data.frame(d$X), binomial(), lib_ps)
  data.frame(dataset = d$name, pi = g(as.data.frame(d$X)),
             group = factor(d$A, levels = c(0, 1), labels = c("control", "treated")))
}
set.seed(1)
dl <- list(load_ihdp(1), load_lalonde("psid"))
## Twins is not analysed here, but its loader calls set.seed() and draws, which fixes the random-number
## state of the fits below. Keep the call: positivity_overlap.csv regenerates exactly only with it (and
## with the Twins file under DATA_DIR).
invisible(load_twins(n_sub = 2000L))
pp <- do.call(rbind, lapply(dl, prop_of))
p <- ggplot(pp, aes(pi, fill = group)) +
  geom_density(alpha = 0.5, colour = NA) +
  facet_wrap(~dataset, scales = "free") +
  labs(title = "Estimated propensity overlap by treatment group",
       subtitle = "Good overlap (IHDP) vs poor overlap (LaLonde/NSW-PSID)",
       x = "estimated P(A = 1 | X)", y = "density", fill = NULL) +
  theme_bw(base_size = 11) + theme(legend.position = "top")
ggsave(file.path(OUT, "positivity_overlap.png"), p, width = 7.5, height = 3.6, dpi = 130)
write.csv(pp, file.path(OUT, "positivity_overlap.csv"), row.names = FALSE, quote = FALSE)  # density source for fig:positivity
cat("wrote positivity_overlap.png + positivity_overlap.csv\n")

## ---- (2) trimming demonstration on LaLonde-PSID ----------------------------
d <- load_lalonde("psid")
est_on <- function(idx) {
  dt <- list(Y = d$Y[idx], A = d$A[idx], X = d$X[idx, , drop = FALSE])
  nu <- nuisance_crossfit(dt$Y, dt$A, dt$X, 5L, default_SL, d$family_y, seed = 1L)
  e  <- est_aipw(dt$Y, dt$A, nu)
  list(n = sum(idx), est = e$estimate, ci = e$ci)
}
## estimated propensity on full data to define the trimming region
nu_full <- nuisance_crossfit(d$Y, d$A, d$X, 5L, default_SL, d$family_y, seed = 1L)
pi <- nu_full$pi
rules <- list("none (all data)" = rep(TRUE, length(pi)),
              "trim to [0.05, 0.95]" = pi >= 0.05 & pi <= 0.95,
              "trim to [0.10, 0.90]" = pi >= 0.10 & pi <= 0.90)
tab <- do.call(rbind, lapply(names(rules), function(nm) {
  r <- est_on(rules[[nm]])
  data.frame(rule = nm, n_kept = r$n, pct_kept = round(100 * r$n / length(pi)),
             AIPW_DML = round(r$est), ci_lo = round(r$ci[1]), ci_hi = round(r$ci[2]),
             covers_1794 = (1794 >= r$ci[1] & 1794 <= r$ci[2]))
}))
write.csv(tab, file.path(OUT, "positivity_trimming_lalonde.csv"), row.names = FALSE)
## $1,794 is the experimental estimate for a different population, so the report below does not use
## covers_1794 to judge the trimmed estimates.
md <- c("# Trimming -- LaLonde/NSW-PSID (poor overlap)", "",
        "Cross-fit AIPW (=DML) re-estimated after restricting to the estimated common-support region.",
        "Trimming changes the population, and with it the estimand; the trimmed intervals are wide.",
        "", "| rule | n kept | % kept | AIPW/DML | 95% CI |",
        "|---|--:|--:|--:|:--:|")
for (i in seq_len(nrow(tab))) md <- c(md, sprintf(
  "| %s | %d | %d%% | %d | (%d, %d) |", tab$rule[i], tab$n_kept[i], tab$pct_kept[i],
  tab$AIPW_DML[i], tab$ci_lo[i], tab$ci_hi[i]))
writeLines(md, file.path(OUT, "positivity_trimming_lalonde.md"))
cat("\n"); print(tab, row.names = FALSE); cat("\nDONE positivity\n")
