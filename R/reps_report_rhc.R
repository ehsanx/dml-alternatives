## =============================================================================
## reps_report_rhc.R -- RHC repetition-study report.
## Reads output/rhc_reps_table.csv (from reps_rhc.R) and writes a compact companion table
## (markdown + LaTeX) at r = 1, 25, 50, 100: each split-based estimator's running median and the SD
## of the single-split estimates across the repetitions (Step 5). RHC has no ground truth, so this is
## about stability, not coverage. Figure 8 is drawn by R/figures/fig_reps.R. Light I/O only.
## Run: Rscript R/reps_report_rhc.R
OUT <- "output"

tab <- read.csv(file.path(OUT, "rhc_reps_table.csv"), stringsAsFactors = FALSE)
morder <- c("DML (cross-fit)", "CV-TMLE (cross-fit)", "DC-AIPW", "DC-TMLE")
mlab   <- c("DML", "CV-TMLE", "DC-AIPW", "DC-TMLE")
mcol   <- c("#0072B2", "#56B4E9", "#E69F00", "#D55E00")   # Okabe-Ito
tab$method <- factor(tab$method, levels = morder)
grid <- sort(unique(tab$r)); R_MAX <- max(grid)
kr   <- intersect(c(1, 25, 50, 100), grid); if (length(kr) < 2) kr <- sort(unique(c(1, R_MAX)))
val  <- function(m, rr, col) { v <- tab[[col]][tab$method == m & tab$r == rr]; if (length(v)) v[1] else NA }

## ---- markdown table --------------------------------------------------------
md <- c(sprintf("# RHC repetition study (n=5735, no ground truth; %d repetitions)", R_MAX), "",
        "Repeated cross-fitting on the single RHC dataset; median-aggregated over the",
        "first r splits. No coverage/bias (observational). `partition_sd` = spread of the",
        "single-split estimates across all repetitions (how much the split matters).", "",
        paste0("| Estimator | ", paste(sprintf("r=%d", kr), collapse = " | "), " | partition_sd |"),
        paste0("|---|", paste(rep("---:", length(kr) + 1), collapse = "|"), "|"))
for (m in morder) md <- c(md, paste0("| ", m, " | ",
  paste(sapply(kr, function(rr) sprintf("%.4f (%.4f, %.4f)", val(m, rr, "estimate"), val(m, rr, "ci_lo"), val(m, rr, "ci_hi"))), collapse = " | "),
  sprintf(" | %.4f |", val(m, kr[1], "partition_sd"))))
writeLines(md, file.path(OUT, "rhc_reps_table.md"))

## ---- LaTeX companion table -------------------------------------------------
tex <- c("\\begin{table}[!ht]", "\\centering",
  sprintf(paste0("\\caption{RHC repetition study (n=5735, 30-day mortality; %d repeated cross-fits, ",
    "median-aggregated as in \\citet{karim2025splits}). RHC has no ground truth, so this reports ",
    "\\emph{stability}, not coverage. ",
    "$\\mathrm{SD}_{\\mathrm{split}}$ is the spread of the single-split estimates across the %d repetitions.}"),
    R_MAX, R_MAX),
  "\\label{tab:rhc-reps}",
  sprintf("\\begin{tabular}{@{}l%sr@{}}", paste(rep("r", length(kr)), collapse = "")),
  "\\toprule",
  paste0("Estimator & ", paste(sprintf("$r{=}%d$", kr), collapse = " & "), " & $\\mathrm{SD}_{\\mathrm{split}}$ \\\\"),
  "\\midrule")
for (k in seq_along(morder)) {
  m <- morder[k]
  cells <- sapply(kr, function(rr) sprintf("%.3f", val(m, rr, "estimate")))
  tex <- c(tex, paste0(mlab[k], " & ", paste(cells, collapse = " & "),
                       sprintf(" & %.4f \\\\", val(m, kr[1], "partition_sd"))))
}
tex <- c(tex, "\\bottomrule", "\\end{tabular}", "\\end{table}")
writeLines(tex, file.path(OUT, "rhc_reps_table.tex"))

## ---- console summary -------------------------------------------------------
cat(sprintf("RHC reps: %d repetitions, r-grid 1..%d\n", R_MAX, R_MAX))
cat("change (r=1 -> r=", R_MAX, ", risk difference):\n", sep = "")
for (m in morder) cat(sprintf("  %-20s % .4f -> % .4f   (change % .4f, partition_sd %.4f)\n",
  m, val(m, 1, "estimate"), val(m, R_MAX, "estimate"),
  val(m, R_MAX, "estimate") - val(m, 1, "estimate"), val(m, 1, "partition_sd")))
cat("wrote rhc_reps_table.{md,tex}\n")
