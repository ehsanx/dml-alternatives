## fig_overlap.R -- Figure 2: estimated propensity overlap in NHEFS (nhefs_overlap.pdf).
## Reads docs/data/nhefs_overlap_data.csv: the scores of the overlap diagnostic's own Super Learner fit
## (SL.glm, SL.glmnet, SL.ranger), written by R/nhefs_analysis.R. That fit was not seeded, so the
## frozen file is the reference.
##   Rscript R/figures/fig_overlap.R [output folder]      (default: output/figures)
##
## Drawn at its printed size (width 0.90 x 5.147 in), lettering 9 pt, with embedded fonts.

## ---- paths: the repository root is the nearest folder above that holds R/estimators.R ----
REPO <- getwd()
while (!file.exists(file.path(REPO, "R", "estimators.R")) && dirname(REPO) != REPO) REPO <- dirname(REPO)
DATA <- Sys.getenv("DML_OUT", file.path(REPO, "docs", "data"))             # frozen inputs, read only
args <- commandArgs(trailingOnly = TRUE)
IMG  <- if (length(args)) args[1] else file.path(REPO, "output", "figures")  # output folder
dir.create(IMG, showWarnings = FALSE, recursive = TRUE)
FONT <- Sys.getenv("DML_FONT", "Arial")   # without Arial, fontconfig substitutes a metric-compatible font

ok   <- c(blue = "#0072B2", orange = "#D55E00")       # Okabe-Ito (colour-blind safe)
cols <- c(ok[["orange"]], ok[["blue"]])               # group 1 orange, group 2 blue
LWD  <- 1.5                                           # outline width: 1.125 pt at final size
PT   <- 9                                             # base lettering (pt) at final size

## A density shape closed along y = 0, not by a chord between the end points of the curve.
dshape <- function(x, y, col, lwd = LWD)
  polygon(c(x[1], x, x[length(x)]), c(0, y, 0),
          col = adjustcolor(col, 0.4), border = col, lwd = lwd)

## ---- Figure 2: NHEFS propensity overlap -------------------------------------------
pp  <- read.csv(file.path(DATA, "nhefs_overlap_data.csv"))
den <- tapply(pp$pi, pp$group, density)
cairo_pdf(file.path(IMG, "nhefs_overlap.pdf"), width = 0.90 * 5.147, height = 2.6,
          pointsize = PT, family = FONT)
par(mar = c(3.0, 3.2, 0.5, 0.6), mgp = c(1.95, 0.5, 0), tcl = -0.3, las = 1)
plot(NA, xlim = range(pp$pi), ylim = c(0, max(sapply(den, function(d) max(d$y)))),
     xlab = "estimated P(quit smoking | covariates)", ylab = "density")
grp <- names(den)
for (i in seq_along(den)) dshape(den[[i]]$x, den[[i]]$y, cols[i])
legend("topright", legend = grp, fill = adjustcolor(cols, 0.4), border = cols, bty = "n",
       inset = c(0.01, 0.02))
invisible(dev.off())
cat("wrote", file.path(IMG, "nhefs_overlap.pdf"), "\n")
