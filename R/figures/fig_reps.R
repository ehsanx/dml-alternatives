## fig_reps.R -- Figure 8: the RHC repetition study (rhc_reps.pdf): estimate and 95% CI width against
## the number of repetitions r, for the four split-based estimators.
## Reads docs/data/rhc_reps_table.csv.
##   Rscript R/figures/fig_reps.R [output folder]         (default: output/figures)
##
## Drawn at the final printed size (5.147 in wide, scale 1) with cairo_pdf, so fonts are embedded and a
## label of cex x pointsize prints at exactly that size.

## ---- paths: the repository root is the nearest folder above that holds R/estimators.R ----
REPO <- getwd()
while (!file.exists(file.path(REPO, "R", "estimators.R")) && dirname(REPO) != REPO) REPO <- dirname(REPO)
DATA <- Sys.getenv("DML_OUT", file.path(REPO, "docs", "data"))             # frozen inputs, read only
args <- commandArgs(trailingOnly = TRUE)
IMG  <- if (length(args)) args[1] else file.path(REPO, "output", "figures")  # output folder
dir.create(IMG, showWarnings = FALSE, recursive = TRUE)
FONT <- Sys.getenv("DML_FONT", "Arial")   # without Arial, fontconfig substitutes a metric-compatible font

## ---- data ---------------------------------------------------------------------
tab <- read.csv(file.path(DATA, "rhc_reps_table.csv"), stringsAsFactors = FALSE)
morder <- c("DML (cross-fit)", "CV-TMLE (cross-fit)", "DC-AIPW", "DC-TMLE")
mlab   <- c("DML", "CV-TMLE", "DC-AIPW", "DC-TMLE")
mcol   <- c("#0072B2", "#56B4E9", "#E69F00", "#D55E00")   # Okabe-Ito
stopifnot(all(morder %in% tab$method))
grid <- sort(unique(tab$r))
traces <- function(col) lapply(morder, function(m) {
  s <- tab[tab$method == m, ]; s <- s[order(s$r), ]; list(r = s$r, y = s[[col]]) })
EST <- traces("estimate"); CIW <- traces("ci_width")
yrange <- function(S) range(unlist(lapply(S, `[[`, "y")))

## ---- final-size geometry ------------------------------------------------------
## \linewidth = 372 pt = 5.147 in = 370.6 bp. cairo_pdf truncates the page to whole bp, so draw on
## exactly 370 bp: LaTeX then scales by 370.6/370 = 1.0016 (sizes print >= the sizes set here).
FIG_W  <- 370 / 72     # in (5.139); include at width=\linewidth
FIG_H  <- 198 / 72     # in (2.75)
LEG_H  <- 0.25         # in; legend row under both panels
PS     <- 9            # pt; axis titles, panel titles, legend
CEX_TK <- 8 / PS       # tick labels 8 pt
LWD    <- 1.5          # trace width 1.5/96 in = 1.1 pt
LH     <- PS * 1.2 / 72                  # one margin line, in (10.8 pt)
TCL    <- -0.30        # tick length, lines
MGP_X  <- 0.25         # x tick-label line
MGP_Y  <- 0.50         # y tick-label line (labels horizontal, right-aligned here)
GAP_Y  <- 3.0 / 72     # in; clear gap between the widest y tick label and the y-axis title
DESC   <- 0.22 * PS / 72                 # in; descender of the 9 pt title glyphs "(", ")", "p"
XLAB_LINE <- 1.30      # x-axis title line
MAIN_LINE <- 0.35      # panel title line
MAR_TOP   <- 1.40
MAR_BOT   <- 2.40
XLAB <- "number of repetitions  r"
MAIN <- c("RHC: estimate stabilising with r", "RHC: 95% CI width vs r")
YLAB <- c("ATE estimate (risk difference)", "95% CI width")

## the y ticks R draws for each panel (axTicks after the same plot window; no file written)
yticks <- function(yl) { pdf(NULL); on.exit(dev.off())
  plot.new(); plot.window(xlim = range(grid), ylim = yl, log = "x"); axTicks(2) }
YL <- list(yrange(EST), yrange(CIW))
YT <- lapply(YL, yticks)

out <- file.path(IMG, "rhc_reps.pdf")
cairo_pdf(out, width = FIG_W, height = FIG_H, pointsize = PS, family = FONT)
ylb <- par("ylbias")   # R's baseline offset for margin text (0.2 line)
## y-title line per panel: tick-label line + widest tick label + clear gap (+ descender)
ylab_line <- sapply(YT, function(yt)
  MGP_Y + (max(strwidth(format(yt), units = "inches", cex = CEX_TK)) + GAP_Y + DESC) / LH - ylb)
mar_left  <- ylab_line + ylb + 0.75 * PS / 72 / LH + 0.25   # title cap height + 0.25 line of paper
mar_right <- c(1.00, 0.40)                                  # gutter between panels | right edge
## equal plot-region widths: share what is left after the margins
plot_w  <- (FIG_W - LH * sum(mar_left + mar_right)) / 2
col_w   <- plot_w + LH * (mar_left + mar_right)
layout(matrix(c(1, 2, 3, 3), 2, byrow = TRUE), widths = col_w, heights = c(FIG_H - LEG_H, LEG_H))
## a 2 x 2 layout silently sets cex = 0.83; restore 1 so text and margin lines keep their sizes
par(oma = c(0, 0, 0, 0), cex = 1, las = 1, lend = "butt")

panel <- function(i, S) {
  par(mar = c(MAR_BOT, mar_left[i], MAR_TOP, mar_right[i]), cex = 1)
  plot(NA, xlim = range(grid), ylim = YL[[i]], log = "x", axes = FALSE, ann = FALSE)
  xt <- axTicks(1)
  axis(1, at = xt, labels = format(xt, trim = TRUE), cex.axis = CEX_TK, mgp = c(0, MGP_X, 0), tcl = TCL)
  axis(2, at = YT[[i]], labels = format(YT[[i]]), cex.axis = CEX_TK, mgp = c(0, MGP_Y, 0), tcl = TCL, las = 1)
  box()
  mtext(YLAB[i], side = 2, line = ylab_line[i], las = 0)
  mtext(XLAB, side = 1, line = XLAB_LINE)
  mtext(MAIN[i], side = 3, line = MAIN_LINE, font = 2)
  for (k in seq_along(S)) lines(S[[k]]$r, S[[k]]$y, col = mcol[k], lwd = LWD)
}
panel(1, EST)
panel(2, CIW)

## legend: one row, centred under the span of the two plot regions
par(mar = c(0, mar_left[1], 0, mar_right[2]), cex = 1)
plot.new()
legend("center", legend = mlab, col = mcol, lwd = LWD, horiz = TRUE, bty = "n",
       seg.len = 2.2, x.intersp = 0.5, text.width = NA, xpd = NA)
invisible(dev.off())
cat(sprintf("wrote %s (%.3f x %.3f in; plot regions %.2f in wide; y-title lines %.2f / %.2f)\n",
            out, FIG_W, FIG_H, plot_w, ylab_line[1], ylab_line[2]))
