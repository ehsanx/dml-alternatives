## fig_ihdp.R -- Figure 7: IHDP calibration (ihdp_calibration.pdf), in a 2 x 2 layout:
##   top:    mean bias | RMSE
##   bottom: 95% CI coverage | mean 95% CI width
## for the eight estimators. Bias and coverage carry 95% Monte Carlo intervals whose standard errors
## treat the 50 outcome realizations as the independent units (grade_table.csv's bias_mcse and
## cov_mcse: the SD of the realization-level means over sqrt(50)). RMSE and mean width are points
## without uncertainty bars. g-computation and IPTW report no interval, so their coverage and width
## are marked "not available" rather than drawn at zero.
##
## Reads docs/data/grade_table.csv, whose columns mean_bias, bias_mcse, rmse, coverage, cov_mcse and
## mean_CI_width summarize grade_raw.csv (8 estimators x 50 realizations x 20 fold partitions, each
## scored as a single split, r = 1).
##   Rscript R/figures/fig_ihdp.R [output folder]         (default: output/figures)
##
## Drawn at its final printed size (5.147 in wide, 8 pt base lettering); cairo_pdf embeds the fonts.

## ---- paths: the repository root is the nearest folder above that holds R/estimators.R ----
REPO <- getwd()
while (!file.exists(file.path(REPO, "R", "estimators.R")) && dirname(REPO) != REPO) REPO <- dirname(REPO)
DATA <- Sys.getenv("DML_OUT", file.path(REPO, "docs", "data"))             # frozen inputs, read only
args <- commandArgs(trailingOnly = TRUE)
IMG  <- if (length(args)) args[1] else file.path(REPO, "output", "figures")  # output folder
dir.create(IMG, showWarnings = FALSE, recursive = TRUE)
FONT <- Sys.getenv("DML_FONT", "Arial")   # without Arial, fontconfig substitutes a metric-compatible font

## --- data: selection, order, labels and groups ---
g <- read.csv(file.path(DATA, "grade_table.csv"), stringsAsFactors = FALSE)
key <- c("g-computation","IPTW","AIPW","TMLE",
         "AIPW = DML","TMLE (CV-TMLE)","AIPW (DC-AIPW)","TMLE (DC-TMLE)")
lab <- c("g-computation","IPTW","AIPW (full)","TMLE (full)",
         "DML (cross-fit)","CV-TMLE (cross-fit)","DC-AIPW","DC-TMLE")
grp <- c("single","single","full","full","cross","cross","double","double")
g <- g[match(key, g$estimator), ]
stopifnot(!anyNA(g$estimator), all(g$n_panels == 1000))
bias <- g$mean_bias; bias_mc <- g$bias_mcse; rmse <- g$rmse
cvg  <- g$coverage;  cvg_mc  <- g$cov_mcse;  wid <- g$mean_CI_width
ok <- !is.na(cvg)                                        # rows that report an interval
stopifnot(identical(ok, !is.na(wid)), identical(which(!ok), 1:2))

cols <- c(single = "grey55", full = "#D55E00", cross = "#0072B2", double = "#009E73")
mcol <- unname(cols[grp]); n <- length(lab); y <- n:1

## --- geometry (inches, final print size) ---
PS <- 8                         # base pointsize: cex = 1 prints at 8 pt
pt <- function(p) p / PS        # cex giving p points
pitch <- 0.165                  # row pitch (in)
ylim  <- c(0.45, n + 0.55)
ttl_h <- 0.20; top_h <- 0.43; plot_h <- diff(ylim) * pitch; axis_h <- 0.22
blk_gap <- 0.05; leg_min <- 0.21
## cairo_pdf truncates the page to whole points; request just above the needed height and lay
## everything out on the size the device actually reports.
H_req <- (ceiling((ttl_h + 2 * (top_h + plot_h + axis_h) + blk_gap + leg_min) * 72) + 0.25) / 72

f <- file.path(IMG, "ihdp_calibration.pdf")
cairo_pdf(f, width = 5.147, height = H_req, pointsize = PS, family = FONT)
W <- dev.size("in")[1]; H <- dev.size("in")[2]

leg_h <- H - (ttl_h + 2 * (top_h + plot_h + axis_h) + blk_gap)   # legend band takes the slack
lab_w <- 1.135                  # row-label column: 0.025 pad + 1.04 (longest label) + 0.07
val_w <- 0.29                   # coverage value column right of the coverage panel
gap   <- 0.10                   # white space before the right-hand panels
rmar  <- 0.03
Lw <- (W - lab_w - val_w - gap - rmar) / 2               # left and right panels share the width
lx <- c(lab_w, lab_w + Lw)
rx <- c(lx[2] + val_w + gap, W - rmar)
py2 <- c(leg_h + axis_h, leg_h + axis_h + plot_h)        # bottom row: coverage | width
py1 <- c(py2[2] + top_h + blk_gap + axis_h, 0)           # top row: bias | RMSE
py1[2] <- py1[1] + plot_h
caph <- 0.16                                             # MC-interval cap half-height (rows)
vtop <- n + 0.42                                         # reference lines stop just above row 1
pcex <- 1.0                                              # point size (radius ~3 pt)
in2x <- function(d) d * diff(par("usr")[1:2]) / par("pin")[1]   # inches -> user x
in2y <- function(d) d * diff(par("usr")[3:4]) / par("pin")[2]   # inches -> user y

ebar <- function(lo, hi, yy, col = "grey35") {
  segments(lo, yy, hi, yy, col = col, lwd = 1.1)
  segments(c(lo, hi), yy - caph, c(lo, hi), yy + caph, col = col, lwd = 1.1)
}
## text on a white backing (masks the row gridline under in-panel labels)
btext <- function(x, yy, s, adj = c(0, 0.5), cex = 1, font = 1, col = "grey25", pad = 0.015) {
  w <- strwidth(s, cex = cex, font = font); h <- strheight("Ay", cex = cex, font = font)
  x0 <- x - adj[1] * w
  rect(x0 - in2x(pad), yy - h / 2 - in2y(pad), x0 + w + in2x(pad), yy + h / 2 + in2y(pad),
       col = "white", border = NA)
  text(x, yy, s, adj = adj, cex = cex, font = font, col = col)
}
panel <- function(x01, y01, xlim) {   # open a panel whose plot region is x01 x y01 (inches)
  par(new = TRUE, plt = c(x01 / W, y01 / H))
  plot.new(); plot.window(xlim = xlim, ylim = ylim, xaxs = "i", yaxs = "i")
  segments(par("usr")[1], y, par("usr")[2], y, col = "grey93")
}
ax <- function(at, labels, minor = NULL) {
  axis(1, at = at, labels = labels, cex.axis = 1, lwd = 0.8, lwd.ticks = 0.8,
       tcl = -0.3, mgp = c(2, 0.32, 0), gap.axis = 0.25)
  if (length(minor)) axis(1, at = minor, labels = FALSE, lwd = 0.8, lwd.ticks = 0.8, tcl = -0.18)
}
rowlab <- function() text(par("usr")[1] - in2x(0.07), y, lab, adj = c(1, 0.5), col = mcol, font = 2, xpd = NA)
na_note <- function(xc) btext(xc, y[!ok], "not available (no interval)", adj = c(0.5, 0.5),
                              cex = pt(7.5), font = 3, col = "grey50")
mn <- function(v) sub("^-", "−", format(v, nsmall = 1))  # true minus sign

par(mai = c(0, 0, 0, 0), las = 1, cex = 1, xpd = FALSE)
plot.new()                                               # the page

## Mean bias
panel(lx, py1, c(-0.44, 0.18))
segments(0, ylim[1], 0, vtop, col = "grey35", lwd = 1)
segments(0, y, bias, y, col = mcol, lwd = 2.2)
ebar(bias - 1.96 * bias_mc, bias + 1.96 * bias_mc, y)
points(bias, y, pch = 19, col = mcol, cex = pcex)
rowlab()
ax(seq(-0.4, 0.1, 0.1), c(mn(-0.4), mn(-0.3), mn(-0.2), mn(-0.1), "0", "0.1"))

## RMSE
panel(rx, py1, c(0, 0.98))
segments(0, y, rmse, y, col = mcol, lwd = 2.2)
points(rmse, y, pch = 19, col = mcol, cex = pcex)
btext(rmse + in2x(0.07), y, sprintf("%.2f", rmse))
ax(seq(0, 0.8, 0.2), sprintf("%.1f", seq(0, 0.8, 0.2)))

## 95% CI coverage
panel(lx, py2, c(0.81, 1.015))
segments(0.95, ylim[1], 0.95, vtop, lty = 2, col = "grey35", lwd = 1)
segments(0.95, y[ok], cvg[ok], y[ok], col = mcol[ok], lwd = 2.2)
ebar(cvg[ok] - 1.96 * cvg_mc[ok], cvg[ok] + 1.96 * cvg_mc[ok], y[ok])
points(cvg[ok], y[ok], pch = 19, col = mcol[ok], cex = pcex)
text(par("usr")[2] + in2x(0.05), y[ok], sprintf("%.2f", cvg[ok]), adj = c(0, 0.5),
     col = "grey25", xpd = NA)
na_note(mean(par("usr")[1:2]))
rowlab()
ax(c(0.85, 0.90, 0.95, 1.00), c("0.85", "0.90", "0.95", "1.00"))

## Mean 95% CI width (descriptive)
panel(rx, py2, c(0, 1.9))
segments(0, y[ok], wid[ok], y[ok], col = mcol[ok], lwd = 2.2)
points(wid[ok], y[ok], pch = 19, col = mcol[ok], cex = pcex)
btext(wid[ok] + in2x(0.07), y[ok], sprintf("%.2f", wid[ok]))
na_note(mean(par("usr")[1:2]))
ax(seq(0, 1.5, 0.5), c("0", "0.5", "1.0", "1.5"))

## overlay in inches: overall title, panel titles/subtitles, legend
par(new = TRUE, plt = c(0, 1, 0, 1))
plot.new(); plot.window(xlim = c(0, W), ylim = c(0, H), xaxs = "i", yaxs = "i")
text(W / 2, H - 0.145, "IHDP calibration (50 realizations, n = 747; single split r = 1)",
     font = 2, cex = pt(9.5), adj = c(0.5, 0))
ttl <- function(x01, py, main, sub1, sub2 = NULL) {
  xc <- mean(x01); yt <- py[2] + 0.31; ys1 <- py[2] + 0.165; ys2 <- py[2] + 0.05
  text(xc, yt, main, font = 2, cex = pt(9), adj = c(0.5, 0))
  text(xc, ys1, sub1, cex = pt(7.5), col = "grey45", adj = c(0.5, 0))
  if (!is.null(sub2)) text(xc, ys2, sub2, cex = pt(7.5), col = "grey45", adj = c(0.5, 0))
}
ttl(lx, py1, "Mean bias", "(0 = unbiased;", "bars = 95% Monte Carlo CI)")
ttl(rx, py1, "RMSE", "(lower = more accurate)")
ttl(lx, py2, "95% CI coverage", "nominal 0.95 (dashed);", "bars = 95% Monte Carlo CI")
ttl(rx, py2, "Mean 95% CI width", "(descriptive mean;", "no uncertainty bars)")
lg <- c("singly robust","full-sample DR","cross-fit","double cross-fit")
legend(W / 2, leg_h / 2 + 0.01, xjust = 0.5, yjust = 0.5, horiz = TRUE, bty = "n",
       cex = 1, pt.cex = pcex, pch = 19, x.intersp = 0.6,
       text.width = strwidth(lg) + 0.09,      # per-entry widths (+0.09 in between entries)
       col = cols[c("single","full","cross","double")], legend = lg)
invisible(dev.off())
cat(sprintf("wrote %s  (%.3f x %.3f in = %g x %g bp)\n", f, W, H, round(W * 72), round(H * 72)))
