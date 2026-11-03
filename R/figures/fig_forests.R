## fig_forests.R -- Figure 5 (nhefs_forest.pdf) and Figure 6 (rhc_forest.pdf): the NHEFS and RHC
## estimator forest plots.
## Reads docs/data/nhefs_panel_table.csv (Figure 5), real_application.csv and rhc_reps_table.csv
## (Figure 6).
##   Rscript R/figures/fig_forests.R [output folder]      (default: output/figures)
##
## Drawn at the final printed size (width 5.147 in, scale 1) with cairo_pdf, so fonts are embedded and
## cex x 8 pt is the printed size. The text columns are sized from the measured label widths.
## Vertical guides: NHEFS, the Hernan-Robins band and the two dotted r=100 span lines; RHC, the no-effect
## line and the r=100 span band.

## ---- paths: the repository root is the nearest folder above that holds R/estimators.R ----
REPO <- getwd()
while (!file.exists(file.path(REPO, "R", "estimators.R")) && dirname(REPO) != REPO) REPO <- dirname(REPO)
DATA <- Sys.getenv("DML_OUT", file.path(REPO, "docs", "data"))             # frozen inputs, read only
args <- commandArgs(trailingOnly = TRUE)
IMG  <- if (length(args)) args[1] else file.path(REPO, "output", "figures")  # output folder
dir.create(IMG, showWarnings = FALSE, recursive = TRUE)
FONT <- Sys.getenv("DML_FONT", "Arial")   # without Arial, fontconfig substitutes a metric-compatible font

C_FULL <- "#444444"; C_R1 <- "#D55E00"; C_R100 <- "#0072B2"   # Okabe-Ito orange / blue
C_PLUGIN_LAB <- "grey35"   # g-computation / IPTW row labels

## ---- data: NHEFS -------------------------------------------------------------------
p <- read.csv(file.path(DATA, "nhefs_panel_table.csv"), stringsAsFactors = FALSE)
gf <- function(e)      p[p$estimator == e & is.na(p$r), ]
gs <- function(e, rr)  p[p$estimator == e & !is.na(p$r) & p$r == rr, ]

## ---- data: RHC ---------------------------------------------------------------------
f <- read.csv(file.path(DATA, "real_application.csv"), stringsAsFactors = FALSE)
f <- f[grepl("^RHC", f$dataset), ]
gfr <- function(e) f[f$estimator == e, ]
r  <- read.csv(file.path(DATA, "rhc_reps_table.csv"), stringsAsFactors = FALSE)
gr <- function(m, rr) r[r$method == m & r$r == rr, ]

one <- function(d, what) { if (nrow(d) != 1) stop("expected one row for ", what, ", got ", nrow(d)); d }

## rows: y positions, labels, colours and plotting symbols
build_rows <- function(full, split, full_keys, split_keys) {
  rows <- list()
  add <- function(y, label, est, lo, hi, col, pch, hasci)
    rows[[length(rows) + 1]] <<- list(y = y, label = label, est = est, lo = lo, hi = hi,
                                      col = col, pch = pch, hasci = hasci)
  add(15.0, "g-computation", one(full(full_keys[1]), full_keys[1])$estimate, NA, NA, C_FULL, 1, FALSE)
  add(14.0, "IPTW",          one(full(full_keys[2]), full_keys[2])$estimate, NA, NA, C_FULL, 1, FALSE)
  x <- one(full(full_keys[3]), full_keys[3]); add(13.0, "AIPW (one-step)", x$estimate, x$ci_lo, x$ci_hi, C_FULL, 18, TRUE)
  x <- one(full(full_keys[4]), full_keys[4]); add(12.0, "TMLE",            x$estimate, x$ci_lo, x$ci_hi, C_FULL, 18, TRUE)
  pair <- function(y, lab, m) {
    a <- one(split(m, 1), m);   add(y,       sprintf("%s — r=1",   lab), a$estimate, a$ci_lo, a$ci_hi, C_R1,   16, TRUE)
    b <- one(split(m, 100), m); add(y - 0.7, sprintf("%s — r=100", lab), b$estimate, b$ci_lo, b$ci_hi, C_R100, 15, TRUE)
  }
  pair(9.6, "DML (cross-fit AIPW)", split_keys[1])
  pair(7.9, "CV-TMLE",              split_keys[2])
  pair(6.2, "DC-AIPW",              split_keys[3])
  pair(4.5, "DC-TMLE",              split_keys[4])
  rows
}
NH_SPLIT  <- c("DML (cross-fit AIPW)", "CV-TMLE", "DC-AIPW", "DC-TMLE")
RHC_SPLIT <- c("DML (cross-fit)", "CV-TMLE (cross-fit)", "DC-AIPW", "DC-TMLE")
rows_nh  <- build_rows(gf,  gs, c("g-computation", "IPTW", "AIPW (one-step)", "TMLE"), NH_SPLIT)
rows_rhc <- build_rows(gfr, gr, c("g-computation", "IPTW", "AIPW", "TMLE"),            RHC_SPLIT)
dr100_nh  <- sapply(NH_SPLIT,  function(m) gs(m, 100)$estimate)   # span of the four r=100 DR estimates
dr100_rhc <- sapply(RHC_SPLIT, function(m) gr(m, 100)$estimate)

fmt_nh  <- function(z) if (z$hasci) sprintf("%.2f  [%.2f, %.2f]", z$est, z$lo, z$hi) else sprintf("%.2f  (point est.)", z$est)
fmt_rhc <- function(z) if (z$hasci) sprintf("%.3f  [%.3f, %.3f]", z$est, z$lo, z$hi) else sprintf("%.3f  (point est.)", z$est)

HDR_FULL  <- "Full-sample nuisance (single fit)"
HDR_SPLIT <- "Repeated cross-fitting:  r = 1 (single split)  vs  r = 100 (median)"
HDR_EST   <- "Estimate [95% CI]"

## ---- final-size geometry ----------------------------------------------------------
FIG_W    <- 5.147        # in; = \linewidth (372 TeX pt) -> include at width=\linewidth, scale 1
                         # (cairo_pdf truncates to 370 bp = 5.139 in: include scale 1.0016, 8 pt -> 8.01 pt)
PS       <- 8            # pt; cex 1 = 8 pt in print
CEX_ANN  <- 8 / PS       # benchmark-band label and "no effect" (grey40: secondary by colour)
CEX_XLAB <- 9 / PS       # axis title
U_PT     <- 18           # printed height of one y unit (rows 1 unit apart; r=1/r=100 pairs 0.7 = 12.6 pt)
PAD_IN   <- 2 / 72       # outer padding, left and right
OFF_IN   <- 4 / 72       # gap between the plot region and the label / estimate columns
Y_T      <- 16.65        # column-header row: "Estimate [95% CI]" and the guide labels
Y_H1     <- 15.75        # "Full-sample nuisance (single fit)"
Y_H2     <- 10.55        # "Repeated cross-fitting: ..."
BLOCKS   <- list(c(11.60, 15.40), c(3.40, 10.00))   # data-row extents of the two blocks (guides)
Y_BOT    <- 3.00         # x axis
ANN_DY   <- 9.6 / U_PT   # line step of the two-line benchmark label (1.2 x 8 pt)
MAI_BOT  <- 0.42         # in; ticks, tick labels and axis title
MAI_TOP  <- 2 / 72       # in
TCL      <- -0.35        # tick length, lines
MGP_X    <- 0.40         # tick-label line
XLAB_LINE <- 1.75        # axis-title line
YADJ     <- 1 / 3        # vertical text adjustment R uses for pos = 2 / 4 (yCharOffset)
LWD_CI <- 1.7; LWD_WH <- 1.3; LWD_PT <- 1.3                   # 1.3 / 1.0 / 1.0 pt
CEX_DIAMOND <- 1.45; CEX_PT <- 1.10                           # pch 18; pch 1 / 15 / 16
WH <- 0.12                                                    # whisker half-height, y units

## text-column widths measured with the output device's own font metrics.
## Measured on the first output file, which draw_forest() then overwrites, so nothing else is written.
measure <- function(file) {
  cairo_pdf(file, width = FIG_W, height = 4, pointsize = PS, family = FONT)
  on.exit(dev.off())
  par(mar = c(0, 0, 0, 0)); plot.new()
  w <- function(s, font = 1) max(strwidth(s, units = "inches", font = font))
  c(lab = w(sapply(c(rows_nh, rows_rhc), `[[`, "label")),
    est = max(w(c(sapply(rows_nh, fmt_nh), sapply(rows_rhc, fmt_rhc))), w(HDR_EST, font = 2)))
}
TW <- measure(file.path(IMG, "nhefs_forest.pdf"))
MAI_L <- PAD_IN + TW[["lab"]] + OFF_IN       # same plot region in both figures
MAI_R <- OFF_IN + TW[["est"]] + PAD_IN

## ---- one forest panel -------------------------------------------------------------
draw_forest <- function(file, rows, fmt, xlim, ticks, xlab, y_top, guides, guide_labels) {
  fig_h <- MAI_BOT + MAI_TOP + (y_top - Y_BOT) * U_PT / 72
  cairo_pdf(file, width = FIG_W, height = fig_h, pointsize = PS, family = FONT)
  par(mai = c(MAI_BOT, MAI_L, MAI_TOP, MAI_R), las = 1, lend = "butt", xpd = NA)
  plot(NA, xlim = xlim, ylim = c(Y_BOT, y_top), yaxs = "i", xlab = "", ylab = "", axes = FALSE)
  usr <- par("usr")
  dx  <- OFF_IN / par("pin")[1] * diff(usr[1:2])           # OFF_IN in x user units
  x_left <- grconvertX(PAD_IN, from = "inches", to = "user")  # flush-left edge for the headers

  for (b in BLOCKS) guides(b)       # guides over the data rows of each block only
  guide_labels()                    # their labels, in the column-header row

  text(x_left, Y_H1, HDR_FULL,  adj = c(0, YADJ), font = 3, col = "grey25")
  text(x_left, Y_H2, HDR_SPLIT, adj = c(0, YADJ), font = 3, col = "grey25")
  text(usr[2] + dx, Y_T, HDR_EST, adj = c(0, YADJ), font = 2)

  for (z in rows) {
    if (z$hasci) {
      segments(z$lo, z$y, z$hi, z$y, col = z$col, lwd = LWD_CI)
      segments(c(z$lo, z$hi), z$y - WH, c(z$lo, z$hi), z$y + WH, col = z$col, lwd = LWD_WH)
    }
    points(z$est, z$y, pch = z$pch, col = z$col, bg = z$col,
           cex = if (z$pch == 18) CEX_DIAMOND else CEX_PT, lwd = LWD_PT)
    text(usr[1] - dx, z$y, z$label, adj = c(1, YADJ),
         col = if (z$pch == 1) C_PLUGIN_LAB else "black")
    text(usr[2] + dx, z$y, fmt(z), adj = c(0, YADJ), col = z$col)
  }
  axis(1, at = ticks, pos = Y_BOT, tcl = TCL, mgp = c(0, MGP_X, 0))
  mtext(xlab, side = 1, line = XLAB_LINE, cex = CEX_XLAB)
  dsz <- dev.size("in")             # cairo_pdf truncates the page to whole points
  invisible(dev.off())
  cat(sprintf("wrote %s  (%.3f x %.3f in = %.0f x %.0f pt; plot region %.2f in wide; label column %.0f pt, estimate column %.0f pt)\n",
              file, dsz[1], dsz[2], 72 * dsz[1], 72 * dsz[2], dsz[1] - MAI_L - MAI_R,
              72 * TW[["lab"]], 72 * TW[["est"]]))
}

## ---- Figure 5: NHEFS ------------------------------------------------------------------
draw_forest(
  file.path(IMG, "nhefs_forest.pdf"), rows_nh, fmt_nh,
  xlim = c(2.30, 4.55), ticks = seq(2.5, 4.5, 0.5),
  xlab = "average treatment effect  (weight change, kg)",
  y_top = Y_T + ANN_DY + 0.32,
  guides = function(b) {
    ## Hernan-Robins benchmark band (~3.4-3.5 kg)
    rect(3.4, b[1], 3.5, b[2], col = "#88888822", border = NA)
    ## consensus band: span of the four r=100 doubly-robust split estimates
    segments(range(dr100_nh), b[1], range(dr100_nh), b[2], lty = 3, col = "#0072B2", lwd = 1)
  },
  guide_labels = function() {
    text(3.45, Y_T + ANN_DY, "Hernán–Robins", adj = c(0.5, YADJ), cex = CEX_ANN, col = "grey40")
    text(3.45, Y_T,          "3.4–3.5 kg",          adj = c(0.5, YADJ), cex = CEX_ANN, col = "grey40")
  })

## ---- Figure 6: RHC ----------------------------------------------------------------------
draw_forest(
  file.path(IMG, "rhc_forest.pdf"), rows_rhc, fmt_rhc,
  xlim = c(-0.006, 0.086), ticks = seq(0, 0.08, 0.02),
  xlab = "risk difference   (ATE of RHC on 30-day mortality)",
  y_top = Y_T + 0.32,
  guides = function(b) {
    rect(min(dr100_rhc), b[1], max(dr100_rhc), b[2], col = "#0072B215", border = NA)
    segments(0, b[1], 0, b[2], lty = 2, col = "grey40")
  },
  guide_labels = function()
    text(0, Y_T, "no effect", adj = c(0.5, YADJ), cex = CEX_ANN, col = "grey40"))
