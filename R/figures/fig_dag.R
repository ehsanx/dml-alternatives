## fig_dag.R -- Figure 1: the hypothesized DAG for the NHEFS analysis (nhefs_dag.pdf).
## A schematic: no data are read.
##   Rscript R/figures/fig_dag.R [output folder]          (default: output/figures)
##
## Drawn at its final printed size, so 1 pt in R is 1 pt on the page; type sizes are in printed points.
## cairo_pdf embeds the fonts. Device width: cairo_pdf truncates the page to whole PostScript points,
## so the figure is opened at (N + 0.5)/72 in and gets an N-bp page. N = 266 is the largest whole bp
## not exceeding the include width (0.72 x 372 TeX pt = 266.84 bp), so printed type is never smaller
## than the nominal size.

## ---- paths: the repository root is the nearest folder above that holds R/estimators.R ----
REPO <- getwd()
while (!file.exists(file.path(REPO, "R", "estimators.R")) && dirname(REPO) != REPO) REPO <- dirname(REPO)
DATA <- Sys.getenv("DML_OUT", file.path(REPO, "docs", "data"))             # frozen inputs, read only
args <- commandArgs(trailingOnly = TRUE)
IMG  <- if (length(args)) args[1] else file.path(REPO, "output", "figures")  # output folder
dir.create(IMG, showWarnings = FALSE, recursive = TRUE)
FONT <- Sys.getenv("DML_FONT", "Arial")   # without Arial, fontconfig substitutes a metric-compatible font

## ---- device and text helpers --------------------------------------------------
PS <- 9                                   # device pointsize; cex = pt / PS
## open a figure of wbp x hbp PostScript points; user coordinates = inches from bottom-left
open_fig <- function(file, wbp, hbp) {
  cairo_pdf(file, width = (wbp + 0.5) / 72, height = (hbp + 0.5) / 72,
            pointsize = PS, family = FONT)
  par(mar = c(0, 0, 0, 0), oma = c(0, 0, 0, 0), xpd = NA)
  plot.new()
  plot.window(xlim = c(0, wbp / 72), ylim = c(0, hbp / 72), xaxs = "i", yaxs = "i")
}
## plain (single-font) text with its baseline at y
btext <- function(x, y, s, pt, f = 1, adj = 0.5, col = "black")
  text(x, y, s, adj = c(adj, 0), cex = pt / PS, font = f, col = col)
cap_h <- function(pt) 0.716 * pt / 72     # Arial cap height (in)
arr <- function(x0, y0, x1, y1, col = "grey30", lwd = 1.3, len = 0.062)
  arrows(x0, y0, x1, y1, length = len, angle = 30, lwd = lwd, col = col)

## ---- Figure 1: X = the confounder set (named covariates listed), A = exposure, Y = outcome ----
local({
  wbp <- 266; hbp <- 147
  open_fig(file.path(IMG, "nhefs_dag.pdf"), wbp, hbp)
  W <- wbp / 72; H <- hbp / 72
  P <- list(node = 15, label = 9, cov = 8)       # printed pt
  labA <- "quit smoking"; labY <- "weight change (kg)"
  wA <- strwidth(labA, "user", cex = P$label / PS); wY <- strwidth(labY, "user", cex = P$label / PS)
  ## symmetric triangle (A and Y equidistant from X), shifted so the left and right
  ## margins of the node captions are equal (0.1 in)
  Xx <- W / 2 + (wA - wY) / 4; d <- Xx - wA / 2 - 0.10
  Ax <- Xx - d; Yx <- Xx + d
  ## X node: confounder set with the named covariates (wording = Figure 1 caption)
  cov <- c("age, sex, race, education, baseline weight,",
           "smoking intensity and duration,",
           "and physical activity")
  yX <- H - 0.035 - cap_h(P$node)                # baseline of the X letter
  btext(Xx, yX, "X", P$node, f = 2)
  yBC <- yX - 0.205
  btext(Xx, yBC, "baseline confounders", P$label, f = 3)
  yc <- yBC - 0.165 - 0.137 * (seq_along(cov) - 1)
  for (i in seq_along(cov)) btext(Xx, yc[i], cov[i], P$cov, col = "grey30")
  ## A and Y nodes
  yL <- 0.31                                      # baseline of the node letters
  yC <- yL - 0.205                                # baseline of the node captions
  btext(Ax, yL, "A", P$node, f = 2); btext(Ax, yC, labA, P$label)
  btext(Yx, yL, "Y", P$node, f = 2); btext(Yx, yC, labY, P$label)
  ## edges, with gaps so arrowheads clear the node labels
  ytop <- min(yc) - 0.075                         # just below the covariate block
  yhead <- yL + cap_h(P$node) + 0.07              # just above the node letters
  arr(Xx - 0.52, ytop, Ax + 0.13, yhead, lwd = 1.5, len = 0.07)        # X -> A
  arr(Xx + 0.52, ytop, Yx - 0.13, yhead, lwd = 1.5, len = 0.07)        # X -> Y
  yAY <- yC + 0.045
  arr(Ax + wA / 2 + 0.10, yAY, Yx - wY / 2 - 0.10, yAY, lwd = 1.5, len = 0.07)   # A -> Y
  dev.off()
})
cat("wrote", file.path(IMG, "nhefs_dag.pdf"), "\n")
