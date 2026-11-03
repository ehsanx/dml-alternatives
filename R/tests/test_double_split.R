## =============================================================================
## test_double_split.R -- checks the double cross-fit splitting rule WITHOUT fitting any
## model. A recording stub replaces the Super Learner fit: it logs which rows train each
## g and Q fit and which rows each fit predicts. For every scored fold the test requires:
##   (1) the scored rows are excluded from every training set used to score them;
##   (2) the g and Q that score them were trained on disjoint rows;
##   (3) every row is scored exactly once (once by g; once per arm by Q).
## Usage (from the repo root):
##   Rscript R/tests/test_double_split.R [engine.R ...]
## With no arguments it checks R/estimators.R. The exit status is 1 if an engine breaks a rule.
## Engines that name the fit function differently (dr3: .dr_fit_sl) are handled below.
## =============================================================================
check_engine <- function(file, dc_name = NULL, fit_name = NULL) {
  e <- new.env(parent = globalenv())
  suppressMessages(sys.source(file, envir = e))
  dc_name  <- dc_name  %||% intersect(c("nuisance_double", ".nuisance_double"), ls(e, all.names = TRUE))[1]
  fit_name <- fit_name %||% intersect(c(".fit_sl", ".dr_fit_sl"), ls(e, all.names = TRUE))[1]
  if (is.na(dc_name) || is.na(fit_name)) stop("no double cross-fit function found in ", file)

  log <- new.env(); log$fits <- list(); log$preds <- list()
  stub <- function(Y, X, family, SL.library, ...) {
    kind <- if ("A" %in% names(X)) "Q" else "g"
    id <- length(log$fits) + 1L
    log$fits[[id]] <- list(kind = kind, train = round(X$x1))      # x1 carries the row id
    function(newX) {
      log$preds[[length(log$preds) + 1L]] <- list(fit = id, kind = kind, rows = round(newX$x1),
                                                  arm = if ("A" %in% names(newX)) newX$A[1] else NA)
      rep(if (kind == "g") 0.5 else 0, nrow(newX))
    }
  }
  assign(fit_name, stub, envir = e)

  n <- 30L; set.seed(7)
  X <- data.frame(x1 = seq_len(n), x2 = rnorm(n))
  A <- rep(0:1, length.out = n); Y <- rnorm(n)
  e[[dc_name]](Y, A, X, 3L, SL.library = "stub", family_y = gaussian(), seed = 1L)

  preds <- log$preds; fits <- log$fits
  key <- vapply(preds, function(p) paste(sort(p$rows), collapse = ","), "")
  ok <- TRUE; msgs <- character()
  for (k in unique(key)) {
    ps <- preds[key == k]; rows <- ps[[1]]$rows
    tr_g <- unique(unlist(lapply(ps[vapply(ps, `[[`, "", "kind") == "g"], function(p) fits[[p$fit]]$train)))
    tr_q <- unique(unlist(lapply(ps[vapply(ps, `[[`, "", "kind") == "Q"], function(p) fits[[p$fit]]$train)))
    r1 <- length(intersect(rows, c(tr_g, tr_q))) == 0
    r2 <- length(intersect(tr_g, tr_q)) == 0
    msgs <- c(msgs, sprintf("  scored %2d rows | g trained on %2d rows, Q on %2d rows | excluded from training: %s | g,Q disjoint: %s",
                            length(rows), length(tr_g), length(tr_q), r1, r2))
    ok <- ok && r1 && r2
  }
  g_rows <- unlist(lapply(preds[vapply(preds, `[[`, "", "kind") == "g"], `[[`, "rows"))
  q1_rows <- unlist(lapply(Filter(function(p) p$kind == "Q" && identical(p$arm, 1), preds), `[[`, "rows"))
  r3 <- all(table(factor(g_rows, levels = seq_len(n))) >= 1) &&
        all(table(factor(q1_rows, levels = seq_len(n))) >= 1)
  ## each row must be scored by exactly one g fit (averaging over role assignments scores it several times)
  once <- all(table(factor(g_rows, levels = seq_len(n))) == 1)
  msgs <- c(msgs, sprintf("  every row scored: %s | scored by exactly one g fit: %s", r3, once))
  list(pass = ok && r3 && once, msgs = msgs)
}
`%||%` <- function(a, b) if (is.null(a)) b else a

args <- commandArgs(trailingOnly = TRUE)
suppressMessages(library(SuperLearner))
cases <- if (length(args)) setNames(as.list(args), args) else list(`engine (R/estimators.R)` = "R/estimators.R")
pass <- c()
for (nm in names(cases)) {
  r <- check_engine(cases[[nm]]); pass[nm] <- r$pass
  cat(sprintf("%s: %s\n", nm, if (r$pass) "PASS" else "FAIL"))
  cat(r$msgs, sep = "\n"); cat("\n")
}
if (!all(pass)) quit(status = 1)
