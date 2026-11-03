## make_figures.R -- draw all 10 manuscript figures as PDFs.
##   Rscript R/make_figures.R [output folder]          (default: output/figures)
## Figures 1, 2 and 5-8 are drawn by the R scripts in R/figures/; Figures 3, 4, 9 and 10 by
## python/paper_figures.py (packages in requirements.txt). Set PYTHON to choose the interpreter.
## All inputs are the frozen CSVs in docs/data/; DML_OUT points the scripts to another folder.
## Compare the result with the manuscript's copies:
##   python tools/compare_pdf.py output/figures figures/reference

## the repository root: two levels above this script
self <- sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE))
REPO <- if (length(self)) normalizePath(file.path(dirname(self[1]), "..")) else getwd()
stopifnot(file.exists(file.path(REPO, "R", "estimators.R")))
setwd(REPO)
args <- commandArgs(trailingOnly = TRUE)
IMG <- normalizePath(if (length(args)) args[1] else file.path(REPO, "output", "figures"), mustWork = FALSE)
dir.create(IMG, showWarnings = FALSE, recursive = TRUE)

rscript <- file.path(R.home("bin"), "Rscript")
for (s in c("fig_dag", "fig_overlap", "fig_forests", "fig_ihdp", "fig_reps")) {
  status <- system2(rscript, c(shQuote(file.path("R", "figures", paste0(s, ".R"))), shQuote(IMG)))
  if (status != 0) stop(s, ".R failed")
}

## Python: $PYTHON if set; otherwise the first of python3 and python whose matplotlib is the version
## pinned in requirements.txt (another version draws the same figures, but not byte for byte)
want <- sub("^matplotlib==", "", grep("^matplotlib==", readLines("requirements.txt"), value = TRUE))
mpl <- function(py) if (!nzchar(py)) NA_character_ else suppressWarnings(tryCatch(
  system2(py, c("-c", shQuote("import matplotlib; print(matplotlib.__version__)")), stdout = TRUE, stderr = FALSE)[1],
  error = function(e) NA_character_))
cand <- if (nzchar(Sys.getenv("PYTHON"))) Sys.getenv("PYTHON") else unname(Sys.which(c("python3", "python")))
ver <- vapply(cand, mpl, "")
ok <- which(!is.na(ver))
pick <- if (nzchar(Sys.getenv("PYTHON"))) ok else c(ok[ver[ok] == want], ok)
if (length(pick)) {
  py <- cand[pick[1]]
  if (ver[pick[1]] != want)
    message("matplotlib ", ver[pick[1]], " in ", py, "; requirements.txt pins ", want,
            ": Figures 3, 4, 9 and 10 may not match figures/reference byte for byte")
  if (system2(py, c(shQuote(file.path("python", "paper_figures.py")), shQuote(IMG))) != 0)
    stop("python/paper_figures.py failed")
} else {
  cat("No Python with matplotlib found. Draw Figures 3, 4, 9 and 10 with:\n",
      " python python/paper_figures.py", IMG, "\n")
}
cat(length(list.files(IMG, "\\.pdf$")), "figure PDFs in", IMG, "\n")
