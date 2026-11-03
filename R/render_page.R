## =============================================================================
## render_page.R -- build the companion page, docs/index.html, from the frozen CSVs.
##   Rscript R/render_page.R
## 1. Draws the 10 manuscript figures into output/figures/ (R/make_figures.R: the R scripts in
##    R/figures/, then python/paper_figures.py), and writes the numeric tables into output/tables/
##    (R/make_tables.R).
## 2. Compares them with the manuscript's copies in figures/reference/ (tools/compare_pdf.py) and keeps
##    the report in output/figures/compare.txt, which the page prints.
## 3. Rasterizes each PDF at 200 dpi into output/figures/png/. The page shows these images, so it shows
##    exactly the PDFs, not a redraw.
## 4. Renders docs/index.qmd with Quarto (1.7.31) into the self-contained docs/index.html.
## No model is fitted; the page reads docs/data/ only.
self <- sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE))
REPO <- if (length(self)) normalizePath(file.path(dirname(self[1]), "..")) else getwd()
stopifnot(file.exists(file.path(REPO, "R", "estimators.R")))
setwd(REPO)
Sys.unsetenv("DML_OUT")            # the page is always built from the frozen CSVs in docs/data/
IMG <- file.path("output", "figures")
unlink(IMG, recursive = TRUE)

rscript <- file.path(R.home("bin"), "Rscript")
if (system2(rscript, c(file.path("R", "make_figures.R"), IMG)) != 0) stop("R/make_figures.R failed")
if (system2(rscript, file.path("R", "make_tables.R")) != 0) stop("R/make_tables.R failed")
pdfs <- list.files(IMG, "\\.pdf$", full.names = TRUE)
if (length(pdfs) != 10) stop("expected 10 figure PDFs in ", IMG, ", found ", length(pdfs))

## the comparison report, from the first Python that has the packages in requirements.txt
cmp <- "Python not found: the figures were not compared with figures/reference/"
for (py in Filter(nzchar, c(Sys.getenv("PYTHON"), unname(Sys.which(c("python3", "python")))))) {
  out <- suppressWarnings(tryCatch(system2(py, c(file.path("tools", "compare_pdf.py"), IMG,
    file.path("figures", "reference")), stdout = TRUE, stderr = TRUE), error = function(e) NULL))
  if (any(grepl("of [0-9]+ identical$", out))) { cmp <- out; break }
}
writeLines(cmp, file.path(IMG, "compare.txt"))
cat(cmp, sep = "\n")

dir.create(file.path(IMG, "png"), showWarnings = FALSE)
for (f in pdfs)   # a fixed file name, not pdf_convert's sprintf template, so silence its "arguments not used" warning
  suppressWarnings(pdftools::pdf_convert(f, format = "png", pages = 1, dpi = 200, verbose = FALSE,
                                         filenames = file.path(IMG, "png", sub("\\.pdf$", ".png", basename(f)))))

if (system2("quarto", c("render", file.path("docs", "index.qmd"))) != 0) stop("quarto render failed")
