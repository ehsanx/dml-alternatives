## =============================================================================
## check_data.R -- compare the raw data files with data/checksums.csv before any refit.
## Reads the files under DML_DATA_DIR (default data/) and reports, for each one, whether it is
## missing or has the expected md5 (listed with its expected rows and columns). NHEFS needs no file
## (causaldata package). Exit status 1 if any file is missing or differs.
## Run from the repository root:  Rscript R/check_data.R
DATA_DIR <- Sys.getenv("DML_DATA_DIR", "data")
ck <- read.csv(file.path("data", "checksums.csv"), stringsAsFactors = FALSE)
ck$status <- vapply(seq_len(nrow(ck)), function(i) {
  f <- file.path(DATA_DIR, ck$file[i])
  if (!file.exists(f)) "missing" else if (unname(tools::md5sum(f)) == ck$md5[i]) "ok" else "md5 differs"
}, "")
print(ck[, c("dataset", "file", "rows", "columns", "status")], row.names = FALSE)
bad <- sum(ck$status != "ok")
cat(sprintf("\n%d of %d files match data/checksums.csv\n", nrow(ck) - bad, nrow(ck)))
if (bad && !interactive()) quit(status = 1)
