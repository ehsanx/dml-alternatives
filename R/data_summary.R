## =============================================================================
## data_summary.R -- sample sizes and design facts that the manuscript quotes, from the raw data.
## One row per dataset: n, treated units, the outcome mean (a proportion for a binary outcome) and
## the number of covariate columns the estimators see (RHC factors one-hot coded). For IHDP it also
## records whether the covariates and treatment are identical across all 50 outcome realizations.
## Needs the raw data under DML_DATA_DIR (see README.md); NHEFS comes from causaldata.
## Writes data_summary.csv to DML_OUT (default output/); the frozen copy is docs/data/data_summary.csv,
## which R/verify_numbers.R reads.  Run from the repository root:  Rscript R/data_summary.R
source("R/utils_data.R")
OUT <- Sys.getenv("DML_OUT", "output"); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

row_of <- function(d, name = d$name)
  data.frame(dataset = name, n = length(d$Y), n_treated = sum(d$A == 1), outcome_mean = mean(d$Y),
             n_covariate_columns = ncol(d$X), design_fixed_all_50 = NA, stringsAsFactors = FALSE)

ih <- lapply(1:50, load_ihdp)
same_design <- function(d) identical(d$A, ih[[1]]$A) && identical(unname(as.matrix(d$X)), unname(as.matrix(ih[[1]]$X)))
ihdp <- row_of(ih[[1]], "IHDP (realizations 1-50)")
ihdp$outcome_mean <- NA                        # the outcome differs by realization
ihdp$design_fixed_all_50 <- all(vapply(ih, same_design, TRUE))

tab <- rbind(row_of(load_nhefs()), row_of(load_rhc()), row_of(load_lalonde("psid")), ihdp)
write.csv(tab, file.path(OUT, "data_summary.csv"), row.names = FALSE)
print(tab, row.names = FALSE)
