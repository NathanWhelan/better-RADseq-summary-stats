#!/usr/bin/env Rscript
###############################################################################
#
#  diversity_core.R -- command-line entry point for
#  RADdiversity::diversity_core_selftest().
#
#  USAGE
#     Rscript diversity_core.R --selftest
#
#  Checks every estimator in R/estimators.R against brute-force Monte Carlo,
#  and the Stacks formulas against real published Stacks output. Needs no
#  data and no packages beyond RADdiversity itself.
#
#  This script is a thin wrapper: the estimator functions and the self-test
#  live in the RADdiversity package (R/estimators.R, R/selftest.R), and can
#  also be called directly from an R/RStudio session after
#  library(RADdiversity).
#
###############################################################################

## ---------------------------------------------------------------------------
## Load the package: installed copy if available, else source R/*.R straight
## from this script's own directory (so `bash test/run_tests.sh` and any
## other dev workflow keep working on a fresh clone with nothing installed).
## ---------------------------------------------------------------------------
.self <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
.d <- if (is.na(.self)) NA_character_ else dirname(.self)
if (!requireNamespace("RADdiversity", quietly = TRUE)) {
  if (is.na(.d) || !file.exists(file.path(.d, "R"))) {
    stop("RADdiversity is not installed, and this script's own directory ",
         "could not be resolved to fall back to its R/ sources (invoke via ",
         "`Rscript /path/to/diversity_core.R ...`, or run ",
         "`R CMD INSTALL .` / `devtools::install()` from the package root ",
         "first).")
  }
  for (f in list.files(file.path(.d, "R"), pattern = "\\.R$", full.names = TRUE))
    source(f)
} else {
  library(RADdiversity)
}

if (identical(commandArgs(trailingOnly = TRUE)[1], "--selftest")) {
  ok <- diversity_core_selftest()
  quit(status = if (ok) 0 else 1, save = "no")
}
