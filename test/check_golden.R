#!/usr/bin/env Rscript
###############################################################################
#
#  check_golden.R -- numeric regression check for run_tests.sh.
#
#  The rest of the suite only asserts exit code (ok/err). That catches a
#  script crashing but not a script running to completion and printing a
#  WRONG number -- a swapped column, a shifted index, a broken join would
#  all still exit 0. This compares a script's TSV output against a
#  checked-in "golden" copy, column-for-column, numeric columns within a
#  small tolerance and everything else exactly.
#
#  USAGE
#     Rscript check_golden.R <actual.tsv> <golden.tsv> [tol]
#     Rscript check_golden.R <actual.tsv> <golden.tsv> --update
#
#     tol       numeric tolerance (default 1e-6). Output columns are already
#               rounded to 4 decimal places by the scripts that produce them,
#               so a real regression differs by orders of magnitude more than
#               this; it exists only to absorb text<->double round-tripping.
#     --update  overwrite golden.tsv with actual.tsv instead of comparing.
#               Use this ONLY after confirming by hand that a change in the
#               numbers is intentional (a formula fix, a new column) and not
#               a regression -- never to make a failing check pass unseen.
#
#  Exits 0 and prints one OK line on a match; exits 1 and prints every
#  differing cell (column, row, golden value, actual value) otherwise.
#
###############################################################################

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  cat("Usage: Rscript check_golden.R <actual.tsv> <golden.tsv> [tol | --update]\n")
  q(status = 1)
}
actual_f <- args[1]; golden_f <- args[2]
update   <- length(args) >= 3 && args[3] == "--update"
tol      <- if (length(args) >= 3 && !update) suppressWarnings(as.numeric(args[3])) else 1e-6
if (is.na(tol)) tol <- 1e-6

if (!file.exists(actual_f)) stop("Actual output not found: ", actual_f,
  "\n  The script that should have produced it may have failed upstream.")

if (update) {
  dir.create(dirname(golden_f), showWarnings = FALSE, recursive = TRUE)
  file.copy(actual_f, golden_f, overwrite = TRUE)
  cat("UPDATED golden: ", golden_f, " <- ", actual_f, "\n", sep = "")
  q(status = 0)
}

if (!file.exists(golden_f))
  stop("No golden file at ", golden_f, ". Create one deliberately with:\n",
       "  Rscript check_golden.R ", actual_f, " ", golden_f, " --update")

a <- read.delim(actual_f, check.names = FALSE, stringsAsFactors = FALSE)
g <- read.delim(golden_f, check.names = FALSE, stringsAsFactors = FALSE)

problems <- character(0)

if (!identical(names(a), names(g)))
  problems <- c(problems, sprintf(
    "column names/order differ:\n    golden: %s\n    actual: %s",
    paste(names(g), collapse = ", "), paste(names(a), collapse = ", ")))

if (nrow(a) != nrow(g))
  problems <- c(problems, sprintf("row count differs: golden %d, actual %d",
                                  nrow(g), nrow(a)))

## Compare whatever columns and rows both files have, even after a mismatch
## above, so one wrong column doesn't hide every other difference.
common_cols <- intersect(names(g), names(a))
n <- min(nrow(g), nrow(a))
for (col in common_cols) {
  gv <- g[[col]][seq_len(n)]; av <- a[[col]][seq_len(n)]
  if (is.numeric(gv) || is.numeric(av)) {
    gv <- suppressWarnings(as.numeric(gv)); av <- suppressWarnings(as.numeric(av))
    bad <- which(!(is.na(gv) & is.na(av)) & (is.na(gv) != is.na(av) |
                  (!is.na(gv) & abs(gv - av) > tol)))
  } else {
    bad <- which(as.character(gv) != as.character(av))
  }
  for (i in bad)
    problems <- c(problems, sprintf("  [row %d] %s: golden = %s, actual = %s",
                                    i, col, format(gv[i]), format(av[i])))
}

if (length(problems)) {
  cat("MISMATCH: ", actual_f, " vs ", golden_f, "\n", sep = "")
  cat(paste0("  ", problems, collapse = "\n"), "\n", sep = "")
  cat("  If this difference is an intentional, verified change (not a\n")
  cat("  regression), update the golden file with:\n")
  cat("    Rscript check_golden.R ", actual_f, " ", golden_f, " --update\n", sep = "")
  q(status = 1)
} else {
  cat(sprintf("OK: %s matches %s (%d cols, %d rows)\n",
              actual_f, golden_f, length(common_cols), n))
  q(status = 0)
}
