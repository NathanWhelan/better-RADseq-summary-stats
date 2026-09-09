#!/usr/bin/env Rscript
###############################################################################
#
#  het_between_pops.R -- command-line entry point for
#  RADdiversity::het_between_pops().
#
#  USAGE
#     Rscript het_between_pops.R <vcf> <popmap.tsv> [min_call] [outdir]
#     Rscript het_between_pops.R --selftest
#
#     min_call  minimum per-individual genotyping rate for a locus to be used
#               (default 0.9). outdir defaults to "." (the current directory).
#
#  This script is a thin wrapper: all the logic lives in the RADdiversity
#  package's het_between_pops() / het_between_pops_selftest() functions
#  (R/het_between_pops.R, R/selftest.R), which can also be called directly
#  from an R/RStudio session after library(RADdiversity).
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
         "`Rscript /path/to/het_between_pops.R ...`, or run ",
         "`R CMD INSTALL .` / `devtools::install()` from the package root ",
         "first).")
  }
  for (f in list.files(file.path(.d, "R"), pattern = "\\.R$", full.names = TRUE))
    source(f)
} else {
  library(RADdiversity)
}

.a0 <- commandArgs(trailingOnly = TRUE)
if (length(.a0) >= 1 && .a0[1] == "--selftest") {
  het_between_pops_selftest()
  quit(status = 0, save = "no")
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  cat("Usage: Rscript het_between_pops.R <vcf> <popmap.tsv> [min_call] [outdir]\n")
  cat("       Rscript het_between_pops.R --selftest\n\n")
  cat("  Tests whether heterozygosity differs between populations, using the\n")
  cat("  INDIVIDUAL as the unit of replication. Needs a VCF: per-individual\n")
  cat("  heterozygosity cannot be recovered from populations.sumstats.tsv.\n")
  quit(status = 1, save = "no")
}
vcf_file <- args[1]
popmap_f <- args[2]
min_call <- if (length(args) >= 3) args[3] else 0.9
outdir   <- if (length(args) >= 4) args[4] else "."

## Any error becomes a printed message and a non-zero exit status, since
## het_between_pops() itself must never call quit()/q() (that would kill an
## interactive RStudio session too).
res <- tryCatch(
  het_between_pops(vcf_file, popmap_f, min_call = min_call, outdir = outdir),
  error = function(e) {
    message("Error: ", conditionMessage(e))
    quit(status = 1, save = "no")
  }
)
