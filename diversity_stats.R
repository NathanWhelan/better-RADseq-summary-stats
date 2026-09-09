#!/usr/bin/env Rscript
###############################################################################
#
#  diversity_stats.R -- command-line entry point for RADdiversity::diversity_stats().
#
#  USAGE
#     Rscript diversity_stats.R <vcf> <popmap.tsv> --g=N [--nboot=N] [--boot=MODE] [--sites=N] [--min-n=N]
#     Rscript diversity_stats.R <vcf> <popmap.tsv> --g=N ... --complete-case
#
#     --g=N     REQUIRED. Rarefaction size in GENE COPIES (10 diploids = 20).
#               No default -- must be <= 2x the smallest population. Choose
#               this deliberately rather than let it silently track whichever
#               population happens to be smallest.
#     --nboot=N bootstrap replicates (default 10000; 0 = none)
#     --boot=MODE  which axis the bootstrap resamples: `loci` (RAD loci,
#               the DEFAULT -- matches hierfstat), `individuals` (individuals
#               within each population, matches diveRsity), or `both` (both
#               at once, the Owen & Eckles 2012 crossed-factor scheme).
#               `individuals`/`both` are COMPARISON modes only -- see
#               README.md, "Bootstrap mode".
#     --sites=N total sequenced sites -- the `Sites` column of the "All
#               positions (variant and fixed)" block of
#               populations.sumstats_summary.tsv. Adds per-sequenced-site Ho
#               and He. Omit for per-record values only.
#     --min-n=N DEFAULT (available-data) mode only: minimum typed individuals
#               a population needs at a locus to use that locus for that
#               population. Default 2.
#     --complete-case   use the OLD, stricter rule instead: a record is used
#               only if EVERY individual of EVERY population is genotyped
#               there (Schmidt et al. 2021 recommendation (b)).
#
#  Run it TWICE. populations.snps.vcf gives Ho, He, pct_poly and the
#  per-sequenced-site values; populations.haps.vcf gives FIS, Ar and privAr.
#
#  This script is a thin wrapper: all the logic lives in the RADdiversity
#  package's diversity_stats() function (R/diversity_stats.R), which can also
#  be called directly from an R/RStudio session after library(RADdiversity).
#
###############################################################################

print_usage <- function() {
  cat("Usage: Rscript diversity_stats.R <vcf> <popmap.tsv> --g=N [--nboot=N] [--boot=MODE] [--sites=N] [--min-n=N]\n")
  cat("       Rscript diversity_stats.R <vcf> <popmap.tsv> --g=N ... --complete-case\n\n")
  cat("  Works on populations.snps.vcf or populations.haps.vcf.\n")
  cat("  --g=N      REQUIRED. Rarefaction size in GENE COPIES (10 diploids = 20).\n")
  cat("  --nboot=N  bootstrap replicates (default 10000, 0 = none)\n")
  cat("  --boot=MODE  loci | individuals | both (default loci; the other two\n")
  cat("             are comparison-only -- they undercover, see\n")
  cat("             README.md, \"Bootstrap mode\").\n")
  cat("  --sites=N  sequenced nucleotide sites, for autosomal heterozygosity\n")
  cat("  --min-n=N  available-data mode only: min typed individuals per population\n")
  cat("             per locus (default 2). Ignored with --complete-case.\n")
  cat("  --complete-case   old rule: require every individual of every\n")
  cat("             population genotyped at a record before using it at all.\n")
}

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
         "`Rscript /path/to/diversity_stats.R ...`, or run ",
         "`R CMD INSTALL .` / `devtools::install()` from the package root ",
         "first).")
  }
  for (f in list.files(file.path(.d, "R"), pattern = "\\.R$", full.names = TRUE))
    source(f)
} else {
  library(RADdiversity)
}

## ---------------------------------------------------------------------------
## Parse argv into flags (--name=value or bare --complete-case) and
## positional args, rather than relying on fixed argument positions -- so
## `g` can be set without also having to supply `sites` as a placeholder.
## ---------------------------------------------------------------------------
args_raw <- commandArgs(trailingOnly = TRUE)
flag_names <- c("nboot", "sites", "g", "min-n", "boot")
flags <- list()
positional <- character(0)
for (a in args_raw) {
  if (identical(a, "--complete-case")) next
  m <- regmatches(a, regexec("^--([a-z-]+)=(.*)$", a))[[1]]
  if (length(m) == 3) {
    if (!(m[2] %in% flag_names)) {
      print_usage()
      stop("Unrecognized flag: --", m[2], "=... . Recognized: --",
           paste(flag_names, collapse = "=, --"), "=, and --complete-case.")
    }
    flags[[m[2]]] <- m[3]
  } else if (startsWith(a, "--") && sub("^--", "", a) %in% flag_names) {
    ## e.g. bare `--g 20` (space instead of `=`) -- a common typo of the
    ## `--g=20` form. Name the fix instead of a generic "unrecognized flag".
    print_usage()
    stop("--", sub("^--", "", a), " requires a value: --", sub("^--", "", a), "=N")
  } else if (startsWith(a, "--")) {
    print_usage()
    stop("Unrecognized flag: ", a)
  } else {
    positional <- c(positional, a)
  }
}
complete_case <- "--complete-case" %in% args_raw

if (length(positional) != 2) {
  print_usage()
  quit(status = 1, save = "no")
}
## Parse and validate every flag before touching the filesystem: fail fast on
## a missing/bad flag rather than reading (potentially large) files first
## only to error out anyway.
if (is.null(flags$g)) {
  print_usage()
  stop("Missing required --g=N (rarefaction size in gene copies). See usage above.")
}

## ---------------------------------------------------------------------------
## Call the package function; any error becomes a printed message and a
## non-zero exit status, since diversity_stats() itself must never call
## quit()/q() (that would kill an interactive RStudio session too).
## ---------------------------------------------------------------------------
res <- tryCatch(
  diversity_stats(
    vcf_file      = positional[1],
    popmap_f      = positional[2],
    g             = flags$g,
    nboot         = if (!is.null(flags$nboot)) flags$nboot else 10000L,
    boot          = if (!is.null(flags$boot)) flags$boot else "loci",
    sites         = if (!is.null(flags$sites)) flags$sites else 0,
    min_n         = if (!is.null(flags[["min-n"]])) flags[["min-n"]] else 2L,
    complete_case = complete_case
  ),
  error = function(e) {
    message("Error: ", conditionMessage(e))
    quit(status = 1, save = "no")
  }
)
