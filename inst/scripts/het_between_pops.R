#!/usr/bin/env Rscript
## Command-line entry point for RADdiversity::het_between_pops().
##   Rscript het_between_pops.R <vcf> <popmap.tsv> [--min-call=X] [--outdir=DIR]
##   Rscript het_between_pops.R --selftest
## Run it with no arguments for the full usage. All the logic, including the
## argument parsing, lives in the installed RADdiversity package (R/cli.R).
if (!requireNamespace("RADdiversity", quietly = TRUE))
  stop("RADdiversity is not installed. From R: install.packages(\"RADdiversity\")",
       call. = FALSE)
quit(status = RADdiversity:::.cli_main("het_between_pops", commandArgs(trailingOnly = TRUE)),
     save = "no")
