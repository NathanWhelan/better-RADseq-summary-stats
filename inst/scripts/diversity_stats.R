#!/usr/bin/env Rscript
## Command-line entry point for RADdiversity::diversity_stats().
##   Rscript diversity_stats.R <vcf> <popmap.tsv> --g=N [options]
## Run it with no arguments for the full usage. All the logic, including the
## argument parsing, lives in the installed RADdiversity package (R/cli.R).
if (!requireNamespace("RADdiversity", quietly = TRUE))
  stop("RADdiversity is not installed. From R: ",
       "remotes::install_github(\"NathanWhelan/better-RADseq-summary-stats\")",
       call. = FALSE)
quit(status = RADdiversity:::.cli_main("diversity_stats", commandArgs(trailingOnly = TRUE)),
     save = "no")
