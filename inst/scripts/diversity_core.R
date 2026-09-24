#!/usr/bin/env Rscript
## Command-line entry point for RADdiversity::diversity_core_selftest().
##   Rscript diversity_core.R --selftest
## Checks every estimator against brute-force Monte Carlo and against
## published Stacks output. The logic lives in the installed RADdiversity
## package (R/cli.R, R/selftest.R).
if (!requireNamespace("RADdiversity", quietly = TRUE))
  stop("RADdiversity is not installed. From R: ",
       "remotes::install_github(\"NathanWhelan/better-RADseq-summary-stats\")",
       call. = FALSE)
quit(status = RADdiversity:::.cli_main("diversity_core", commandArgs(trailingOnly = TRUE)),
     save = "no")
