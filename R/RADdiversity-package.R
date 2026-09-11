#' RADdiversity: Population Genetic Statistics from RAD-seq Data
#'
#' Summary statistics for RAD-seq data, designed around Stacks 2
#' `populations` output:
#'
#' * [diversity_stats()]: Ho, He (Nei & Chesser 1983), FIS as a ratio of sums,
#'   rarefied allelic and private allelic richness, per-sequenced-site Ho/He.
#' * [het_between_pops()]: do populations differ in heterozygosity, with the
#'   individual (not the locus) as the unit of replication.
#' * [differentiation_stats()]: Weir & Cockerham FST, Jost's D and Weir &
#'   Goudet's beta.
#' * [kinship_check()] and [hwe_test()]: a relatedness screen and an exact
#'   Hardy-Weinberg test (a report, never a filter).
#' * [read_stacks_vcf()], the `filter_*()` functions and the `write_*()`
#'   functions: read, filter and export.
#'
#' Standard errors and confidence intervals resample whole RAD loci. The
#' command-line scripts are in `system.file("scripts", package =
#' "RADdiversity")`; the statistical rationale is in the package README.
#'
#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @importFrom stats cor cor.test median p.adjust quantile rnorm runif setNames sd t.test var wilcox.test
#' @importFrom utils head read.delim write.table
## usethis namespace: end
NULL

## Not exported. Captures the caller's current RNG state so it can be
## restored via on.exit() -- e.g. `restore <- .save_rng_state();
## on.exit(restore(), add = TRUE)` before calling set.seed(). Unlike the
## command-line scripts (each a fresh `Rscript` process with no "before"
## RNG state to preserve), diversity_stats(), het_between_pops() and both
## selftest functions are meant to be called interactively in a live R
## session -- without this, set.seed() would silently overwrite the
## caller's own random-number stream for the rest of that session.
.save_rng_state <- function() {
  if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
    old <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    function() assign(".Random.seed", old, envir = .GlobalEnv)
  } else {
    ## No prior RNG state (nothing has drawn a random number yet this
    ## session) -- restore to that same "unset" state rather than leaving
    ## a .Random.seed behind that didn't exist before.
    function() if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
      rm(".Random.seed", envir = .GlobalEnv)
  }
}
