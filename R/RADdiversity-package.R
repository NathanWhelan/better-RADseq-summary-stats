#' RADdiversity: Population Diversity Statistics from Stacks RADseq Output
#'
#' Computes Ho, He, FIS (Nei & Chesser 1983, ratio of sums), rarefied allelic
#' richness and rarefied private allelic richness from Stacks 2 `populations`
#' VCF output, with block-bootstrap confidence intervals and delete-one-block
#' jackknife standard errors over RAD loci ([diversity_stats()]). Also
#' provides a correctly calibrated test of whether two populations differ in
#' mean heterozygosity, using the individual (not the locus) as the unit of
#' replication ([het_between_pops()]). See `README.md` in the package source
#' for the full statistical rationale.
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
