###############################################################################
#
#  R/individual_inbreeding.R -- a method-of-moments inbreeding coefficient for
#  each individual.
#
#    F_i = 1 - O_i / E_i
#
#  O_i  the number of loci at which individual i is heterozygous
#  E_i  the heterozygosity expected at those same loci under random mating in
#       i's own population: the sum over the loci i is typed at of Nei &
#       Chesser's (1983) Hs, from that population's allele frequencies
#
#  This is the form PLINK's --het reports (Purcell et al. 2007: observed vs
#  expected homozygosity gives the same F), with one deliberate difference in
#  the expectation. PLINK uses the plug-in 1 - sum(p^2), or with its
#  small-sample option the 2n/(2n-1) correction -- the estimator behind
#  Stacks' `Pi`, which is biased whenever F != 0. Hs is unbiased at any F, and
#  with it the population mean of F_i equals that population's FIS from
#  diversity_stats() (1 - sum(Ho)/sum(Hs)) whenever every individual is typed
#  at every locus: E_i is then the same for everyone, so
#  mean_i(1 - O_i/E) = 1 - sum_l Ho_l / sum_l Hs_l.
#
#  F_i is relative to the individual's OWN population, so comparing F between
#  populations (het_between_pops()) asks whether they differ in INBREEDING;
#  comparing individual heterozygosity asks whether they differ in DIVERSITY.
#
###############################################################################

## Not exported. F_i for every individual of ONE population. `a1`/`a2` are
## loci x individuals allele matrices already restricted to the loci to use.
## Loci where Hs is undefined (fewer than 2 typed individuals) are skipped.
## Returns a data frame, one row per individual (columns of `a1`).
.ind_F <- function(a1, a2) {
  nr <- nrow(a1); nc <- ncol(a1)
  typed <- !is.na(a1)
  het   <- typed & (a1 != a2)                  # FALSE where untyped
  n_l   <- rowSums(typed)
  ho_l  <- rowSums(het) / n_l
  cnt  <- .allele_counts(a1, a2)
  hs_l <- hs_nei_chesser(rowSums((cnt / (2 * n_l))^2), ho_l, n_l)   # NA if n < 2
  use  <- is.finite(hs_l)
  obs  <- colSums(het & use)
  expd <- colSums(ifelse(typed & use, hs_l, 0))
  data.frame(n_loci = colSums(typed & use), obs_het = obs, exp_het = expd,
             F = ifelse(expd > 0, 1 - obs / expd, NA_real_), row.names = NULL)
}

#' Individual inbreeding coefficients
#'
#' A method-of-moments inbreeding coefficient for every individual:
#' `F = 1 - O/E`, where `O` is the number of loci at which the individual is
#' heterozygous and `E` the heterozygosity expected at those same loci under
#' random mating in its own population -- the sum of Nei & Chesser's (1983)
#' unbiased gene diversity over the loci the individual is typed at. This is
#' the statistic PLINK's `--het` reports, with an expectation that stays
#' unbiased when F is not 0; as a result, when every individual is typed at
#' every locus the population mean of `F` equals that population's FIS from
#' [diversity_stats()].
#'
#' `F` is relative to each individual's own population. Comparing it between
#' populations ([het_between_pops()] does this) asks whether they differ in
#' inbreeding; comparing individual heterozygosity asks whether they differ in
#' diversity. Negative values mean more heterozygosity than random mating
#' predicts.
#'
#' @param H A list as returned by [read_stacks_vcf()] (optionally filtered).
#' @param pops A named list of sample-ID vectors, one per population, as
#'   returned by [read_popmap()].
#' @param min_call Use a locus for a population only if at least this fraction
#'   of that population's individuals is genotyped there. Default `0.9`.
#' @return A data frame, one row per individual: `sample`, `population`,
#'   `n_loci` (loci used), `obs_het` (heterozygous loci), `exp_het` (expected
#'   number) and `F`.
#' @references
#' Nei, M. & Chesser, R.K. (1983) Estimation of fixation indices and gene
#' diversities. *Annals of Human Genetics* 47:253-259.
#'
#' Purcell, S. et al. (2007) PLINK: a tool set for whole-genome association
#' and population-based linkage analyses. *American Journal of Human
#' Genetics* 81:559-575.
#' @examples
#' H <- read_stacks_vcf(system.file("extdata", "small.haps.vcf",
#'                                  package = "RADdiversity"), verbose = FALSE)
#' pops <- read_popmap(system.file("extdata", "small_popmap.tsv",
#'                                 package = "RADdiversity"), H$samples)
#' individual_inbreeding(H, pops, min_call = 0.5)
#' @export
individual_inbreeding <- function(H, pops, min_call = 0.9) {
  H <- .resolve_H(H, verbose = FALSE)
  if (!is.list(pops) || is.null(names(pops)))
    stop("pops must be a named list of sample IDs per population (see read_popmap()).")
  if (!(is.numeric(min_call) && length(min_call) == 1L && min_call >= 0 && min_call <= 1))
    stop("min_call must be a single number between 0 and 1.")
  keep <- sweep(.typed_by_pop(H, pops), 2L, lengths(pops), "/") >= min_call - 1e-9
  out <- do.call(rbind, lapply(names(pops), function(p) {
    rows <- keep[, p]
    data.frame(sample = pops[[p]], population = p,
               .ind_F(H$A1[rows, pops[[p]], drop = FALSE],
                      H$A2[rows, pops[[p]], drop = FALSE]),
               row.names = NULL)
  }))
  out
}
