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
  typed <- !is.na(a1)                           # loci x individuals
  het <- typed & (a1 != a2)                     # FALSE where untyped
  ## Per locus: genotyped individuals, observed heterozygosity, Hs.
  n_typed <- rowSums(typed)
  ho_locus <- rowSums(het) / n_typed
  counts <- .allele_counts(a1, a2)
  hs_locus <- hs_nei_chesser(rowSums((counts / (2 * n_typed))^2), ho_locus, n_typed)
  usable <- is.finite(hs_locus)                 # NA when fewer than 2 typed
  ## Per individual: observed heterozygous loci, and the sum of Hs over the
  ## usable loci it is typed at. (hs_locus has one value per row, so
  ## ifelse() recycles it down each individual's column.)
  observed <- colSums(het & usable)
  expected <- colSums(ifelse(typed & usable, hs_locus, 0))
  data.frame(n_loci = as.integer(colSums(typed & usable)), obs_het = as.integer(observed),
             exp_het = expected,
             F = ifelse(expected > 0, 1 - observed / expected, NA_real_), row.names = NULL)
}

#' Individual inbreeding coefficients
#'
#' A method-of-moments inbreeding coefficient for every individual:
#' `F = 1 - O/E`, where `O` is the number of loci at which the individual is
#' heterozygous and `E` the heterozygosity expected at those same loci under
#' random mating in its own population -- the sum of Nei & Chesser's (1983)
#' unbiased gene diversity over the loci the individual is typed at. This is
#' the statistic PLINK's `--het` reports, with two differences in the
#' expectation. PLINK takes allele frequencies from the whole sample (or a
#' file you supply), whereas here they come from each individual's own
#' population. And PLINK uses the plug-in 1 - sum(p^2), whereas here the
#' estimator stays unbiased when F is not 0. As a result, when every
#' individual is typed at every locus, the population mean of `F` equals that
#' population's FIS from [diversity_stats()].
#'
#' `F` is relative to each individual's own population. Comparing it between
#' populations ([het_between_pops()] does this) asks whether they differ in
#' inbreeding; comparing individual heterozygosity asks whether they differ in
#' diversity. Negative values mean more heterozygosity than random mating
#' predicts.
#'
#' @param vcf Path to a VCF file, or the object returned by [read_stacks_vcf()]
#'   (optionally filtered).
#' @param popmap Path to a popmap file, or the list returned by
#'   [read_popmap()].
#' @param min_call Use a locus for a population only if at least this fraction
#'   of that population's individuals is genotyped there. Default `0.9`.
#' @param verbose Print progress messages (when `vcf` or `popmap` is a path).
#'   Default `TRUE`.
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
#' vcf    <- system.file("extdata", "small.haps.vcf", package = "RADdiversity")
#' popmap <- system.file("extdata", "small_popmap.tsv", package = "RADdiversity")
#' individual_inbreeding(vcf, popmap, min_call = 0.5, verbose = FALSE)
#' @export
individual_inbreeding <- function(vcf, popmap, min_call = 0.9, verbose = TRUE) {
  .check_number(min_call, "min_call", min = 0, max = 1)
  .check_flag(verbose, "verbose")
  H <- .resolve_H(vcf, verbose = verbose)
  pops <- .resolve_pops(popmap, H$samples, verbose = verbose)
  ## A locus is used for population p when p itself genotypes >= min_call of
  ## its individuals there.
  locus_sets <- .population_locus_sets(H, pops, min_call)
  do.call(rbind, lapply(names(pops), function(p) {
    loci <- locus_sets[, p]
    data.frame(sample = pops[[p]], population = p,
               .ind_F(H$A1[loci, pops[[p]], drop = FALSE],
                      H$A2[loci, pops[[p]], drop = FALSE]),
               row.names = NULL)
  }))
}
