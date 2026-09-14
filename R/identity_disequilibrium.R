###############################################################################
#
#  R/identity_disequilibrium.R -- g2, the standard measure of variance in
#  inbreeding among individuals (identity disequilibrium).
#
#  WHY IT MATTERS HERE. The package's locus-based standard errors and
#  intervals (diversity_stats()) treat the sampled individuals as fixed:
#  they answer "what if I had typed different loci in these same animals?".
#  That is also enough for population-level inference as long as individuals
#  are exchangeable -- with many loci, the variation that comes from which
#  individuals were sampled averages away (Nei & Roychoudhury 1974; Nei 1978).
#  It stops averaging away when individuals differ in inbreeding (or include
#  relatives): heterozygosity is then correlated across loci WITHIN an
#  individual, and no number of loci removes that. g2 measures exactly this
#  correlation, and it is 0 when individuals do not differ. The overdispersion
#  factor het_between_pops() prints is the same idea in a cruder form.
#
#  THE ESTIMATOR. With h_il = 1 if individual i is heterozygous at locus l,
#  g2 is defined (David et al. 2007) as
#
#       g2 = E[h_l h_l'] / (E[h_l] E[h_l']) - 1          for loci l != l'
#
#  i.e. how much more often an individual is heterozygous at two loci at once
#  than two different individuals are. Estimated here directly from that
#  definition: the numerator averages h_il h_il' over every pair of distinct
#  loci typed in the same individual, the denominator averages h_il h_jl'
#  over every pair of distinct loci typed in two different individuals. With
#  r_i = het loci of individual i, t_i = its typed loci, c_l = het individuals
#  at locus l, u_l = typed individuals at locus l:
#
#       same-individual:  P_s = sum_i (r_i^2 - r_i)       Q_s = sum_i (t_i^2 - t_i)
#       different:        P_d = (sum c)^2 - sum c^2 - P_s  Q_d = (sum u)^2 - sum u^2 - Q_s
#       g2 = (P_s / Q_s) / (P_d / Q_d) - 1
#
#  With no missing data Q_d / Q_s = n - 1 and this is the closed form
#  (n - 1) P_s / P_d - 1 of Hoffman et al. (2014); with missing data it uses
#  exactly the locus pairs that were typed, rather than an approximation.
#  Written from the definition above -- no code from other implementations.
#
#  Uncertainty: a bootstrap over individuals for the CI, and a permutation
#  test of g2 = 0 that shuffles heterozygosity among each locus's typed
#  individuals (keeping every locus's heterozygote count and every
#  individual's missing data) -- the usual choices for g2.
#
###############################################################################

## Not exported. g2 for each column of W, a matrix of individual weights (one
## column per estimate: all 1 for the point estimate, bootstrap draw counts
## for a replicate -- a duplicated individual counts as distinct copies). `h`
## and `m` are loci x individuals 0/1 matrices: heterozygous, typed.
.g2_weighted <- function(h, m, W) {
  het_loci <- colSums(h)        # r_i: heterozygous loci of each individual
  typed_loci <- colSums(m)      # t_i: typed loci of each individual
  P_same <- colSums((het_loci^2 - het_loci) * W)
  Q_same <- colSums((typed_loci^2 - typed_loci) * W)
  het_individuals <- h %*% W    # c_l: heterozygous individuals at each locus
  typed_individuals <- m %*% W  # u_l: typed individuals at each locus
  P_diff <- colSums(het_individuals)^2 - colSums(het_individuals^2) - P_same
  Q_diff <- colSums(typed_individuals)^2 - colSums(typed_individuals^2) - Q_same
  g2 <- (P_same / Q_same) / (P_diff / Q_diff) - 1
  g2[!is.finite(g2)] <- NA_real_
  g2
}

## Not exported. Permutation null for g2: `nperm` data sets in which each
## locus's heterozygotes are re-assigned at random among that locus's typed
## individuals. Heterozygote counts per locus (c), and all typed counts (t,
## u), are unchanged, so only the same-individual sum P_s changes.
.g2_perm <- function(h, m, nperm) {
  n_loci <- nrow(h)
  n_ind <- ncol(h)
  ## Every typed cell, as its (locus, individual) position.
  typed_cell <- which(m == 1)
  locus_of <- (typed_cell - 1L) %% n_loci + 1L
  individual_of <- (typed_cell - 1L) %/% n_loci + 1L
  het_per_locus <- rowSums(h)
  typed_per_locus <- rowSums(m)
  typed_loci <- colSums(m)
  cells_before_locus <- cumsum(c(0, typed_per_locus))[seq_len(n_loci)]
  Q_same <- sum(typed_loci^2 - typed_loci)
  Q_diff <- sum(typed_per_locus)^2 - sum(typed_per_locus^2) - Q_same
  P_total <- sum(het_per_locus)^2 - sum(het_per_locus^2)
  vapply(seq_len(nperm), function(k) {
    ## Shuffle typed cells within each locus (sorting on locus + a uniform
    ## keeps each locus's cells together), then make the first c_l cells of
    ## locus l its heterozygotes.
    o <- order(locus_of + stats::runif(length(locus_of)))
    rank_in_locus <- seq_along(o) - cells_before_locus[locus_of[o]]
    het_loci <- tabulate(individual_of[o][rank_in_locus <= het_per_locus[locus_of[o]]],
                         nbins = n_ind)
    P_same <- sum(het_loci^2 - het_loci)
    (P_same / Q_same) / ((P_total - P_same) / Q_diff) - 1
  }, numeric(1))
}

## Not exported. g2 with a bootstrap CI (and optionally a permutation
## p-value) for ONE population. `a1`/`a2`: loci x individuals allele
## matrices, restricted to the loci to use.
.g2_summary <- function(a1, a2, nboot, nperm) {
  m <- !is.na(a1)
  h <- m & (a1 != a2)
  storage.mode(h) <- "double"
  storage.mode(m) <- "double"
  n_ind <- ncol(h)
  estimate <- .g2_weighted(h, m, matrix(1, n_ind, 1))
  lo <- hi <- se <- p_value <- NA_real_
  if (nboot > 0) {
    ## Each column of `times_drawn` is one bootstrap sample of individuals.
    boot <- .boot_in_batches(n_ind, nboot,
                             function(times_drawn) cbind(.g2_weighted(h, m, times_drawn)),
                             cells_per_replicate = nrow(h))[, 1]
    lo <- stats::quantile(boot, 0.025, na.rm = TRUE, names = FALSE)
    hi <- stats::quantile(boot, 0.975, na.rm = TRUE, names = FALSE)
    se <- stats::sd(boot, na.rm = TRUE)
  }
  if (nperm > 0 && is.finite(estimate)) {
    null <- .g2_perm(h, m, nperm)
    p_value <- (sum(null >= estimate, na.rm = TRUE) + 1) / (nperm + 1)
  }
  data.frame(n_ind = n_ind, n_loci = nrow(h), g2 = estimate, g2_se = se,
             g2_lo = lo, g2_hi = hi, p_value = p_value)
}

#' Identity disequilibrium (g2): do individuals differ in inbreeding?
#'
#' Estimates g2 (David et al. 2007) in each population: how much more often an
#' individual is heterozygous at two loci at once than two different
#' individuals are. g2 is 0 when individuals do not differ in inbreeding and
#' positive when they do.
#'
#' @details
#' **Why it matters for this package.** Standard errors and intervals that
#' resample loci ([diversity_stats()]) hold the sampled individuals fixed.
#' With many loci that is also enough for population-level inference, as
#' long as individuals are exchangeable (Nei & Roychoudhury 1974; Nei 1978).
#' When individuals differ in inbreeding -- or include relatives --
#' heterozygosity is correlated across loci within an individual and the
#' locus-based intervals are too narrow. A g2 confidence interval that
#' excludes 0 is the signal to use the individual-jackknife standard errors
#' of [diversity_stats()] (`se_individuals = TRUE`) and the individual-level
#' tests of [het_between_pops()].
#'
#' **Estimator.** From the definition, using every pair of distinct loci
#' typed in the same individual (numerator) and in two different individuals
#' (denominator); with no missing data this is the closed form of Hoffman et
#' al. (2014). The confidence interval is a bootstrap over individuals; the
#' p-value comes from permutations that shuffle heterozygosity among each
#' locus's typed individuals (`(exceedances + 1) / (nperm + 1)`). Any marker
#' type works -- SNPs or multi-allelic RAD haplotypes. Pooling structured
#' populations inflates g2 (Wahlund effect), hence one estimate per
#' population.
#'
#' **Failed individuals.** An individual that looks like a failed library
#' (genotyped at fewer than 50 records when the rest of its population has
#' that many, or at under half of the records where the rest of its population
#' is genotyped) is kept and named in a warning:
#' allele dropout lowers its heterozygosity at many loci at once, which looks
#' like variance in inbreeding and inflates g2. See "Failed individuals" in
#' [het_between_pops()].
#'
#' @param vcf Path to a VCF file, or the object returned by [read_stacks_vcf()]
#'   (optionally filtered).
#' @param popmap Path to a popmap file, or the list returned by
#'   [read_popmap()].
#' @param nboot Bootstrap replicates over individuals for the CI. Default
#'   `1000`; `0` for none.
#' @param nperm Permutations for the p-value. Default `1000`; `0` for none.
#' @param min_call Use a locus for a population only if at least this fraction
#'   of that population's individuals is genotyped there. Default `0.9`.
#' @param seed Random seed. Default `NULL`: use R's current random-number
#'   stream (call `set.seed()` first for reproducible results). A number makes
#'   the result reproducible on its own and leaves your session's
#'   random-number stream as it was.
#' @param verbose Print progress messages. Default `TRUE`.
#' @return A data frame, one row per population: `population`, `n_ind`,
#'   `n_loci`, `g2`, `g2_se` (bootstrap SD), `g2_lo`, `g2_hi` (95% bootstrap
#'   interval) and `p_value` (one-sided, g2 > 0).
#' @references
#' David, P., Pujol, B., Viard, F., Castella, V. & Goudet, J. (2007) Reliable
#' selfing rate estimates from imperfect population genetic data. *Molecular
#' Ecology* 16:2474-2487.
#'
#' Hoffman, J.I., Simpson, F., David, P., Rijks, J.M., Kuiken, T., Thorne,
#' M.A.S., Lacy, R.C. & Dasmahapatra, K.K. (2014) High-throughput sequencing
#' reveals inbreeding depression in a natural population. *PNAS*
#' 111:3775-3780.
#'
#' Nei, M. (1978) Estimation of average heterozygosity and genetic distance
#' from a small number of individuals. *Genetics* 89:583-590.
#'
#' Nei, M. & Roychoudhury, A.K. (1974) Sampling variances of heterozygosity
#' and genetic distance. *Genetics* 76:379-390.
#' @examples
#' vcf    <- system.file("extdata", "small.haps.vcf", package = "RADdiversity")
#' popmap <- system.file("extdata", "small_popmap.tsv", package = "RADdiversity")
#' set.seed(1)
#' identity_disequilibrium(vcf, popmap, nboot = 200, nperm = 200, min_call = 0.5,
#'                         verbose = FALSE)
#' @export
identity_disequilibrium <- function(vcf, popmap, nboot = 1000L, nperm = 1000L,
                                    min_call = 0.9, seed = NULL, verbose = TRUE) {
  nboot <- .check_count(nboot, "nboot")
  nperm <- .check_count(nperm, "nperm")
  .check_number(min_call, "min_call", min = 0, max = 1)
  .check_seed(seed)
  .check_flag(verbose, "verbose")
  H <- .resolve_H(vcf, verbose = verbose)
  pops <- .resolve_pops(popmap, H$samples, verbose = verbose)
  .warn_failed_individuals(H, pops)
  if (!is.null(seed)) {
    restore_rng <- .save_rng_state()
    on.exit(restore_rng(), add = TRUE)
    set.seed(seed)
  }
  locus_sets <- .population_locus_sets(H, pops, min_call)
  do.call(rbind, lapply(names(pops), function(p) {
    loci <- locus_sets[, p]
    cbind(population = p,
          .g2_summary(H$A1[loci, pops[[p]], drop = FALSE],
                      H$A2[loci, pops[[p]], drop = FALSE], nboot, nperm))
  }))
}
