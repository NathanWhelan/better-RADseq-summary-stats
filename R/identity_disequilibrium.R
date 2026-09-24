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
#  correlation, and it is 0 when individuals do not differ.
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
#  LINKED SNPs. g2 assumes the two loci of a pair are inherited
#  independently. SNPs on one RAD tag are not: an individual carrying two
#  different haplotypes is heterozygous at several of them at once, which
#  looks exactly like identity disequilibrium. So only pairs of records on
#  DIFFERENT RAD loci are counted. With R_Bi = het records of individual i on
#  RAD locus B, and C_B = sum over B's records of c_l:
#
#       P_s = sum_i (r_i^2 - sum_B R_Bi^2)      P_d = (sum c)^2 - sum_B C_B^2 - P_s
#
#  and Q_s, Q_d the same with typed counts. With one record per RAD locus
#  (a haplotype VCF) this is the formula above. Without it, in a simulation
#  with no inbreeding and about 1.8 SNPs per RAD tag, the permutation test
#  rejected g2 = 0 in 54-56% of runs at a nominal 5%, from 500 to 8,000 RAD
#  loci; with it, in 2-4%, and real variance in inbreeding was still found
#  every time (inst/sims/g2_linked_snps.R).
#
#  Uncertainty: a bootstrap over individuals for the CI, and a permutation
#  test of g2 = 0 that shuffles heterozygosity among individuals, one RAD
#  locus at a time (keeping each locus's heterozygote count and every
#  individual's missing data) -- the usual choices for g2.
#
###############################################################################

## Not exported. g2 for each column of W, a matrix of individual weights (one
## column per estimate: all 1 for the point estimate, bootstrap draw counts
## for a replicate -- a duplicated individual counts as distinct copies).
## `h` and `m` are RAD loci x individuals matrices: how many of the
## individual's records on that locus are heterozygous, and typed (0/1 when
## each RAD locus has one record). Only pairs of records on different RAD
## loci count (see "LINKED SNPs" above).
.g2_weighted <- function(h, m, W) {
  het_records <- colSums(h)       # r_i: heterozygous records of each individual
  typed_records <- colSums(m)     # t_i: typed records of each individual
  ## Pairs within one individual, minus those on the same RAD locus.
  P_same <- colSums((het_records^2 - colSums(h^2)) * W)
  Q_same <- colSums((typed_records^2 - colSums(m^2)) * W)
  het_individuals <- h %*% W      # C_B: heterozygous records on each RAD locus
  typed_individuals <- m %*% W    # the same for typed records
  P_diff <- colSums(het_individuals)^2 - colSums(het_individuals^2) - P_same
  Q_diff <- colSums(typed_individuals)^2 - colSums(typed_individuals^2) - Q_same
  g2 <- (P_same / Q_same) / (P_diff / Q_diff) - 1
  g2[!is.finite(g2)] <- NA_real_
  g2
}

## Not exported. Permutation null for g2: `nperm` data sets in which, at each
## RAD locus, the individuals typed at every record of that locus swap their
## whole heterozygosity at random (all records of the locus move together, so
## the linkage within a locus is kept). Individuals typed at only some of a
## locus's records keep theirs. Typed counts, and each locus's heterozygote
## count, are unchanged, so only the same-individual sum P_s changes. `h` and
## `m` are as in .g2_weighted(); `records` is how many records each RAD locus
## has. With one record per locus this is the classic shuffle of each locus's
## heterozygotes among its typed individuals.
.g2_perm <- function(h, m, nperm, records = rep(1, nrow(h))) {
  n_loci <- nrow(h)
  n_ind <- ncol(h)
  ## The cells that are shuffled: individuals typed at every record of a locus.
  cell <- which(m == records)            # `records` recycles down each column
  locus_of <- (cell - 1L) %% n_loci + 1L
  individual_of <- (cell - 1L) %/% n_loci + 1L
  value <- h[cell]
  ## Each locus's values, largest first, in the order the shuffled cells are
  ## filled (so, with 0/1 values, the first c_l shuffled cells become the
  ## heterozygotes).
  values_in_order <- value[order(locus_of, -value)]
  ## What the cells that are never shuffled add to each individual's sums.
  fixed <- h
  fixed[cell] <- 0
  fixed_sum <- colSums(fixed)
  fixed_sum_sq <- colSums(fixed^2)
  per_individual <- function(x, who) {
    out <- numeric(n_ind)
    s <- rowsum(x, who)
    out[as.integer(rownames(s))] <- s
    out
  }
  Q_same <- sum(colSums(m)^2 - colSums(m^2))
  Q_diff <- sum(rowSums(m))^2 - sum(rowSums(m)^2) - Q_same
  P_total <- sum(rowSums(h))^2 - sum(rowSums(h)^2)
  vapply(seq_len(nperm), function(k) {
    ## Sorting on locus + a uniform keeps each locus's cells together and puts
    ## them in random order; the values then fill them in that order.
    o <- order(locus_of + stats::runif(length(locus_of)))
    who <- individual_of[o]
    het_sum <- fixed_sum + per_individual(values_in_order, who)
    het_sum_sq <- fixed_sum_sq + per_individual(values_in_order^2, who)
    P_same <- sum(het_sum^2 - het_sum_sq)
    (P_same / Q_same) / ((P_total - P_same) / Q_diff) - 1
  }, numeric(1))
}

## Not exported. g2 with a bootstrap CI (and optionally a permutation
## p-value) for ONE population. `a1`/`a2`: records x individuals allele
## matrices, restricted to the records to use. `block`: the RAD locus of each
## record (H$locus_raw), so that pairs of records on one RAD locus are left
## out; NULL treats every record as its own locus.
.g2_summary <- function(a1, a2, nboot, nperm, block = NULL) {
  m <- !is.na(a1)
  h <- m & (a1 != a2)
  storage.mode(h) <- "double"
  storage.mode(m) <- "double"
  records <- rep(1, nrow(h))
  if (!is.null(block) && anyDuplicated(block)) {
    block <- match(block, unique(block))
    h <- rowsum(h, block, reorder = FALSE)
    m <- rowsum(m, block, reorder = FALSE)
    records <- as.numeric(tabulate(block))
  }
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
    null <- .g2_perm(h, m, nperm, records)
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
#' resample loci ([diversity_stats()]'s `_se`, `_lo`, `_hi`) hold the sampled
#' individuals fixed. When individuals are alike that costs little, because
#' with many loci the chance of which individuals were caught averages out
#' (Nei & Roychoudhury 1974; Nei 1978). When individuals differ in inbreeding,
#' heterozygosity is correlated across loci within an individual and it does
#' not average out: locus-based SEs for Ho and FIS are then too small. A g2
#' interval above 0 is that signal. Report `_se_combined` from
#' [diversity_stats()] (computed by default; see "How standard errors are
#' calculated" there) and compare populations with
#' [het_between_pops()]. Relatives in a sample also raise g2; screen for them
#' with [kinship_check()].
#'
#' **Estimator.** From the definition, using every pair of loci typed in the
#' same individual (numerator) and in two different individuals
#' (denominator); with no missing data this is the closed form of Hoffman et
#' al. (2014). The confidence interval is a bootstrap over individuals; the
#' p-value comes from permutations that shuffle heterozygosity among
#' individuals, one RAD locus at a time (`(exceedances + 1) / (nperm + 1)`).
#' Pooling structured populations inflates g2 (Wahlund effect), hence one
#' estimate per population.
#'
#' **Test g2 > 0 with `p_value`.** With few individuals the bootstrap interval
#' is too narrow. In simulations with no inbreeding (true g2 = 0), the 95%
#' interval missed 0 in 17-19% of data sets with 8 individuals, 10-15% with
#' 15 and 7% with 30, almost always by lying below 0, while the permutation
#' test rejected in 2-5% (`inst/sims/g2_ci_small_n.R`). inbreedR uses the same
#' bootstrap. So with fewer than about 30 individuals, read the interval as a
#' rough guide to size and use `p_value` to decide whether g2 is above 0.
#'
#' **SNP or haplotype VCF.** Both work. SNPs on one RAD locus are inherited
#' together, so an individual is often heterozygous at several of them at
#' once even when nobody is inbred. Pairs of records on the same RAD locus
#' (`locus_raw`, see [read_stacks_vcf()]) are therefore left out, and the
#' permutations move a RAD locus's records together. Counting those pairs
#' made the test reject a true g2 of 0 in more than half of simulated data
#' sets; without them, 2-4% (`inst/sims/g2_linked_snps.R`).
#'
#' **Failed individuals.** An individual that looks like a failed library
#' (genotyped at fewer than 50 records when the rest of its population has
#' that many, or at under half of the records where the rest of its population
#' is genotyped) is kept and named in a warning: allele dropout lowers its
#' heterozygosity at many loci at once, which looks like variance in
#' inbreeding and inflates g2. See "Failed individuals" in
#' [het_between_pops()]. A population with fewer than 2 individuals, or with
#' no locus at `min_call`, stops the function with a message naming it.
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
#'   `n_loci` (RAD loci used), `g2`, `g2_se` (bootstrap SD), `g2_lo`,
#'   `g2_hi` (95% bootstrap interval) and `p_value` (one-sided, g2 > 0).
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
  .check_run_inputs(vcf, popmap, stem = NULL, outdir = NULL)
  H <- .resolve_H(vcf, verbose = verbose)
  pops <- .resolve_pops(popmap, H$samples, verbose = verbose)
  .stop_tiny_pops(pops, "g2 compares heterozygosity between individuals, so it needs at least 2.")
  .warn_failed_individuals(H, pops)
  if (!is.null(seed)) {
    restore_rng <- .save_rng_state()
    on.exit(restore_rng(), add = TRUE)
    set.seed(seed)
  }
  locus_sets <- .population_locus_sets(H, pops, min_call)
  .check_population_loci(H, pops, locus_sets, min_call)
  do.call(rbind, lapply(names(pops), function(p) {
    loci <- locus_sets[, p]
    cbind(population = p,
          .g2_summary(H$A1[loci, pops[[p]], drop = FALSE],
                      H$A2[loci, pops[[p]], drop = FALSE], nboot, nperm,
                      block = H$locus_raw[loci]))
  }))
}
