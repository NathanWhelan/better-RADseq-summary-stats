###############################################################################
#
#  R/estimators.R  --  the five statistics this package exists to get right:
#                       Ho, He, FIS, rarefied allelic richness, and rarefied
#                       private allelic richness.
#
#  Called by diversity_stats(). Run diversity_core_selftest() to check every
#  estimator here against brute-force Monte Carlo, and the Stacks formulas
#  against real published Stacks output.
#
#  ---------------------------------------------------------------------------
#  1. EXPECTED HETEROZYGOSITY (He, gene diversity, Hs)
#  ---------------------------------------------------------------------------
#  Use Nei & Chesser (1983) / Nei (1987, eq. 7.39):
#
#       Hs = (n / (n - 1)) * (1 - sum_i p_i^2 - Ho / (2n))
#
#  with n the number of DIPLOID INDIVIDUALS genotyped, p_i the sample allele
#  frequencies, and Ho the observed heterozygote frequency.
#
#  WHY NOT (2n/(2n-1)) * (1 - sum p^2), which is what Stacks' `Pi` column is
#  and what most pipelines use? Because that correction assumes the 2n gene
#  copies are an independent sample of gametes, i.e. FIS = 0 -- the very
#  assumption under test. Writing H for the true gene diversity and F for the
#  true inbreeding coefficient:
#
#       E[1 - sum p_hat^2]              = H * (2n - 1 - F) / (2n)
#       E[(2n/(2n-1)) (1 - sum p^2)]    = H * (1 - F/(2n-1))     <- biased low
#       E[Nei-Chesser Hs]               = H                       <- unbiased
#
#  Consequences that matter here:
#    * FIS built on Stacks' Pi is biased TOWARD ZERO by a factor
#      [1 - (1-F)/(2n-1)]: about 3% at n = 15 and 5% at n = 10. The bias
#      depends on n, so it does not cancel when comparing two populations of
#      different size, and it makes FIS drift under rarefaction for a purely
#      arithmetic reason.
#    * FIS built on Stacks' Pi CANNOT REACH -1. With every individual
#      heterozygous the true FIS is -1, but 1 - Ho/Pi returns -(n-1)/n: -0.90
#      at n = 10. The Nei-Chesser version returns exactly -1. The self-test
#      checks this.
#
#  WHICH STACKS COLUMN IS WHICH. `populations.sumstats.tsv` has BOTH `Exp_Het`
#  and `Pi`, and they are different numbers estimating the same thing:
#
#      Exp_Het = 1 - sum p^2                    the plug-in estimate, biased low
#      Pi      = (2n/(2n-1)) (1 - sum p^2)      the same with the standard
#                                               unbiased-pi correction applied
#
#  So Pi = Exp_Het * 2n/(2n-1): a factor of 1.034 at n = 15, 1.053 at n = 10,
#  1.333 at n = 2. Verified against a published Stacks summary file
#  (Exp_Het 0.40000, Num_Indv 2, Pi 0.53333; 0.4 * 4/3 = 0.53333 exactly).
#  Compare anything from this package against `Pi`, never against `Exp_Het`.
#
#  Both estimators are computed here. hs_nei_chesser() is the default and what
#  diversity_stats() reports; hs_stacks_pi() / hs_stacks_pi_counts() reproduce
#  the Stacks quantity so the two can be printed side by side.
#
#  A CAVEAT ON THE NAME. At a BIALLELIC site, per-site pi and expected
#  heterozygosity are the same parameter -- the probability that two randomly
#  drawn gene copies differ. At a MULTI-ALLELIC locus they are not: gene
#  diversity (1 - sum p^2) counts every pair of distinct alleles as equally
#  different, while pi weights each pair by how many nucleotides actually
#  differ. Stacks keeps that distinction in populations.hapstats.tsv, which
#  reports Gene Diversity and Haplotype Diversity as separate columns. So on a
#  haplotype VCF nothing here is nucleotide diversity, whatever it is called.
#
#  ---------------------------------------------------------------------------
#  2. WHAT "He = 0.29" ACTUALLY MEANS (read this before quoting a number)
#  ---------------------------------------------------------------------------
#  Averaging He over the SNPs in a joint, multi-population call set gives
#  heterozygosity PER ASCERTAINED SNP. That number depends on how many SNPs
#  the run happened to call, so it is not comparable to any other study.
#
#  Schmidt et al. (2021) show that SNP-based heterozygosity is biased by sample
#  size and by analysing differentiated populations together, and that the
#  unbiased quantity is AUTOSOMAL heterozygosity: the same sum divided by every
#  sequenced site, monomorphic ones included. Retaining sites that are
#  monomorphic WITHIN a population but variable in another is a step toward
#  that, not the thing itself.
#
#  The conversion is one multiplication, so there is no excuse for not
#  reporting it:
#
#       He_autosomal = He_per_SNP * (n_SNPs_used / n_sites_sequenced)
#
#  `n_sites_sequenced` is the `Sites` column of the "All positions (variant and
#  fixed)" block of populations.sumstats_summary.tsv. Stacks reports both
#  blocks; the "All positions" one is the Schmidt-compliant estimate, and the
#  "Variant positions" one is not.
#
#  autosomal_het() below does the conversion and is used wherever a total-site
#  count is supplied.
#
#  ---------------------------------------------------------------------------
#  3. FIS
#  ---------------------------------------------------------------------------
#  FIS = 1 - sum(Ho) / sum(Hs), a ratio of sums over loci, never a mean of
#  per-locus ratios (Weir & Cockerham 1984; Bhatia et al. 2013). A site
#  monomorphic in a population adds 0 to both sums, so the ratio of sums does
#  not care whether such sites are included -- which is the whole reason to
#  prefer it.
#
#  ---------------------------------------------------------------------------
#  4. RAREFIED ALLELIC RICHNESS
#  ---------------------------------------------------------------------------
#       A_j(g) = sum_i [ 1 - C(N_j - N_ij, g) / C(N_j, g) ]
#
#  the expected number of distinct alleles in g GENE COPIES drawn without
#  replacement (Hurlbert 1971; El Mousadik & Petit 1996; Petit et al. 1998;
#  Kalinowski 2004). g is in gene copies: 10 diploids is g = 20.
#
#  5. RAREFIED PRIVATE ALLELIC RICHNESS
#       P_j(g) = sum_i [ Pr(allele i drawn in j) * prod_{k!=j} Pr(not drawn in k) ]
#  (Kalinowski 2004; Szpiech et al. 2008, implemented in HP-RARE and ADZE).
#
#  Both are analytic expectations over gene copies, which is what HP-RARE and
#  ADZE do. Do not substitute Monte-Carlo subsampling of INDIVIDUALS: that is a
#  different estimand whenever FIS != 0, because the two gene copies inside one
#  individual are not an independent draw.
#
#  Requires only base R.
#
###############################################################################


## ---------------------------------------------------------------------------
## He / Hs
## ---------------------------------------------------------------------------

#' Nei & Chesser (1983) unbiased gene diversity
#'
#' Unbiased expected heterozygosity (gene diversity, Hs) for one population at
#' one locus, vectorised over all three arguments. Unlike the estimator behind
#' Stacks' `Pi` column ([hs_stacks_pi()]), this stays unbiased at any FIS.
#'
#' @param sum_p2 Sum of squared sample allele frequencies.
#' @param ho Observed heterozygote frequency.
#' @param n Number of diploid individuals genotyped.
#' @return Numeric vector of Hs values; `NA` where `n < 2`.
#' @export
hs_nei_chesser <- function(sum_p2, ho, n) {
  out <- (n / (n - 1)) * (1 - sum_p2 - ho / (2 * n))
  out[!is.finite(out) | n < 2] <- NA_real_
  out
}

#' Biallelic wrapper around [hs_nei_chesser()]
#'
#' @param p Frequency of one allele.
#' @param ho Observed heterozygote frequency.
#' @param n Number of diploid individuals genotyped.
#' @return Numeric vector of Hs values.
#' @export
hs_biallelic <- function(p, ho, n) hs_nei_chesser(p^2 + (1 - p)^2, ho, n)

#' The estimator behind Stacks' `Pi` column (biallelic form)
#'
#' Unbiased when FIS = 0, biased low otherwise. Reported alongside
#' [hs_nei_chesser()] for comparison.
#'
#' @param p Frequency of one allele.
#' @param n Number of diploid individuals genotyped.
#' @return Numeric vector.
#' @export
hs_stacks_pi <- function(p, n) {
  out <- 2 * p * (1 - p) * (2 * n) / (2 * n - 1)
  out[!is.finite(out) | n < 1] <- NA_real_
  out
}

#' Gene diversity from gene-copy counts (multi-allelic)
#'
#' The same correction as [hs_stacks_pi()], computed from allele counts so it
#' also works on a haplotype VCF. Named for the correction it applies, not for
#' the software column it happens to match. `N` is the number of gene copies
#' actually observed, which equals `2n` on the complete-data locus set
#' [diversity_stats()] uses, so this agrees exactly with `hs_stacks_pi()`
#' there (checked in [diversity_core_selftest()]).
#'
#' @param counts Vector of gene-copy counts, one per allele.
#' @param n Unused; present for interface symmetry. Deprecated.
#' @return A single numeric value.
#' @export
gene_div_2n_counts <- function(counts, n = NULL) {
  N <- sum(counts)
  if (!is.finite(N) || N < 2) return(NA_real_)
  (N / (N - 1)) * (1 - sum((counts / N)^2))
}

#' Deprecated alias for [gene_div_2n_counts()]
#'
#' @param counts Vector of gene-copy counts, one per allele.
#' @param n Unused; present for interface symmetry. Deprecated.
#' @return A single numeric value.
#' @export
hs_stacks_pi_counts <- gene_div_2n_counts

#' Multi-allelic Nei-Chesser gene diversity from allele counts
#'
#' @param counts Vector of gene-copy counts for one population at one locus.
#' @param ho Observed heterozygote frequency.
#' @param n Number of diploid individuals genotyped.
#' @return A single numeric value.
#' @export
hs_from_counts <- function(counts, ho, n) {
  tot <- sum(counts)
  if (!is.finite(tot) || tot < 2 || n < 2) return(NA_real_)
  hs_nei_chesser(sum((counts / tot)^2), ho, n)
}

#' FIS as a ratio of sums over loci
#'
#' `FIS = 1 - sum(Ho) / sum(He)`, never a mean of per-locus ratios (Weir &
#' Cockerham 1984; Bhatia et al. 2013). A locus monomorphic in a population
#' adds 0 to both sums, so the ratio of sums does not care whether such loci
#' are included.
#'
#' @param ho Per-locus (or per-site) vector of observed heterozygosity.
#' @param he Per-locus (or per-site) vector of expected heterozygosity.
#' @return A single numeric value.
#' @export
fis_ratio_of_sums <- function(ho, he) {
  ok <- is.finite(ho) & is.finite(he)
  sh <- sum(he[ok]); if (!isTRUE(sh > 0)) return(NA_real_)
  1 - sum(ho[ok]) / sh
}

#' Per-ascertained-SNP heterozygosity to autosomal heterozygosity
#'
#' Schmidt et al. (2021) conversion. `n_sites_sequenced` is the `Sites` column
#' of the "All positions (variant and fixed)" block of
#' `populations.sumstats_summary.tsv`.
#'
#' Vectorised over `n_sites_sequenced` (as well as `het_per_snp`), so a
#' per-population sequenced-site count can be passed alongside a
#' per-population `het_per_snp` -- e.g. [diversity_stats()]'s `sites`
#' argument, which lets each population use its own denominator. `NA` is
#' returned element-wise wherever `n_sites_sequenced` is non-finite or
#' smaller than `n_snps_used`, rather than aborting the whole vector.
#'
#' @param het_per_snp Heterozygosity per ascertained SNP.
#' @param n_snps_used Number of SNPs used to compute `het_per_snp`.
#' @param n_sites_sequenced Total sequenced sites (variant and fixed). A
#'   single value or one per element of `het_per_snp`.
#' @return A numeric vector, the same length as `het_per_snp`.
#' @export
autosomal_het <- function(het_per_snp, n_snps_used, n_sites_sequenced) {
  bad <- !is.finite(n_sites_sequenced) | n_sites_sequenced < n_snps_used
  out <- het_per_snp * n_snps_used / n_sites_sequenced
  out[bad] <- NA_real_
  out
}


## ---------------------------------------------------------------------------
## Rarefaction
## ---------------------------------------------------------------------------

#' Probability each allele is sampled in g gene copies
#'
#' Pr(each allele appears at least once in g gene copies drawn without
#' replacement). Uses `lchoose` for numerical stability at large N.
#'
#' @param counts Vector of gene-copy counts, one per allele.
#' @param g Rarefaction size in gene copies.
#' @return Numeric vector, one probability per allele.
#' @export
p_sampled <- function(counts, g) {
  N <- sum(counts)
  if (!is.finite(N) || g > N || g < 1) return(rep(NA_real_, length(counts)))
  1 - exp(lchoose(N - counts, g) - lchoose(N, g))
}

#' Rarefied allelic richness at one locus
#'
#' @param counts Vector of gene-copy counts, one per allele.
#' @param g Rarefaction size in gene copies.
#' @return A single numeric value.
#' @export
rare_richness <- function(counts, g) {
  if (!is.finite(sum(counts)) || sum(counts) < g) return(NA_real_)
  ps <- p_sampled(counts, g)
  if (anyNA(ps)) return(NA_real_) else sum(ps)
}

#' Rarefied private allelic richness for one population
#'
#' @param count_mat Matrix of gene-copy counts, populations x alleles.
#' @param j Row index (population) to compute private richness for.
#' @param g Rarefaction size in gene copies.
#' @return A single numeric value.
#' @export
rare_private <- function(count_mat, j, g) {
  if (!is.matrix(count_mat)) count_mat <- rbind(count_mat)
  ps <- matrix(NA_real_, nrow(count_mat), ncol(count_mat))
  for (rr in seq_len(nrow(count_mat))) ps[rr, ] <- p_sampled(count_mat[rr, ], g)
  if (anyNA(ps)) return(NA_real_)
  term <- ps[j, ]
  for (k in setdiff(seq_len(nrow(count_mat)), j)) term <- term * (1 - ps[k, ])
  sum(term)
}

#' Rarefied private allelic richness for every population
#'
#' @param count_mat Matrix of gene-copy counts, populations x alleles.
#' @param g Rarefaction size in gene copies.
#' @return Numeric vector, length `nrow(count_mat)`.
#' @export
rare_private_all <- function(count_mat, g) {
  vapply(seq_len(nrow(count_mat)), function(j) rare_private(count_mat, j, g),
         numeric(1))
}
