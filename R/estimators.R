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
#  diversity_stats() reports; hs_stacks_pi() / gene_div_2n_counts() reproduce
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
#  size and by analyzing differentiated populations together, and that the
#  unbiased quantity is AUTOSOMAL heterozygosity: the same sum divided by every
#  sequenced site, monomorphic ones included. Retaining sites that are
#  monomorphic WITHIN a population but variable in another is a step toward
#  that, not the thing itself.
#
#  The conversion is one multiplication. For each population separately:
#
#       He_autosomal = He_per_SNP * (n_variant_records / n_sites_sequenced)
#
#  He_per_SNP * n_variant_records rebuilds He summed over the population's
#  SNPs, and dividing by n_sites_sequenced spreads it over every sequenced
#  site (fixed sites add 0 to the sum).
#
#  n_sites_sequenced: every sequenced site, variant or fixed, of the RAD loci
#  in the data, counting a site for a population only where it has at least
#  one genotyped individual. For data exactly as Stacks wrote them this is
#  the population's `Sites` in the "All positions (variant and fixed)" block
#  of populations.sumstats_summary.tsv (SumStatsSummary::accumulate(), Stacks
#  2.68). If whole loci are removed afterwards, their sites must leave this
#  count too, and Stacks' `Sites` no longer fits.
#
#  n_variant_records: the SNP records IN THE DATA where that population has
#  at least one genotyped individual. Not every record in the VCF: with 3 or
#  more populations and `-p` below their number, Stacks blanks a population
#  at the sites where it fails `-r`, and those sites are not in its `Sites`.
#  Not just the records that passed diversity_stats()'s min_n rule either:
#  the mean over the used records stands in for the mean over all of them.
#  And not Stacks' `Variant_Sites`: that also counts SNPs removed after Stacks,
#  and multiplying by it would treat each removed SNP as having the average He
#  of the kept ones -- far too high after a MAC or MAF filter, which removes
#  low-He SNPs. diversity_stats() counts it per population
#  (.autosomal_counts() in R/diversity_stats.R).
#
#  Stacks reports both blocks; the "All positions" one is the
#  Schmidt-compliant estimate, and the "Variant positions" one is not.
#
#  NO ALLELE-FREQUENCY FILTER. Stacks' --min-mac/--min-maf turn a failing SNP
#  into a fixed site that still counts in `Sites`, so its contribution to pi
#  is lost from Stacks' own all-positions Pi and from this conversion alike.
#  Under a neutral site-frequency spectrum, sites with minor allele count <= 2
#  carry about 4/(N - 1) of pi for N pooled gene copies: 21% at 10 diploids,
#  10% at 20, 7% at 30 -- a bias that depends on sample size.
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

#' Low-level estimators: gene diversity, FIS, per-site scaling, rarefaction
#'
#' The building blocks [diversity_stats()] is made of, exported so any number
#' it reports can be checked by hand. Rarefaction sizes are in GENE COPIES:
#' 10 diploid individuals are 20 gene copies.
#'
#' * `hs_nei_chesser()`: unbiased gene diversity (expected heterozygosity,
#'   Hs) for one population at one locus, `(n/(n-1)) * (1 - sum_p2 -
#'   ho/(2n))` (Nei & Chesser 1983). Unbiased at any FIS. Vectorised.
#' * `hs_biallelic()`: the same at a biallelic site, from one allele's
#'   frequency `p`.
#' * `hs_stacks_pi()`: `(2n/(2n-1)) * 2p(1-p)`, the estimator behind Stacks'
#'   `Pi` column. Unbiased only when FIS = 0 and biased low otherwise;
#'   provided so the two can be compared.
#' * `gene_div_2n_counts()`: the same correction as `hs_stacks_pi()`, from
#'   allele counts, so it also works at multi-allelic (haplotype) loci.
#' * `hs_from_counts()`: `hs_nei_chesser()` from allele counts.
#' * `fis_ratio_of_sums()`: `1 - sum(ho) / sum(he)` over loci -- a ratio of
#'   sums, never a mean of per-locus ratios (Weir & Cockerham 1984). A locus
#'   monomorphic in the population adds 0 to both sums, so including it
#'   changes nothing.
#' * `autosomal_het()`: heterozygosity per variant record rescaled to per
#'   sequenced site (Schmidt et al. 2021), `het_per_snp * n_snps_used /
#'   n_sites_sequenced`; `NA` wherever `n_sites_sequenced` is missing or
#'   smaller than `n_snps_used`. Vectorised, so each population can have its
#'   own `n_sites_sequenced`.
#' * `p_sampled()`: for each allele, the probability it appears at least once
#'   in `g` gene copies drawn without replacement.
#' * `rare_richness()`: rarefied allelic richness at one locus, the expected
#'   number of distinct alleles in `g` gene copies (Hurlbert 1971; El
#'   Mousadik & Petit 1996).
#' * `rare_private()`, `rare_private_all()`: rarefied private allelic
#'   richness, the expected number of alleles present in `g` copies from one
#'   population and absent from `g` copies of every other population
#'   (Kalinowski 2004; Szpiech et al. 2008).
#'
#' The derivations are in the header of `R/estimators.R`;
#' [diversity_core_selftest()] checks each estimator against brute-force
#' Monte Carlo and against published Stacks output.
#'
#' @param sum_p2 Sum of squared sample allele frequencies.
#' @param ho Observed heterozygote frequency. For `fis_ratio_of_sums()`, a
#'   vector of per-locus (or per-site) values.
#' @param n Number of diploid individuals genotyped.
#' @param p Frequency of one allele at a biallelic site.
#' @param counts Gene-copy counts, one per allele, for one population at one
#'   locus.
#' @param he Per-locus (or per-site) expected heterozygosity, same loci as
#'   `ho`.
#' @param het_per_snp Heterozygosity per variant record.
#' @param n_snps_used Number of SNP records in the data at which the
#'   population has at least one genotyped individual (not only those used to
#'   estimate `het_per_snp`). One value, or one per element of `het_per_snp`.
#' @param n_sites_sequenced Sequenced sites, variant and fixed, of the RAD loci
#'   in the data, counting a site only where the population has at least one
#'   genotyped individual. For unfiltered Stacks output this is `Sites` in the
#'   "All positions (variant and fixed)" block of
#'   `populations.sumstats_summary.tsv`; see `sites` in [diversity_stats()]
#'   for filtered data. One value, or one per element of `het_per_snp`.
#' @param g Rarefaction size in gene copies.
#' @param count_mat Matrix of gene-copy counts, populations (rows) x alleles.
#' @param j Row (population) of `count_mat` to compute private richness for.
#' @return A numeric vector: one value per input element for the vectorised
#'   `hs_nei_chesser()`, `hs_biallelic()`, `hs_stacks_pi()` and
#'   `autosomal_het()`; one per allele for `p_sampled()`; one per population
#'   for `rare_private_all()`; otherwise a single value. `NA` wherever the
#'   quantity is undefined (fewer than 2 individuals, or fewer than `g` gene
#'   copies).
#' @references
#' Nei, M. & Chesser, R.K. (1983) Estimation of fixation indices and gene
#' diversities. *Annals of Human Genetics* 47:253-259.
#'
#' Weir, B.S. & Cockerham, C.C. (1984) Estimating F-statistics for the
#' analysis of population structure. *Evolution* 38:1358-1370.
#'
#' Hurlbert, S.H. (1971) The nonconcept of species diversity: a critique and
#' alternative parameters. *Ecology* 52:577-586.
#'
#' El Mousadik, A. & Petit, R.J. (1996) High level of genetic differentiation
#' for allelic richness among populations of the argan tree. *Theoretical and
#' Applied Genetics* 92:832-839.
#'
#' Kalinowski, S.T. (2004) Counting alleles with rarefaction: private alleles
#' and hierarchical sampling designs. *Conservation Genetics* 5:539-543.
#'
#' Szpiech, Z.A., Jakobsson, M. & Rosenberg, N.A. (2008) ADZE: a rarefaction
#' approach for counting alleles private to combinations of populations.
#' *Bioinformatics* 24:2498-2504.
#'
#' Schmidt, T.L., Jasper, M.-E., Weeks, A.R. & Hoffmann, A.A. (2021) Unbiased
#' population heterozygosity estimates from genome-wide sequence data.
#' *Methods in Ecology and Evolution* 12:1888-1898.
#' @examples
#' # One biallelic locus: 10 diploids, allele frequency 0.3, Ho = 0.35.
#' hs_biallelic(0.3, ho = 0.35, n = 10)    # Nei & Chesser
#' hs_stacks_pi(0.3, n = 10)               # what Stacks' `Pi` column would give
#'
#' # The same locus from allele counts (6 and 14 gene copies):
#' hs_from_counts(c(6, 14), ho = 0.35, n = 10)
#' gene_div_2n_counts(c(6, 14))
#'
#' # FIS as a ratio of sums over three loci (the third is monomorphic):
#' fis_ratio_of_sums(ho = c(0.30, 0.10, 0), he = c(0.42, 0.15, 0))
#'
#' # Rarefied allelic richness, 4 alleles, rarefied to 10 gene copies:
#' rare_richness(c(12, 5, 2, 1), g = 10)
#' # Rarefied private allelic richness, two populations at one locus:
#' rare_private_all(rbind(popA = c(10, 5, 5, 0), popB = c(12, 8, 0, 0)), g = 10)
#'
#' # He per variant record -> per sequenced site:
#' autosomal_het(0.29, n_snps_used = 10348, n_sites_sequenced = 1e6)
#' @name estimators
NULL

#' @rdname estimators
#' @export
hs_nei_chesser <- function(sum_p2, ho, n) {
  out <- (n / (n - 1)) * (1 - sum_p2 - ho / (2 * n))
  out[!is.finite(out) | n < 2] <- NA_real_
  out
}

#' @rdname estimators
#' @export
hs_biallelic <- function(p, ho, n) hs_nei_chesser(p^2 + (1 - p)^2, ho, n)

#' @rdname estimators
#' @export
hs_stacks_pi <- function(p, n) {
  out <- 2 * p * (1 - p) * (2 * n) / (2 * n - 1)
  out[!is.finite(out) | n < 1] <- NA_real_
  out
}

## Named for the correction it applies, not for the software column it happens
## to match. N = gene copies actually observed at the locus (2 x the typed
## individuals), so this equals hs_stacks_pi() at that many individuals
## (checked in diversity_core_selftest()).
#' @rdname estimators
#' @export
gene_div_2n_counts <- function(counts) {
  N <- sum(counts)
  if (!is.finite(N) || N < 2) return(NA_real_)
  (N / (N - 1)) * (1 - sum((counts / N)^2))
}

#' @rdname estimators
#' @export
hs_from_counts <- function(counts, ho, n) {
  tot <- sum(counts)
  if (!is.finite(tot) || tot < 2 || n < 2) return(NA_real_)
  hs_nei_chesser(sum((counts / tot)^2), ho, n)
}

#' @rdname estimators
#' @export
fis_ratio_of_sums <- function(ho, he) {
  ok <- is.finite(ho) & is.finite(he)
  sum_he <- sum(he[ok])
  ## No expected heterozygosity at all (every locus monomorphic): undefined.
  if (!isTRUE(sum_he > .zero_tol)) return(NA_real_)
  1 - sum(ho[ok]) / sum_he
}

#' @rdname estimators
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

## Not exported. The rarefaction formula itself, for many records at once:
## `counts` is a records x alleles matrix of gene-copy counts, and the result
## has the same shape, holding for each allele the probability that it
## appears at least once in `g` gene copies drawn without replacement from
## that record's N copies:
##
##     1 - C(N - N_i, g) / C(N, g)
##
## computed on the log scale (lchoose) so it stays finite at large N. A row
## with fewer than `g` copies is all NA (the draw is impossible). An allele
## with count 0 gets exactly 0, so zero-padded columns change nothing.
## p_sampled(), rare_richness(), rare_private() and diversity_stats() all use
## this one function, so the formula exists in exactly one place.
.p_sampled_mat <- function(counts, g) {
  N <- rowSums(counts)
  out <- 1 - exp(lchoose(N - counts, g) - lchoose(N, g))
  out[!is.finite(N) | N < g, ] <- NA_real_
  out
}

#' @rdname estimators
#' @export
p_sampled <- function(counts, g) {
  if (!is.finite(g) || g < 1) return(rep(NA_real_, length(counts)))
  .p_sampled_mat(matrix(counts, nrow = 1L), g)[1L, ]
}

#' @rdname estimators
#' @export
rare_richness <- function(counts, g) {
  if (!is.finite(sum(counts)) || sum(counts) < g) return(NA_real_)
  p_drawn <- p_sampled(counts, g)
  if (anyNA(p_drawn)) NA_real_ else sum(p_drawn)
}

#' @rdname estimators
#' @export
rare_private <- function(count_mat, j, g) {
  if (!is.matrix(count_mat)) count_mat <- rbind(count_mat)
  ## Row = population, column = allele: Pr(allele drawn in g copies).
  p_drawn <- .p_sampled_mat(count_mat, g)
  ## "Private" needs every population's sample defined at the same depth g.
  if (anyNA(p_drawn)) return(NA_real_)
  ## Drawn in population j AND not drawn in each other population.
  term <- p_drawn[j, ]
  for (k in setdiff(seq_len(nrow(count_mat)), j)) term <- term * (1 - p_drawn[k, ])
  sum(term)
}

#' @rdname estimators
#' @export
rare_private_all <- function(count_mat, g) {
  vapply(seq_len(nrow(count_mat)), function(j) rare_private(count_mat, j, g),
         numeric(1))
}
