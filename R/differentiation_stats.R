###############################################################################
#
#  R/differentiation_stats.R -- how much do these populations DIFFER from
#  each other (as opposed to diversity_stats(), which asks how variable each
#  population is on its own)?
#
#  Two numbers answer that question, and this file computes both:
#
#    FST   Weir & Cockerham (1984), the most widely reported measure. It
#          ranges from about 0 (no differentiation) to 1 (no shared alleles).
#          Its ceiling shrinks as loci get more variable, which can make FST
#          look small on very diverse markers such as RAD haplotypes. That is
#          a property of FST itself, not of this implementation.
#    D     Jost's D (Jost 2008, Molecular Ecology 17:4015), built to avoid
#          that ceiling, and the better choice for highly variable markers.
#
#  When hierfstat is installed, Weir & Goudet's (2017) beta is also reported
#  (hierfstat::pairwise.betas()), alongside FST and D.
#
#  ENGINE AND VALIDATION. FST and FIS come from this package's own
#  implementation of Weir & Cockerham's (1984) variance components a, b and c
#  (also in Weir 1996, Genetic Data Analysis II), NOT from hierfstat: the
#  bootstrap and jackknife re-sum per-locus components thousands of times,
#  and hierfstat has no function that accepts a resampled locus set.
#  tests/testthat/test-differentiation-stats.R compares this implementation
#  with hierfstat::wc() and hierfstat::pairwise.WCfst(); on records typed in
#  at least two of the compared populations (with more than one individual
#  per population on average) they agree to machine precision. They differ,
#  deliberately, at records typed in only ONE of the compared populations,
#  or with a single individual in each typed population: this package skips
#  those for FST (VCFtools and Stacks' own pairwise Fst also skip records
#  typed in one population), while hierfstat keeps their within-population
#  variance and so pulls FST toward 0 (see fst_block() below). hierfstat is
#  used at run time only for beta and, with hierfstat_check = TRUE, a printed
#  comparison of the global FST.
#
#  Jost's D follows the five-step formula of the mmod package's HsHt()
#  function (Winter 2012; MIT license, see inst/NOTICE), rewritten for this
#  package's genotype matrices. Several loci are combined by averaging Hs and
#  Ht over loci FIRST and computing one D from the averages -- mmod's
#  `global.het` summary; mmod also offers the harmonic mean of per-locus D
#  (`global.harm_mean`), which is not computed here. Jost's bias-corrected Hs
#  uses the 2N/(2N - 1) factor, which assumes random mating within
#  populations (the same assumption behind Stacks' `Pi`; see R/estimators.R).
#  With a heterozygote deficit Hs is underestimated by about F/(2N - 1), so D
#  is biased upward, most noticeably when Ht - Hs is small.
#
#  WHICH RECORDS D USES. A pairwise D uses every record typed in both
#  populations of the pair. The global D uses only records typed in EVERY
#  population: a record missing some populations measures differentiation
#  among a different set of populations, and averaging it in biases the
#  global D (see d_block() below). The result table's `D_records` column and
#  print() say how many records that is.
#
#  PRECISION: a bootstrap over whole RAD loci (the locus_raw groups) for a 95%
#  interval, and a delete-one-locus jackknife for a standard error, both from
#  per-locus sums (R/resampling.R). hierfstat::boot.ppfst() is not used
#  because it resamples SNP rows independently, ignoring linkage within a
#  RAD tag.
#
#  REFERENCES
#    Weir, B.S. & Cockerham, C.C. (1984) Estimating F-statistics for the
#      analysis of population structure. Evolution 38: 1358-1370.
#    Weir, B.S. & Goudet, J. (2017) A unified characterization of population
#      structure and relatedness. Genetics 206: 2085-2103.
#    Jost, L. (2008) GST and its relatives do not measure differentiation.
#      Molecular Ecology 17: 4015-4026.
#    Winter, D.J. (2012) mmod: an R library for the calculation of
#      population differentiation statistics. Molecular Ecology Resources
#      12: 1158-1160.
#
###############################################################################

## ---------------------------------------------------------------------------
## Shared inputs: counts per population, computed once
## ---------------------------------------------------------------------------

## Not exported. The per-population counts every FST and D calculation needs,
## for every record of H (see .pop_counts() in R/vcf_io.R):
##   counts      list (one per population) of records x alleles gene-copy counts
##   het_counts  the same, counting only heterozygous individuals' copies: each
##               heterozygote adds one to each of its two alleles
##   n_typed     records x populations: genotyped individuals
.differentiation_counts <- function(H, pops) {
  k <- max(1L, H$n_alleles)
  het_counts <- lapply(pops, function(ids) {
    a1 <- H$A1[, ids, drop = FALSE]
    a2 <- H$A2[, ids, drop = FALSE]
    homozygous <- !is.na(a1) & a1 == a2
    a1[homozygous] <- NA_integer_
    a2[homozygous] <- NA_integer_
    .allele_counts(a1, a2, k)
  })
  list(counts = .pop_counts(H, pops, k = k), het_counts = het_counts,
       n_typed = .typed_by_pop(H, pops))
}

## ---------------------------------------------------------------------------
## Weir & Cockerham (1984)
## ---------------------------------------------------------------------------

## Not exported. For one set of populations (all of them for the global FST,
## or two for a pairwise FST), Weir & Cockerham's three variance components
## at every record:
##   comp_a  variance AMONG the populations                 ("a")
##   comp_b  variance among INDIVIDUALS within populations  ("b")
##   comp_c  variance WITHIN individuals (heterozygosity)   ("c")
## FST is the share of the total (a + b + c) that is among populations; see
## .fst_from_sums(). Components are kept per record so the jackknife and
## bootstrap can re-sum them without recounting genotypes.
##
## `dc` is .differentiation_counts() restricted to those populations. All
## quantities below are records x alleles matrices (or per-record vectors,
## which R recycles down each allele column); only populations with at least
## one genotyped individual at a record take part at that record.
.wc_components <- function(dc) {
  n_typed <- dc$n_typed                                   # records x populations
  has_data <- n_typed > 0
  n_pops_with_data <- rowSums(has_data)                   # r
  n_total <- rowSums(n_typed)                             # sum of n_i
  n_mean <- n_total / n_pops_with_data                    # n-bar

  p_bar <- Reduce(`+`, dc$counts) / (2 * n_total)         # pooled allele frequency
  h_bar <- Reduce(`+`, dc$het_counts) / n_total           # pooled heterozygote frequency

  ## s2: the sample-size-weighted variance of the populations' allele
  ## frequencies around p_bar -- how different the populations are.
  weighted_sq_dev <- 0
  for (p in seq_along(dc$counts)) {
    p_i <- dc$counts[[p]] / (2 * n_typed[, p])
    term <- n_typed[, p] * (p_i - p_bar)^2
    term[!has_data[, p], ] <- 0
    weighted_sq_dev <- weighted_sq_dev + term
  }
  s2 <- weighted_sq_dev / ((n_pops_with_data - 1) * n_mean)
  s2[!(n_pops_with_data > 1), ] <- 0

  ## n_c, the sample-size correction for unequal population sizes.
  n_c <- (n_total - rowSums(n_typed^2) / n_total) / (n_pops_with_data - 1)
  r_ratio <- (n_pops_with_data - 1) / n_pops_with_data

  ## a needs at least two populations with data and n-bar > 1; b needs
  ## n-bar > 1. A record that fails contributes 0 to that component (as
  ## hierfstat::wc() also treats it).
  comp_a <- (n_mean / n_c) *
    (s2 - (1 / (n_mean - 1)) * (p_bar * (1 - p_bar) - r_ratio * s2 - h_bar / 4))
  comp_a[!(n_pops_with_data > 1 & !is.na(n_mean) & n_mean > 1), ] <- 0
  ## NOTE: b uses n-bar - 1 (the plain mean sample size), NOT n_c - 1. This is
  ## the published estimator; the hierfstat comparison in the tests would
  ## catch a swap.
  comp_b <- (n_mean / (n_mean - 1)) *
    (p_bar * (1 - p_bar) - r_ratio * s2 - ((2 * n_mean - 1) / (4 * n_mean)) * h_bar)
  comp_b[!(!is.na(n_mean) & n_mean > 1), ] <- 0
  comp_c <- h_bar / 2

  ## Sum over alleles; a record with fewer than 2 genotyped individuals in
  ## total is undefined.
  undefined <- n_pops_with_data < 1 | n_total < 2
  totals <- list(comp_a = rowSums(comp_a), comp_b = rowSums(comp_b), comp_c = rowSums(comp_c))
  out <- lapply(totals, function(v) {
    v[undefined] <- NA_real_
    v
  })
  ## How many of the compared populations have data at each record, and their
  ## mean sample size: FST uses a record only where at least 2 populations
  ## have data and n-bar > 1 (see fst_block() in differentiation_stats()).
  out$n_pops_with_data <- n_pops_with_data
  out$n_mean <- n_mean
  out
}

## Not exported. FST and the "metapopulation" FIS from a, b and c SUMMED over
## the loci an estimate uses (all loci for the point estimate, a resampled set
## for a bootstrap replicate, all but one for a jackknife replicate). A ratio
## of sums, never a mean of per-locus ratios:
##   FST = sum(a) / sum(a + b + c)
##   FIS = sum(b) / sum(b + c)
## `n` is the number of loci with the components defined. Vectorised: one
## element per estimate. No defined locus, or a total of ~0 (every locus
## monomorphic), gives NA.
.fst_from_sums <- function(n, a_sum, b_sum, c_sum) {
  total <- a_sum + b_sum + c_sum
  ifelse(n < 0.5 | abs(total) <= .zero_tol, NA_real_, a_sum / total)
}
.fis_from_sums <- function(n, b_sum, c_sum) {
  total <- b_sum + c_sum
  ifelse(n < 0.5 | abs(total) <= .zero_tol, NA_real_, b_sum / total)
}

## ---------------------------------------------------------------------------
## Jost's D (2008)
## ---------------------------------------------------------------------------

## Not exported. For one set of populations, Jost's bias-corrected
## within-population (Hs_est) and total (Ht_est) gene diversity at every
## record: the two ingredients of D (see .jost_d_from_sums()). Unlike Weir &
## Cockerham, Jost's formula weights populations EQUALLY. The five steps
## (harmonic-mean n -> HpS -> Hs_est -> HpT -> Ht_est) follow the mmod
## package's HsHt() function (Winter 2012; see inst/NOTICE). A record with
## fewer than 2 populations genotyped is NA. `n_pops_with_data` is returned
## too: the global D uses only records where EVERY population has data (see
## d_block() in differentiation_stats()).
.jost_hsht <- function(dc) {
  n_typed <- dc$n_typed
  has_data <- n_typed > 0
  n_pops_with_data <- rowSums(has_data)

  ## Harmonic mean of population sizes (pulled down by any small population).
  harmonic_n <- n_pops_with_data / rowSums(ifelse(has_data, 1 / n_typed, 0))
  ## Naive within-population diversity, averaged over populations, and the
  ## populations' average allele frequencies.
  within_div <- 0
  mean_freq <- 0
  for (p in seq_along(dc$counts)) {
    freq <- dc$counts[[p]] / (2 * n_typed[, p])
    freq[!has_data[, p], ] <- 0
    within_div <- within_div + ifelse(has_data[, p], 1 - rowSums(freq^2), 0)
    mean_freq <- mean_freq + freq
  }
  Hp_S <- within_div / n_pops_with_data
  mean_freq <- mean_freq / n_pops_with_data

  Hs_est <- (2 * harmonic_n / (2 * harmonic_n - 1)) * Hp_S
  Hp_T <- 1 - rowSums(mean_freq^2)
  Ht_est <- Hp_T + Hs_est / (2 * harmonic_n * n_pops_with_data)
  too_few <- n_pops_with_data < 2
  Hs_est[too_few] <- NA_real_
  Ht_est[too_few] <- NA_real_
  list(Hs = Hs_est, Ht = Ht_est, n_pops_with_data = n_pops_with_data)
}

## Not exported. Jost's D from Hs_est and Ht_est SUMMED over `n` loci: average
## Hs and Ht over loci first, then compute one D. `n_pops` is the FIXED
## number of populations in the comparison (2 for a pair), not how many had
## data at a particular record. Vectorised over estimates.
.jost_d_from_sums <- function(n, hs_sum, ht_sum, n_pops) {
  hs_mean <- ifelse(n < 0.5, NA_real_, hs_sum / n)
  ht_mean <- ht_sum / n
  ifelse(is.na(hs_mean) | hs_mean >= 1, NA_real_,
         (ht_mean - hs_mean) / (1 - hs_mean) * (n_pops / (n_pops - 1)))
}

## ---------------------------------------------------------------------------
## The exported function
## ---------------------------------------------------------------------------

#' Population differentiation: FST, Weir & Goudet's beta, and Jost's D
#'
#' Computes how much populations differ from each other (as opposed to
#' [diversity_stats()], which asks how variable each population is). Reports
#' Weir & Cockerham's (1984) FST and FIS and Jost's (2008) D, for all
#' populations together (global) and for each pair, and, when the `hierfstat`
#' package is installed, Weir & Goudet's (2017) pairwise beta. Printing the
#' result shows the tables; `summary()` adds notes on interpretation.
#'
#' @details
#' **Why report both FST and D.** FST's ceiling shrinks as a locus's own
#' diversity grows, so two datasets with identical population structure can
#' show very different FST just because one used more variable markers (Jost
#' 2008). RAD haplotype loci are such markers, so on a haplotype VCF prefer D
#' (or beta) for comparing datasets.
#'
#' **Bootstrap and jackknife.** Confidence intervals come from this package's
#' bootstrap over whole RAD loci, not from `hierfstat::boot.ppfst()`, which
#' resamples SNP rows independently and so ignores the linkage between SNPs
#' on one RAD tag. The FST/FIS formula is Weir & Cockerham's (1984) published
#' estimator; the package tests check it against `hierfstat::wc()`.
#'
#' **Records used for FST.** FST (global and pairwise) and the global FIS use a
#' record only where at least two of the compared populations have a
#' genotyped individual. A record typed in only one of them carries no
#' information about differences between populations, but including it
#' would add its within-population variance to FST's denominator and pull
#' FST toward 0. Such records are common when a Stacks run has 3 or more
#' populations and `-p` below their number. VCFtools and Stacks' own
#' pairwise FST skip them too; `hierfstat::wc()` and `pairwise.WCfst()` do
#' not, so on such data their FST is lower than this function's. Records
#' where every typed population has a single genotyped individual are also
#' skipped for FST and FIS, because no variance among individuals can be
#' estimated there. Missing individuals within a typed population are
#' handled by the estimator itself and are not affected. The number of
#' skipped records is in `settings$fst_records_skipped`.
#'
#' **Records used for global D.** The global D uses ONLY records typed in
#' every population; each pairwise D uses the records typed in both
#' populations of its pair. Jost's D for k populations measures
#' differentiation among those k populations, and a record missing some of
#' them measures it among a different set: averaging such records in biased
#' the global D low (0.875 instead of 1 for three completely differentiated
#' populations, one of them untyped at half the records). The global table's
#' `D_records` column, `print()` and `summary()` give the number of records
#' used. With 3 or more populations and much missing data this can be a
#' small share of the records; report pairwise D as well. If no record is
#' typed in every population, the global D is `NA` with a warning.
#'
#' **Combining loci for D.** Hs and Ht are averaged over loci and one D is
#' computed from the averages (mmod's `global.het`), not averaged over
#' per-locus D values. Jost's bias correction of Hs, 2N/(2N - 1), assumes
#' random mating within each population (the same assumption behind Stacks'
#' `Pi`), so with a heterozygote deficit Hs is slightly underestimated. D is
#' then biased upward, most noticeably when differentiation is weak.
#'
#' **References.** Weir, B.S. & Cockerham, C.C. (1984) Estimating
#' F-statistics for the analysis of population structure. *Evolution*
#' 38:1358-1370. \doi{10.1111/j.1558-5646.1984.tb05657.x} --
#' Weir, B.S. & Goudet, J. (2017) A unified characterization of population
#' structure and relatedness. *Genetics* 206:2085-2103.
#' \doi{10.1534/genetics.116.198424} -- Jost, L. (2008) GST and
#' its relatives do not measure differentiation. *Molecular Ecology*
#' 17:4015-4026. \doi{10.1111/j.1365-294X.2008.03887.x} --
#' Winter, D.J. (2012) mmod: an R library for the calculation of population
#' differentiation statistics. *Molecular Ecology Resources* 12:1158-1160.
#' \doi{10.1111/j.1755-0998.2012.03174.x} --
#' Goudet, J. (2005) HIERFSTAT, a package for R to compute and test
#' hierarchical F-statistics. *Molecular Ecology Notes* 5:184-186.
#'
#' @inheritParams diversity_stats
#' @param nboot Bootstrap replicates over RAD loci. Default `10000`; `0` for
#'   none (the jackknife standard errors are always computed).
#' @param beta Also compute Weir & Goudet's pairwise beta with
#'   `hierfstat::pairwise.betas()`, when hierfstat is installed. Default
#'   `TRUE`. On large datasets this one call takes most of the run time (about
#'   95% at 100,000 SNPs); set `FALSE` to skip it. FST and D do not need it.
#' @param hierfstat_check If `TRUE` and `hierfstat` is installed, also run
#'   `hierfstat::wc()` and report how far its global FST is from this
#'   package's. Default `FALSE` (the package tests make this comparison, and
#'   `wc()` is slow on large datasets).
#' @return An object of class `raddiv_differentiation`, a list of:
#'   \describe{
#'     \item{global}{One row: `FST`, `FIS` and `D`, each with jackknife SE
#'       (`_se`) and 95% bootstrap interval (`_lo`, `_hi`), and `D_records`,
#'       the number of records the global D uses: only those typed in every
#'       population (see Details).}
#'     \item{pairwise}{One row per pair of populations: FST, beta and D, with
#'       standard errors and bootstrap intervals. Each pair's D uses the
#'       records typed in both of its populations.}
#'     \item{pairwise_fst, pairwise_beta, pairwise_D}{The pairwise point
#'       estimates as symmetric matrices (`NA` on the diagonal).
#'       `pairwise_beta` is `NULL` when hierfstat is not installed or
#'       `beta = FALSE`.}
#'     \item{settings}{The settings of this run.}
#'   }
#'   Values are stored at full precision. With `outdir`, `global` and
#'   `pairwise` are written to `differentiation_global.<stem>.tsv` and
#'   `differentiation_pairwise.<stem>.tsv`.
#' @examples
#' vcf    <- system.file("extdata", "small.haps.vcf", package = "RADdiversity")
#' popmap <- system.file("extdata", "small_popmap.tsv", package = "RADdiversity")
#' set.seed(1)
#' res <- differentiation_stats(vcf, popmap, nboot = 100, verbose = FALSE)
#' res
#' res$global$D
#' @export
differentiation_stats <- function(vcf, popmap, nboot = 10000L, beta = TRUE,
                                  hierfstat_check = FALSE, outdir = NULL, stem = NULL,
                                  seed = NULL, verbose = TRUE) {
  ## ---- 1. Check the arguments ----------------------------------------------
  nboot <- .check_count(nboot, "nboot")
  .check_flag(beta, "beta")
  .check_flag(hierfstat_check, "hierfstat_check")
  .check_flag(verbose, "verbose")
  .check_seed(seed)
  .check_run_inputs(vcf, popmap, stem, outdir)
  if (!is.null(seed)) {
    restore_rng <- .save_rng_state()
    on.exit(restore_rng(), add = TRUE)
    set.seed(seed)
  }

  ## ---- 2. Read the data -----------------------------------------------------
  H <- .resolve_H(vcf, verbose = verbose)
  pops <- .resolve_pops(popmap, H$samples, verbose = verbose)
  pop_names <- names(pops)
  n_pops <- length(pops)
  if (n_pops < 2) stop("Need at least 2 populations.", call. = FALSE)
  tiny <- pop_names[lengths(pops) < 2]
  if (length(tiny))
    stop("Population(s) with fewer than 2 individuals: ", paste(tiny, collapse = ", "),
         "\n  FST/D need at least 2 individuals per population to define a ",
         "within-population allele frequency.", call. = FALSE)

  n_rec <- nrow(H$A1)
  locus_index <- match(H$locus_raw, unique(H$locus_raw))
  n_loci <- max(locus_index)
  .inform(verbose, sprintf("  %s records on %s RAD loci (%.2f per locus)",
                           .big(n_rec), .big(n_loci), n_rec / n_loci))
  pairs <- utils::combn(pop_names, 2, simplify = FALSE)
  pair_labels <- make.unique(vapply(pairs, paste, character(1), collapse = "__"))
  ## Statistics are looked up by the pair's NUMBER, not its label: two pairs
  ## can share a label (populations "a", "a__b", "b__c" and "c" give
  ## "a__b__c" twice), and a lookup by label would return one pair's values
  ## for both.
  pair_keys <- paste0("pair", seq_along(pairs))
  .inform(verbose, sprintf("Computing FST/D for %d population(s), %d pair(s) ...",
                           n_pops, length(pairs)))

  ## ---- 3. Per-record components ---------------------------------------------
  ## Counts are made once; the global comparison uses every population and
  ## each pair uses its two. Components do not change between resampling
  ## replicates -- only which records are summed does.
  dc <- .differentiation_counts(H, pops)
  subset_counts <- function(pair) {
    list(counts = dc$counts[pair], het_counts = dc$het_counts[pair],
         n_typed = dc$n_typed[, pair, drop = FALSE])
  }
  wc_global <- .wc_components(dc)
  jost_global <- .jost_hsht(dc)
  wc_pairs <- lapply(pairs, function(pair) .wc_components(subset_counts(pair)))
  jost_pairs <- lapply(pairs, function(pair) .jost_hsht(subset_counts(pair)))

  ## ---- 4. Pieces summed per RAD locus ---------------------------------------
  ## Each comparison contributes a small block of columns. A record where a
  ## component is undefined adds 0 to every sum and to the count `n`.
  ##   FST block (4 columns): n, a, b, c   (records where a, b, c are defined,
  ##                                        >= 2 compared populations are typed,
  ##                                        and n-bar > 1)
  ##   FIS block (3 columns): n, b, c      (global only; where b, c are defined
  ##                                        and n-bar > 1)
  ##   D block   (3 columns): n, Hs, Ht    (pairwise: records where Hs and Ht
  ##                                        are defined; global: only records
  ##                                        typed in EVERY population)
  ##
  ## WHY FST SKIPS A RECORD WHERE ONLY ONE COMPARED POPULATION IS TYPED. Such a
  ## record says nothing about differences between populations (a = 0), but
  ## its within-population variance (b + c) would still be added to the
  ## denominator, pulling FST toward 0 in proportion to how many such records
  ## there are. They are common: with 3 or more populations and Stacks' -p
  ## below their number, a population failing -r at a site is blanked there
  ## while the site is kept for the others. hierfstat::wc() and
  ## pairwise.WCfst() keep b and c at these records (and so dilute FST);
  ## VCFtools --weir-fst-pop and Stacks' own pairwise Fst skip them, as this
  ## package does. Ordinary missing data -- some individuals missing, both
  ## populations typed -- is unaffected: Weir & Cockerham's unequal-sample-size
  ## terms handle it.
  ##
  ## WHY FST AND FIS ALSO SKIP A RECORD WITH n-bar <= 1 (every typed population
  ## has exactly one genotyped individual). There, b (variance among
  ## individuals within populations) cannot be estimated and is set to 0, and
  ## a is 0, but c (heterozygosity) is still added, pulling FST and FIS toward
  ## 0. hierfstat::wc() keeps such records (dropping its own undefined a and
  ## b there but keeping c); this package skips them, since no variance among
  ## individuals can be estimated from one individual per population.
  informative <- function(w) w$n_pops_with_data >= 2 & !is.na(w$n_mean) & w$n_mean > 1
  fst_block <- function(w) {
    ok <- is.finite(w$comp_a) & is.finite(w$comp_b) & is.finite(w$comp_c) & informative(w)
    cbind(ok, ifelse(ok, w$comp_a, 0), ifelse(ok, w$comp_b, 0), ifelse(ok, w$comp_c, 0))
  }
  ## Records with genotypes that FST skips, per comparison (reported in
  ## settings and summary()).
  fst_skipped <- function(w) sum(is.finite(w$comp_c) & !informative(w))
  fis_block <- function(w) {
    ok <- is.finite(w$comp_b) & is.finite(w$comp_c) & !is.na(w$n_mean) & w$n_mean > 1
    cbind(ok, ifelse(ok, w$comp_b, 0), ifelse(ok, w$comp_c, 0))
  }
  ## WHY THE GLOBAL D USES ONLY RECORDS TYPED IN EVERY POPULATION. Jost's D
  ## for k populations measures differentiation among those k populations,
  ## and its final step multiplies by k/(k - 1). A record typed in only k' < k
  ## of them gives Hs and Ht for a different, smaller set of populations, so
  ## averaging it in mixes quantities: with three completely differentiated
  ## populations (true D = 1) and one of them untyped at half the records,
  ## global D came out as 0.875. Requiring all k populations makes every
  ## record estimate the same quantity. A pairwise D (k = 2) already requires
  ## both populations, so it is unaffected.
  d_block <- function(j, need_all = NULL) {
    ok <- is.finite(j$Hs) & is.finite(j$Ht)
    if (!is.null(need_all)) ok <- ok & j$n_pops_with_data == need_all
    cbind(ok, ifelse(ok, j$Hs, 0), ifelse(ok, j$Ht, 0))
  }
  d_global_block <- d_block(jost_global, need_all = n_pops)
  d_records_global <- as.integer(sum(d_global_block[, 1]))
  blocks <- c(list(fst_block(wc_global), fis_block(wc_global), d_global_block),
              lapply(wc_pairs, fst_block), lapply(jost_pairs, d_block))
  ## Where each block's columns sit in the combined matrix.
  block_end <- cumsum(vapply(blocks, ncol, integer(1)))
  block_columns <- Map(function(end, width) (end - width + 1L):end,
                       block_end, vapply(blocks, ncol, integer(1)))
  n_pairs <- length(pairs)
  columns <- list(fst_global = block_columns[[1]], fis_global = block_columns[[2]],
                  d_global = block_columns[[3]],
                  fst_pair = block_columns[3L + seq_len(n_pairs)],
                  d_pair = block_columns[3L + n_pairs + seq_len(n_pairs)])
  pieces <- do.call(cbind, blocks)
  storage.mode(pieces) <- "double"
  S <- rowsum(pieces, locus_index, reorder = TRUE)        # RAD loci x pieces

  stat_names <- c("FST", "FIS", "D", paste0("pFST_", pair_keys), paste0("pD_", pair_keys))
  stats_from_sums <- function(totals) {
    totals <- matrix(totals, ncol = ncol(S))
    col <- function(j) totals[, j]
    fst_of <- function(cols) .fst_from_sums(col(cols[1]), col(cols[2]), col(cols[3]), col(cols[4]))
    d_of <- function(cols, k) .jost_d_from_sums(col(cols[1]), col(cols[2]), col(cols[3]), k)
    per_pair <- function(f) matrix(vapply(seq_len(n_pairs), f, numeric(nrow(totals))),
                                   nrow = nrow(totals))
    out <- cbind(fst_of(columns$fst_global),
                 .fis_from_sums(col(columns$fis_global[1]), col(columns$fis_global[2]),
                                col(columns$fis_global[3])),
                 d_of(columns$d_global, n_pops),
                 per_pair(function(k) fst_of(columns$fst_pair[[k]])),
                 per_pair(function(k) d_of(columns$d_pair[[k]], 2)))
    colnames(out) <- stat_names
    out
  }

  ## ---- 5. Estimates, jackknife and bootstrap --------------------------------
  point <- stats_from_sums(colSums(S))[1L, ]
  se <- .jack_block_sums(S, stats_from_sums)
  boot_matrix <- NULL
  if (nboot > 0) {
    .inform(verbose, sprintf("Bootstrapping %s replicates over %s RAD loci ...",
                             .big(nboot), .big(n_loci)))
    boot_matrix <- .boot_block_sums(S, nboot, stats_from_sums)
  }
  ci <- .percentile_ci(boot_matrix, stat_names)

  ## ---- 6. hierfstat: beta, and the optional cross-check ---------------------
  beta_matrix <- NULL
  if (.hierfstat_available() && (beta || hierfstat_check)) {
    if (hierfstat_check) {
      ## Compared on the records this package's global FST uses (>= 2
      ## populations typed, n-bar > 1); see "Records used for FST" in
      ## ?differentiation_stats for why the others differ.
      fst_rows <- which(informative(wc_global))
      dat_fst <- .to_hierfstat_df(H, pops, fst_rows)
      wc <- if (is.null(dat_fst)) simpleError("a record has more than 99 alleles, which hierfstat cannot read reliably")
            else tryCatch(hierfstat::wc(dat_fst), error = function(e) e)
      if (inherits(wc, "error")) {
        warning("Skipped the hierfstat::wc() cross-check: it failed with \"",
                conditionMessage(wc), "\".", call. = FALSE)
      } else {
        d_fst <- abs(wc$FST - point[["FST"]])
        .inform(verbose, sprintf("  cross-check vs hierfstat::wc(): FST diff %.2e", d_fst))
        if (d_fst > 1e-6) warning("This run's FST disagrees with hierfstat::wc().", call. = FALSE)
      }
    }
    if (beta) {
      .inform(verbose, "Computing Weir & Goudet's beta with hierfstat::pairwise.betas() ",
              "(slow on large datasets; beta = FALSE skips it) ...")
      dat_beta <- .to_hierfstat_df(H, pops)
      beta_matrix <- if (is.null(dat_beta)) NULL
                     else tryCatch(as.matrix(hierfstat::pairwise.betas(dat_beta)), error = function(e) NULL)
      if (is.null(dat_beta)) {
        .inform(verbose, "  beta skipped: a record has more than 99 alleles, which hierfstat cannot ",
                "read reliably; pairwise_beta will be NULL.")
      } else if (is.null(beta_matrix)) {
        .inform(verbose, "  hierfstat::pairwise.betas() failed on this dataset; pairwise_beta will be NULL.")
      } else {
        dimnames(beta_matrix) <- list(pop_names, pop_names)
      }
    }
  } else if (beta) {
    .inform(verbose, "Install hierfstat for Weir & Goudet's beta: install.packages(\"hierfstat\")")
  }

  ## ---- 7. Result tables ------------------------------------------------------
  with_uncertainty <- function(name, label) {
    stats::setNames(data.frame(point[[name]], se[[name]], ci[name, "lo"], ci[name, "hi"]),
                    paste0(label, c("", "_se", "_lo", "_hi")))
  }
  ## D_records: how many records the global D rests on (those typed in every
  ## population), so the table itself says so.
  global <- cbind(with_uncertainty("FST", "FST"), with_uncertainty("FIS", "FIS"),
                  with_uncertainty("D", "D"), D_records = d_records_global)
  pairwise <- do.call(rbind, lapply(seq_len(n_pairs), function(k) {
    p1 <- pairs[[k]][1]
    p2 <- pairs[[k]][2]
    cbind(data.frame(pop1 = p1, pop2 = p2),
          with_uncertainty(paste0("pFST_", pair_keys[k]), "FST"),
          data.frame(beta = if (is.null(beta_matrix)) NA_real_ else beta_matrix[p1, p2]),
          with_uncertainty(paste0("pD_", pair_keys[k]), "D"))
  }))
  matrix_of <- function(prefix) {
    m <- matrix(NA_real_, n_pops, n_pops, dimnames = list(pop_names, pop_names))
    for (k in seq_len(n_pairs)) {
      m[pairs[[k]][1], pairs[[k]][2]] <- point[[paste0(prefix, pair_keys[k])]]
      m[pairs[[k]][2], pairs[[k]][1]] <- point[[paste0(prefix, pair_keys[k])]]
    }
    m
  }

  written <- character(0)
  if (!is.null(outdir)) {
    stem <- .derive_stem(vcf, stem)
    written <- .write_tables(
      list(global = .round_table(global), pairwise = .round_table(pairwise)), outdir,
      c(global = sprintf("differentiation_global.%s.tsv", stem),
        pairwise = sprintf("differentiation_pairwise.%s.tsv", stem)))
  }

  fst_records_skipped <- stats::setNames(
    c(fst_skipped(wc_global), vapply(wc_pairs, fst_skipped, integer(1))),
    c("global", pair_labels))
  if (any(fst_records_skipped > 0))
    .inform(verbose, sprintf("  FST skipped records typed in fewer than 2 compared populations or with one individual per population: %s",
                             paste(sprintf("%s %s", names(fst_records_skipped)[fst_records_skipped > 0],
                                           .big(fst_records_skipped[fst_records_skipped > 0])),
                                   collapse = ", ")))
  .inform(verbose, "  ", .d_records_sentence(d_records_global, n_rec, n_pops))
  if (d_records_global == 0L)
    warning("Global D is NA: no record is genotyped in all ", n_pops, " populations, and the ",
            "global D uses only such records (see ?differentiation_stats, \"Records used ",
            "for global D\"). Pairwise D in $pairwise uses every record typed in both ",
            "populations of a pair.", call. = FALSE)
  settings <- list(n_pops = n_pops, n_records = n_rec, n_loci = n_loci, nboot = nboot,
                   is_haplotype = .is_haplotype_H(H), fst_records_skipped = fst_records_skipped,
                   d_records_global = d_records_global,
                   d_records_skipped_global = n_rec - d_records_global,
                   beta = beta, seed = seed, files = written)
  structure(list(global = global, pairwise = pairwise, pairwise_fst = matrix_of("pFST_"),
                 pairwise_beta = beta_matrix, pairwise_D = matrix_of("pD_"), settings = settings),
            class = "raddiv_differentiation")
}

## Not exported. The sentence that tells the user which records the global D
## rests on, used in the progress messages, print() and summary().
.d_records_sentence <- function(d_records, n_records, n_pops) {
  sprintf(paste("global D uses only the %s of %s records typed in %s;",
                "pairwise D uses the records typed in both populations of the pair"),
          .big(d_records), .big(n_records),
          if (n_pops == 2) "both populations" else sprintf("all %d populations", n_pops))
}

## Not exported. TRUE when the global D rests on fewer than half the records.
.d_records_few <- function(st) {
  !is.null(st$d_records_global) && st$d_records_global < 0.5 * st$n_records
}

#' @rdname differentiation_stats
#' @param x,object A `raddiv_differentiation` object, as returned by
#'   `differentiation_stats()`.
#' @param ... Ignored.
#' @export
print.raddiv_differentiation <- function(x, ...) {
  st <- x$settings
  cat(sprintf("Differentiation: %d populations, %s records on %s RAD loci (%s VCF)\n",
              st$n_pops, .big(st$n_records), .big(st$n_loci),
              if (st$is_haplotype) "haplotype" else "SNP"))
  cat("  FST, FIS: Weir & Cockerham (1984); D: Jost (2008); _se: jackknife over RAD loci",
      if (st$nboot > 0) sprintf("; _lo/_hi: 95%% bootstrap interval, %s replicates", .big(st$nboot)),
      "\n", sep = "")
  if (!is.null(st$d_records_global))
    cat("  ", .d_records_sentence(st$d_records_global, st$n_records, st$n_pops), "\n", sep = "")
  cat("\n$global\n")
  .print_table(.round_table(x$global))
  if (.d_records_few(st))
    cat("NOTE: global D rests on fewer than half the records (only those typed in every population); see summary().\n")
  cat("\n$pairwise\n")
  if (nrow(x$pairwise) <= 10) {
    .print_table(.round_table(x$pairwise))
  } else {
    cat(sprintf("  (first 10 of %d pairs; the full table is x$pairwise)\n", nrow(x$pairwise)))
    .print_table(.round_table(utils::head(x$pairwise, 10)))
  }
  cat("\n")
  if (st$is_haplotype) cat("NOTE: haplotype VCF -- prefer D (or beta) over FST (see summary()).\n")
  if (is.null(x$pairwise_beta))
    cat(if (isFALSE(st$beta)) "NOTE: beta not computed (beta = FALSE).\n"
        else "NOTE: beta needs the hierfstat package.\n")
  cat("summary() prints the full report.\n")
  invisible(x)
}

#' @rdname differentiation_stats
#' @export
summary.raddiv_differentiation <- function(object, ...) {
  structure(list(result = object), class = "summary.raddiv_differentiation")
}

#' @export
print.summary.raddiv_differentiation <- function(x, ...) {
  res <- x$result
  st <- res$settings
  note <- function(...) cat(paste0("  ", c(...), "\n"), sep = "")
  rule <- function(char) cat(strrep(char, 69), "\n", sep = "")
  cat("\n")
  rule("=")
  cat("  DIFFERENTIATION --", st$n_pops, "populations,", .big(st$n_records),
      "records on", .big(st$n_loci), "RAD loci\n")
  cat("  FST/FIS = Weir & Cockerham (1984); D = Jost (2008)\n")
  if (st$nboot > 0)
    cat("  95% CI from ", .big(st$nboot), " bootstrap replicates over RAD loci\n", sep = "")
  rule("=")
  cat("\n")
  .print_table(.round_table(res$global))
  note(.report_text$differentiation_se_columns)
  if (!is.null(st$d_records_global)) {
    note(paste0(.d_records_sentence(st$d_records_global, st$n_records, st$n_pops), "."),
         .report_text$differentiation_d_global_records)
    if (.d_records_few(st))
      note("NOTE: global D rests on fewer than half the records here; report pairwise D as well.")
  }
  skipped <- st$fst_records_skipped
  if (!is.null(skipped) && any(skipped > 0))
    note(sprintf("FST skipped records typed in fewer than 2 compared populations, or with one individual per typed population (%s):",
                 paste(sprintf("%s %s", names(skipped)[skipped > 0], .big(skipped[skipped > 0])),
                       collapse = ", ")),
         .report_text$differentiation_one_pop_records)
  cat("\nPairwise FST / (beta) / D\n")
  .print_table(.round_table(res$pairwise))
  if (!is.null(res$pairwise_beta))
    note("beta = Weir & Goudet (2017), via hierfstat::pairwise.betas() -- a point estimate only.")
  else
    note(if (isFALSE(st$beta)) "beta not computed (beta = FALSE)."
         else "beta not available (hierfstat not installed).")
  cat("\n")
  rule("-")
  if (st$is_haplotype) note(.report_text$differentiation_haplotype)
  rule("-")
  if (length(st$files)) cat("\nWrote ", paste(st$files, collapse = " and "), "\n", sep = "")
  cat("\n")
  invisible(x)
}
