###############################################################################
#
#  R/differentiation_stats.R -- how much do these populations DIFFER from
#  each other (as opposed to diversity_stats(), which asks how variable each
#  population is on its own)?
#
#  Two different numbers answer that question, and this file computes both:
#
#    FST   Weir & Cockerham (1984) -- the classic, most widely reported
#          measure. Ranges roughly 0 (no differentiation) to 1 (populations
#          share no alleles). Its ceiling shrinks as loci get more variable,
#          which can make FST look artificially small on a very diverse
#          haplotype-VCF marker set -- this is a known, much-discussed
#          property of FST itself, not a bug in this implementation.
#    D     Jost's D (Jost 2008, Molecular Ecology 17:4015) -- built to fix
#          exactly that ceiling problem. Recommended (by Jost's own paper,
#          and widely since) as the better choice for highly variable
#          markers such as RAD haplotypes.
#
#  Also reported when hierfstat is installed: Weir & Goudet's (2017) beta, a
#  newer FST-family estimator that does not need to designate one population
#  as the "reference" -- see hierfstat::pairwise.betas() -- offered
#  alongside FST/D, not in place of them.
#
#  ENGINE AND VALIDATION. FST/FIS here are computed with this package's OWN
#  implementation of the Weir & Cockerham (1984) variance-component formula
#  (functions .wc_components()/.fst_from_sums()/.fis_from_sums()
#  below), NOT by calling hierfstat, even when hierfstat is installed --
#  like everything else this package reports. That is a deliberate choice,
#  not an oversight: the bootstrap
#  and jackknife below need to re-sum these per-locus building blocks tens
#  of thousands of times, and hierfstat has no function that accepts a
#  resampled/weighted locus set. The code implements Weir & Cockerham's
#  (1984) published estimators of the a, b and c variance components (also
#  given in Weir 1996, Genetic Data Analysis II) -- see the comment on
#  `b_comp` in .wc_components() below for the one term that is easy to
#  mistranscribe. No hierfstat code is used or copied. hierfstat serves only
#  as an external numerical check: tests/testthat/test-differentiation-stats.R
#  compares this implementation with hierfstat::wc() and
#  hierfstat::pairwise.WCfst() on a stress-test dataset (4 populations, 25
#  loci with 2-4 alleles each, random missingness INCLUDING a population with
#  zero individuals genotyped at a locus), and they agree to machine
#  precision (~1e-16) for the global value and every pairwise value.
#  hierfstat is also used at run time, when installed, for Weir & Goudet's
#  beta (which this package does not reimplement) and -- only with
#  hierfstat_check = TRUE -- a printed cross-check of this run's own global
#  FST against hierfstat::wc()'s.
#
#  Jost's D has no hierfstat equivalent. `.jost_hsht()` below is an
#  ADAPTATION -- not an independent re-derivation from the Jost (2008)
#  paper -- of the `mmod` package's `HsHt()` function (Winter 2012;
#  functions HsHt() and D.per.locus() in that package's own R/HsHt.R and
#  R/D_Jost.R): the same five-quantity formula sequence (harmN / HpS /
#  Hs_est / HpT / Ht_est, in that order), reimplemented against this
#  package's `H$A1`/`H$A2` integer-genotype matrices in place of `mmod`'s
#  `genind`-object/allele-table representation. `mmod` is MIT-licensed
#  (same license as this package), which permits this; the required
#  attribution is given here, in `.jost_hsht()`'s own comment, and in
#  `inst/NOTICE`. This package's default multi-locus combination rule
#  (average Hs and Ht across loci FIRST, then compute one D from the
#  averaged values) matches `mmod`'s own default ("global.het" with
#  hsht_mean = "arithmetic") rather than the alternative of averaging many
#  separate per-locus D values -- a real, documented ambiguity in how
#  Jost's D is usually reported, so this choice is stated here explicitly
#  rather than left implicit.
#
#  PRECISION: this package's own RAD-tag-BLOCK bootstrap (resampling whole
#  RAD loci, i.e. locus_raw groups, with replacement -- the same mechanism
#  diversity_stats() uses and validates, NOT hierfstat::boot.ppfst(), which
#  resamples individual SNP rows independently and so ignores RAD-tag
#  linkage) for a 95% CI, plus a matching delete-one-block jackknife for a
#  standard error that needs no bootstrap replicates at all.
#
#  REFERENCES
#    Weir, B.S. & Cockerham, C.C. (1984) Estimating F-statistics for the
#      analysis of population structure. Evolution 38: 1358-1370.
#      https://doi.org/10.1111/j.1558-5646.1984.tb05657.x
#    Weir, B.S. & Goudet, J. (2017) A unified characterization of
#      population structure and relatedness. Genetics 206: 2085-2103.
#      https://doi.org/10.1534/genetics.116.198424
#    Jost, L. (2008) GST and its relatives do not measure differentiation.
#      Molecular Ecology 17: 4015-4026.
#      https://doi.org/10.1111/j.1365-294X.2008.03887.x
#    Winter, D.J. (2012) mmod: an R library for the calculation of
#      population differentiation statistics. Molecular Ecology Resources
#      12: 1158-1160. https://doi.org/10.1111/j.1755-0998.2012.03174.x
#      (MIT-licensed; see inst/NOTICE for the attribution .jost_hsht()
#      below owes it.)
#
###############################################################################

## ---------------------------------------------------------------------------
## Weir & Cockerham (1984) building blocks
## ---------------------------------------------------------------------------

## Not exported. For ONE specific set of populations (`pops`, a named list of
## sample-ID vectors -- either every population in the dataset, for a GLOBAL
## FST, or just two of them, for a PAIRWISE FST between that pair), this
## computes three numbers PER RECORD (PER LOCUS): `comp_a`, `comp_b`, `comp_c`.
## These are the three pieces of variance Weir & Cockerham's method splits
## every locus's genetic variation into:
##   comp_a  variance AMONG the populations listed in `pops`     ("a")
##   comp_b  variance among INDIVIDUALS WITHIN those populations ("b")
##   comp_c  variance WITHIN INDIVIDUALS, i.e. heterozygosity     ("c")
## FST is then just how big a share of the total (a+b+c) is the "among
## populations" share (a) -- see .fst_from_sums() below, which is
## where these three numbers actually turn into one FST value. Kept
## separate, per-record, so the SAME three numbers can be resummed however
## many times a bootstrap or jackknife replicate needs, without redoing this
## (somewhat involved) per-genotype counting each time.
##
## HOW TO READ THIS FUNCTION IF YOU DON'T PROGRAM IN R: the outer `for (j ...)`
## loop below walks one RAD locus (one row of H$A1/H$A2) at a time. Inside
## it, `cmat` counts how many copies of each allele were seen in each
## population (like locus_allele_stats(), but split out by population
## instead of pooled), and `hmat` counts how many INDIVIDUALS in each
## population were heterozygous and carrying that allele. Everything below
## `has_data` is arithmetic on those two small count tables.
.wc_components <- function(H, pops) {
  A1 <- H$A1; A2 <- H$A2; n_alleles <- H$n_alleles
  r <- length(pops); n_rec <- nrow(A1)
  comp_a <- comp_b <- comp_c <- numeric(n_rec)
  for (j in seq_len(n_rec)) {
    na_j <- n_alleles[j]
    cmat <- matrix(0, r, na_j)  # gene-copy counts: cmat[pop, allele]
    hmat <- matrix(0, r, na_j)  # heterozygote-CARRIER counts: hmat[pop, allele]
    n_i  <- integer(r)          # typed (non-missing) individuals per population
    for (i in seq_len(r)) {
      a <- A1[j, pops[[i]]]; b <- A2[j, pops[[i]]]
      ok <- !is.na(a)
      n_i[i] <- sum(ok)
      if (n_i[i] == 0L) next    # nobody in this population was typed here
      a <- a[ok]; b <- b[ok]
      cmat[i, ] <- tabulate(c(a, b), nbins = na_j)
      ## A heterozygote (a != b) carries ONE copy of allele `a` and ONE copy
      ## of allele `b` -- both alleles' heterozygote-carrier counts go up by
      ## one for that individual, which is why both tabulate() calls below
      ## are added together rather than averaged.
      het <- a != b
      if (any(het))
        hmat[i, ] <- hmat[i, ] + tabulate(a[het], nbins = na_j) + tabulate(b[het], nbins = na_j)
    }
    has_data <- n_i > 0L
    r_eff <- sum(has_data)      # how many of these populations have ANY data here
    if (r_eff < 1L || sum(n_i) < 2L) { comp_a[j] <- comp_b[j] <- comp_c[j] <- NA_real_; next }
    ni <- n_i[has_data]; cm <- cmat[has_data, , drop = FALSE]; hm <- hmat[has_data, , drop = FALSE]
    nt <- sum(ni); nbar <- nt / r_eff        # total, and average, typed individuals
    p_bar <- colSums(cm) / (2 * nt)          # pooled allele frequency, per allele
    h_bar <- colSums(hm) / nt                # pooled heterozygote-carrier frequency
    s2 <- if (r_eff > 1L) {
      ## Sample-size-WEIGHTED variance of each population's own allele
      ## frequency around the pooled frequency p_bar -- the "how different
      ## are these populations from each other" quantity FST is ultimately
      ## built from.
      p_i <- sweep(cm, 1, 2 * ni, "/")
      colSums(ni * sweep(p_i, 2, p_bar, "-")^2) / ((r_eff - 1) * nbar)
    } else 0
    ## The "among populations" component (a) needs at least two POPULATIONS
    ## with data at this record, and needs nbar > 1 (a single genotyped
    ## individual can't inform a variance-among-populations term). A record
    ## that fails this contributes exactly 0 here. hierfstat::wc() ends up
    ## treating such records the same way (its value there is not-a-number,
    ## which its na.rm = TRUE sums drop), so the two agree on them. `comp_b`/`comp_c` below
    ## are NOT gated the same way: a record where only one population out of
    ## several has any data still has a perfectly well-defined "within
    ## individuals" (heterozygosity) component, so it still counts there.
    a_comp <- if (r_eff > 1L && nbar > 1) {
      nc <- (nt - sum(ni^2) / nt) / (r_eff - 1)
      (nbar / nc) * (s2 - (1 / (nbar - 1)) *
        (p_bar * (1 - p_bar) - ((r_eff - 1) / r_eff) * s2 - h_bar / 4))
    } else 0
    ## NOTE ON THIS FORMULA: this term uses (n̄ - 1) [the plain average
    ## sample size], NOT (ñ - 1) [the same corrected sample size used for
    ## `a_comp` above] -- an easy pair to swap by mistake. This is Weir &
    ## Cockerham's (1984) b estimator as published; the numerical check
    ## against hierfstat::wc() in tests/testthat/test-differentiation-stats.R
    ## (agreement to ~1e-16) would catch the swap.
    b_comp <- if (nbar > 1) {
      (nbar / (nbar - 1)) * (p_bar * (1 - p_bar) - ((r_eff - 1) / r_eff) * s2 -
        ((2 * nbar - 1) / (4 * nbar)) * h_bar)
    } else 0
    c_comp <- h_bar / 2
    comp_a[j] <- sum(a_comp); comp_b[j] <- sum(b_comp); comp_c[j] <- sum(c_comp)
  }
  list(comp_a = comp_a, comp_b = comp_b, comp_c = comp_c)
}

## Not exported. FST, and the "metapopulation" FIS wc() reports, from a, b
## and c SUMMED over whichever loci an estimate uses (every locus for the
## point estimate, a resampled multiset for a bootstrap replicate, all but
## one locus for a jackknife replicate -- see R/resampling.R). A RATIO OF
## SUMS, never a mean of per-locus ratios -- the same principle this package
## applies to FIS in R/estimators.R's fis_ratio_of_sums(), for the same
## reason (a locus with a tiny denominator would otherwise be able to
## dominate the average):
##   FST = sum(a) / sum(a + b + c)  -- the share of variance AMONG populations
##   FIS = sum(b) / sum(b + c)      -- how much of the within-population part
##                                     is among individuals rather than within
## `n` is how many loci had all components defined. Vectorised: one element
## per estimate. No defined locus, or a total of ~0 (every locus monomorphic),
## gives NA; the 1e-12 guard (rather than == 0) absorbs the last-bit rounding
## a jackknife leaves when it subtracts one locus from the total.
.fst_from_sums <- function(n, a, b, c) {
  d <- a + b + c
  ifelse(n < 0.5 | abs(d) <= 1e-12, NA_real_, a / d)
}
.fis_from_sums <- function(n, b, c) {
  d <- b + c
  ifelse(n < 0.5 | abs(d) <= 1e-12, NA_real_, b / d)
}

## ---------------------------------------------------------------------------
## Jost's D (2008) building blocks
## ---------------------------------------------------------------------------

## Not exported. For one specific set of populations `pops`, computes Jost's
## bias-corrected within-population gene diversity (Hs_est) and total gene
## diversity (Ht_est) PER RECORD -- the two ingredients D is built from (see
## .jost_d_from_sums() below). Deliberately a DIFFERENT weighting
## scheme than .wc_components() above: Jost's own formula averages
## populations EQUALLY (a plain mean across populations), not weighted by
## how many individuals each one has -- that is not a simplification on
## this package's part, it is what the cited formula does.
##
## ATTRIBUTION: the five-step formula below (harmN -> HpS -> Hs_est -> HpT
## -> Ht_est, in that order) is ADAPTED from the `mmod` R package's
## `HsHt()` function (Winter 2012; R/HsHt.R in that package's own source),
## not independently re-derived from the Jost (2008) paper -- read
## honestly, this is a translation of mmod's algorithm from its
## `genind`-object/allele-table representation into this package's
## `H$A1`/`H$A2` integer-genotype-matrix representation, not a fresh
## derivation. `mmod` is MIT-licensed (Copyright (c) 2014 David Winter),
## the same license as this package, which permits this adaptation;
## see `inst/NOTICE` for the required attribution text.
.jost_hsht <- function(H, pops) {
  A1 <- H$A1; A2 <- H$A2; n_alleles <- H$n_alleles
  r <- length(pops); n_rec <- nrow(A1)
  Hs <- Ht <- rep(NA_real_, n_rec)
  for (j in seq_len(n_rec)) {
    na_j <- n_alleles[j]
    freqs <- matrix(NA_real_, r, na_j)   # freqs[pop, allele]
    n_i <- integer(r)
    for (i in seq_len(r)) {
      a <- A1[j, pops[[i]]]; b <- A2[j, pops[[i]]]
      ok <- !is.na(a)
      n_i[i] <- sum(ok)
      if (n_i[i] == 0L) next
      freqs[i, ] <- tabulate(c(a[ok], b[ok]), nbins = na_j) / (2 * n_i[i])
    }
    has_data <- n_i > 0L
    r_eff <- sum(has_data)
    if (r_eff < 2L) next   # Jost's D needs at least 2 populations WITH DATA here
    fq <- freqs[has_data, , drop = FALSE]; ni <- n_i[has_data]
    ## Harmonic mean of population sizes -- Jost's own choice of "typical"
    ## sample size for the small-sample bias correction below (a harmonic
    ## mean is pulled down hard by any one small population, which is the
    ## intended, conservative behaviour here).
    harmN <- 1 / mean(1 / ni)
    HpS <- mean(1 - rowSums(fq^2))         # naive within-pop diversity, averaged over pops
    Hs_est <- (2 * harmN / (2 * harmN - 1)) * HpS
    HpT <- 1 - sum(colMeans(fq)^2)         # naive total diversity, from pop-averaged frequencies
    Ht_est <- HpT + Hs_est / (2 * harmN * r_eff)
    Hs[j] <- Hs_est; Ht[j] <- Ht_est
  }
  list(Hs = Hs, Ht = Ht)
}

## Not exported. Combines Hs_est/Ht_est SUMMED over whichever records an
## estimate uses (`n` of them with both defined) into one D value, following
## Jost's (2008) own "average-Hs-and-Ht-first" convention (matching the
## `mmod` package's default -- see the file header comment above): `n_pop` is
## the FIXED number of populations in this particular comparison (2 for one
## pairwise D, or the full population count for the global D) -- NOT how many
## of them happened to have data at any one record, which is why it is passed
## in rather than read off the per-record data. Vectorised over estimates,
## like .fst_from_sums().
.jost_d_from_sums <- function(n, hs, ht, n_pop) {
  hs_bar <- ifelse(n < 0.5, NA_real_, hs / n)
  ht_bar <- ht / n
  ifelse(is.na(hs_bar) | hs_bar >= 1, NA_real_,
         (ht_bar - hs_bar) / (1 - hs_bar) * (n_pop / (n_pop - 1)))
}

## ---------------------------------------------------------------------------
## The exported function
## ---------------------------------------------------------------------------

#' Population differentiation: FST, Weir & Goudet's beta, and Jost's D
#'
#' Computes how much these populations differ from each other -- as opposed
#' to [diversity_stats()], which asks how variable each population is on its
#' own. Reports Weir & Cockerham's (1984) FST and FIS (global, i.e. across
#' every population at once, and pairwise, i.e. one value per pair of
#' populations), Jost's (2008) D (global and pairwise, generally the better
#' choice on highly variable markers such as RAD haplotypes -- see
#' `@details`), and, when the `hierfstat` package is installed, Weir &
#' Goudet's (2017) pairwise beta as a third, reference-population-free
#' alternative.
#'
#' @details
#' **Why report both FST and D.** FST has a ceiling that shrinks as a
#' locus's own diversity grows, so two datasets with genuinely identical
#' population structure can show very different FST values just because one
#' used more variable markers (this is a well known, much-discussed property
#' of FST itself, described at length in Jost (2008), the paper this
#' function's D estimator comes from). RAD-haplotype loci are exactly the
#' kind of highly variable marker where this matters, so on a haplotype VCF,
#' prefer D (or beta) over FST for comparing datasets -- the same caution
#' [diversity_stats()] already gives for allelic richness and gene diversity
#' on haplotype data.
#'
#' **Bootstrap and jackknife, not hierfstat's.** Confidence intervals here
#' always come from this package's own block bootstrap over RAD loci (the
#' same mechanism [diversity_stats()] validates and uses), never from
#' `hierfstat::boot.ppfst()`, which resamples individual SNP rows
#' independently and so ignores the linkage between SNPs on the same RAD
#' tag. The internal FST/FIS formula is Weir & Cockerham's (1984) published
#' estimator; the package tests check it numerically against
#' `hierfstat::wc()` (see the comment at the top of
#' `R/differentiation_stats.R`).
#'
#' **References.** Weir, B.S. & Cockerham, C.C. (1984) Estimating
#' F-statistics for the analysis of population structure. *Evolution*
#' 38:1358-1370. <https://doi.org/10.1111/j.1558-5646.1984.tb05657.x> --
#' Weir, B.S. & Goudet, J. (2017) A unified characterization of population
#' structure and relatedness. *Genetics* 206:2085-2103.
#' <https://doi.org/10.1534/genetics.116.198424> -- Jost, L. (2008) GST and
#' its relatives do not measure differentiation. *Molecular Ecology*
#' 17:4015-4026. <https://doi.org/10.1111/j.1365-294X.2008.03887.x> --
#' Winter, D.J. (2012) mmod: an R library for the calculation of population
#' differentiation statistics. *Molecular Ecology Resources* 12:1158-1160.
#' <https://doi.org/10.1111/j.1755-0998.2012.03174.x> (this package's
#' `.jost_hsht()` is an attributed adaptation of mmod's `HsHt()`; see
#' `inst/NOTICE`).
#'
#' @inheritParams diversity_stats
#' @param vcf_file Path to a Stacks VCF, or an already-parsed `H` list (see
#'   [diversity_stats()] for the same option there).
#' @param popmap_f Path to a two-column, no-header popmap TSV (`sample_id
#'   <TAB> population`).
#' @param nboot Bootstrap replicates over RAD loci. Default `10000L`; `0` =
#'   none (only the jackknife standard error is then available).
#' @param hierfstat_check If `TRUE` and `hierfstat` is installed, also run
#'   `hierfstat::wc()` and report how far its global FST is from this
#'   package's. Default `FALSE` (the package tests already make this
#'   comparison, and `wc()` was most of the run time on large datasets).
#' @return An object of class `raddiv_differentiation` (printing it shows the
#'   report): a list with
#'   \describe{
#'     \item{global}{One-row data frame: `FST`, `FST_se`, `FST_lo`, `FST_hi`,
#'       `FIS`, `FIS_se`, `FIS_lo`, `FIS_hi`, `D`, `D_se`, `D_lo`, `D_hi`.}
#'     \item{pairwise}{One row per pair of populations: FST, beta and D, with
#'       their standard errors and bootstrap CIs (also written to
#'       `differentiation_pairwise.<stem>.tsv` when `outdir` is given, as the
#'       global table is to `differentiation_global.<stem>.tsv`).}
#'     \item{pairwise_fst, pairwise_beta, pairwise_D}{The pairwise point
#'       estimates as symmetric matrices (`NA` on the diagonal), one
#'       row/column per population. `pairwise_beta` is `NULL` (with an
#'       explanatory message) when hierfstat is not installed.}
#'   }
#' @examples
#' # A tiny made-up dataset: 2 populations, 4 loci, biallelic.
#' samp <- c("a1","a2","a3","a4","b1","b2","b3","b4")
#' H <- list(
#'   A1 = rbind(locus_1 = c(1,1,1,1, 2,2,2,1), locus_2 = c(1,2,1,2, 1,1,2,2),
#'              locus_3 = c(1,1,2,1, 1,1,1,1), locus_4 = c(2,1,1,2, 1,2,1,1)),
#'   A2 = rbind(locus_1 = c(1,1,2,1, 2,2,1,2), locus_2 = c(2,2,1,1, 1,2,1,2),
#'              locus_3 = c(1,2,2,1, 1,1,2,1), locus_4 = c(2,2,1,2, 1,2,2,1)),
#'   locus = paste0("locus_", 1:4), locus_raw = paste0("locus_", 1:4),
#'   n_alleles = rep(2L, 4), alleles = replicate(4, c("A", "C"), simplify = FALSE),
#'   samples = samp
#' )
#' colnames(H$A1) <- colnames(H$A2) <- samp
#' popmap <- tempfile()
#' writeLines(c("a1\tpopA","a2\tpopA","a3\tpopA","a4\tpopA",
#'              "b1\tpopB","b2\tpopB","b3\tpopB","b4\tpopB"), popmap)
#' res <- differentiation_stats(H, popmap, nboot = 100, verbose = FALSE)
#' res$global
#' @export
differentiation_stats <- function(vcf_file, popmap_f, nboot = 10000L,
                                   outdir = NULL, seed = 2024, verbose = TRUE,
                                   stem = NULL, hierfstat_check = FALSE) {

  nboot <- as.integer(nboot)
  if (is.na(nboot) || nboot < 0) stop("nboot must be a non-negative integer.")
  .check_run_inputs(vcf_file, popmap_f, stem, outdir)

  ## Same RNG-preservation convention as diversity_stats()/het_between_pops():
  ## restore the caller's own random-number state on exit, so calling this
  ## interactively never changes what random numbers the caller's own code
  ## gets afterwards.
  restore_rng <- .save_rng_state()
  on.exit(restore_rng(), add = TRUE)
  set.seed(seed)

  have_hf <- .hierfstat_available()

  if (is.character(vcf_file)) message("Reading ", vcf_file, " ...")
  H <- .resolve_H(vcf_file, verbose = verbose)
  pops <- read_popmap(popmap_f, H$samples, verbose = verbose)
  r <- length(pops)
  if (r < 2) stop("Need at least 2 populations.")
  nmax <- lengths(pops)
  tiny <- names(pops)[nmax < 2]
  if (length(tiny))
    stop("Population(s) with fewer than 2 individuals: ", paste(tiny, collapse = ", "),
         "\n  FST/D need at least 2 individuals per population to define a ",
         "within-population allele frequency.")

  n_rec <- nrow(H$A1)
  rad <- H$locus_raw
  uloc <- unique(rad); li <- match(rad, uloc); nL <- length(uloc)
  message(sprintf("  %s records on %s RAD loci (%.2f per locus)",
                  format(n_rec, big.mark = ","), format(nL, big.mark = ","), n_rec / nL))

  pair_names <- utils::combn(names(pops), 2, simplify = FALSE)
  message(sprintf("Computing FST/D for %d population(s), %d pair(s) ...", r, length(pair_names)))

  ## Per-record variance components / Hs-Ht, computed ONCE for the global
  ## (all-populations) comparison and once per pair -- these do not change
  ## across bootstrap/jackknife replicates, only which RECORDS get summed
  ## does, so this is the expensive part done exactly once.
  wc_global <- .wc_components(H, pops)
  d_global  <- .jost_hsht(H, pops)
  wc_pair <- lapply(pair_names, function(pn) .wc_components(H, pops[pn]))
  d_pair  <- lapply(pair_names, function(pn) .jost_hsht(H, pops[pn]))

  ## Per-record pieces, summed once per RAD locus into `S` (see
  ## R/resampling.R): for each comparison -- the global one, then each pair --
  ## FST needs (n, a, b, c) over the records where a, b and c are all defined,
  ## and D needs (n, Hs, Ht) over the records where both are defined; the
  ## global FIS keeps its own (n, b, c) mask. An undefined record adds 0 to
  ## every sum and to the count. Columns: 1-4 global FST, 5-7 global FIS,
  ## 8-10 global D, then 4 per pair (FST), then 3 per pair (D).
  pair_lab <- vapply(pair_names, paste, character(1), collapse = "__")
  n_pair <- length(pair_names)
  wc_pieces <- function(w) {
    ok <- is.finite(w$comp_a) & is.finite(w$comp_b) & is.finite(w$comp_c)
    cbind(ok, ifelse(ok, w$comp_a, 0), ifelse(ok, w$comp_b, 0), ifelse(ok, w$comp_c, 0))
  }
  d_pieces <- function(d) {
    ok <- is.finite(d$Hs) & is.finite(d$Ht)
    cbind(ok, ifelse(ok, d$Hs, 0), ifelse(ok, d$Ht, 0))
  }
  fis_ok <- is.finite(wc_global$comp_b) & is.finite(wc_global$comp_c)
  rec <- cbind(wc_pieces(wc_global),
               fis_ok, ifelse(fis_ok, wc_global$comp_b, 0), ifelse(fis_ok, wc_global$comp_c, 0),
               d_pieces(d_global),
               do.call(cbind, lapply(wc_pair, wc_pieces)),
               do.call(cbind, lapply(d_pair, d_pieces)))
  storage.mode(rec) <- "double"
  S <- rowsum(rec, li, reorder = TRUE)          # nL x pieces; row b = RAD locus b
  stat_names <- c("FST", "FIS", "D", paste0("pFST_", pair_lab), paste0("pD_", pair_lab))
  stat_mat <- function(tot) {
    tot <- matrix(tot, ncol = ncol(S))
    col <- function(j) tot[, j]
    per_pair <- function(first, f)                 # one column per pair
      matrix(vapply(first, f, numeric(nrow(tot))), nrow = nrow(tot))
    out <- cbind(
      .fst_from_sums(col(1), col(2), col(3), col(4)),
      .fis_from_sums(col(5), col(6), col(7)),
      .jost_d_from_sums(col(8), col(9), col(10), r),
      per_pair(10L + 4L * (seq_len(n_pair) - 1L), function(j)
        .fst_from_sums(col(j + 1L), col(j + 2L), col(j + 3L), col(j + 4L))),
      per_pair(10L + 4L * n_pair + 3L * (seq_len(n_pair) - 1L), function(j)
        .jost_d_from_sums(col(j + 1L), col(j + 2L), col(j + 3L), 2)))
    colnames(out) <- stat_names
    out
  }
  point <- stat_mat(colSums(S))[1, ]

  ## Delete-one-BLOCK jackknife over RAD loci -- deterministic, needs no
  ## nboot, the same as diversity_stats() (R/resampling.R).
  se_jk <- .jack_block_sums(S, stat_mat)

  if (nboot > 0) {
    message(sprintf("Bootstrapping %s replicates over %s RAD loci ...",
                    format(nboot, big.mark = ","), format(nL, big.mark = ",")))
    bt <- t(.boot_block_sums(S, nboot, stat_mat))
    ci <- t(apply(bt, 1, stats::quantile, c(0.025, 0.975), na.rm = TRUE))
  } else {
    ci <- matrix(NA_real_, length(point), 2, dimnames = list(names(point), NULL))
  }

  ## With hierfstat installed: Weir & Goudet's beta (a reported statistic this
  ## package does not reimplement), and -- only with hierfstat_check = TRUE --
  ## a cross-check of this run's own FST against hierfstat::wc() (message the
  ## difference; warn only if it is surprisingly large). This package's own
  ## formula is what gets reported either way.
  if (have_hf) {
    dat <- .to_hierfstat_df(H, pops)
    if (isTRUE(hierfstat_check)) {
      hw <- try(hierfstat::wc(dat), silent = TRUE)
      if (!inherits(hw, "try-error")) {
        d_fst <- abs(hw$FST - point[["FST"]])
        message(sprintf("  cross-check vs hierfstat::wc(): FST diff %.2e", d_fst))
        if (d_fst > 1e-6) warning("This run's internal FST disagrees with hierfstat::wc().")
      }
    }
    beta_mat <- try(hierfstat::pairwise.betas(dat), silent = TRUE)
    if (inherits(beta_mat, "try-error")) {
      message("  hierfstat::pairwise.betas() failed on this dataset; pairwise_beta will be NULL.")
      beta_mat <- NULL
    } else {
      beta_mat <- as.matrix(beta_mat)
      dimnames(beta_mat) <- list(names(pops), names(pops))
    }
  } else {
    message("Install hierfstat for Weir & Goudet's beta: install.packages(\"hierfstat\")")
    beta_mat <- NULL
  }

  ## ---------------------------------------------------------------------------
  ## Result tables. Nothing is printed here: printing the returned object
  ## (print.raddiv_differentiation(), below) gives the report.
  ## ---------------------------------------------------------------------------
  gl <- data.frame(
    FST = round(point[["FST"]], 4), FST_se = round(se_jk[["FST"]], 4),
    FST_lo = round(ci["FST", 1], 4), FST_hi = round(ci["FST", 2], 4),
    FIS = round(point[["FIS"]], 4), FIS_se = round(se_jk[["FIS"]], 4),
    FIS_lo = round(ci["FIS", 1], 4), FIS_hi = round(ci["FIS", 2], 4),
    D   = round(point[["D"]], 4),   D_se   = round(se_jk[["D"]], 4),
    D_lo = round(ci["D", 1], 4),    D_hi   = round(ci["D", 2], 4))

  pfst <- matrix(NA_real_, r, r, dimnames = list(names(pops), names(pops)))
  pD   <- matrix(NA_real_, r, r, dimnames = list(names(pops), names(pops)))
  pair_rows <- vector("list", length(pair_names))
  for (k in seq_along(pair_names)) {
    p1 <- pair_names[[k]][1]; p2 <- pair_names[[k]][2]
    fst_k <- point[[paste0("pFST_", pair_lab[k])]]; d_k <- point[[paste0("pD_", pair_lab[k])]]
    pfst[p1, p2] <- pfst[p2, p1] <- fst_k
    pD[p1, p2]   <- pD[p2, p1]   <- d_k
    beta_k <- if (!is.null(beta_mat)) beta_mat[p1, p2] else NA_real_
    pair_rows[[k]] <- data.frame(
      pop1 = p1, pop2 = p2,
      FST = round(fst_k, 4), FST_se = round(se_jk[[paste0("pFST_", pair_lab[k])]], 4),
      FST_lo = round(ci[paste0("pFST_", pair_lab[k]), 1], 4),
      FST_hi = round(ci[paste0("pFST_", pair_lab[k]), 2], 4),
      beta = round(beta_k, 4),
      D = round(d_k, 4), D_se = round(se_jk[[paste0("pD_", pair_lab[k])]], 4),
      D_lo = round(ci[paste0("pD_", pair_lab[k]), 1], 4),
      D_hi = round(ci[paste0("pD_", pair_lab[k]), 2], 4))
  }
  pair_tab <- do.call(rbind, pair_rows)

  ## Files, only when asked for (`outdir`).
  written <- character(0)
  if (!is.null(outdir)) {
    dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
    stem <- .derive_stem(vcf_file, stem)
    written <- file.path(outdir, sprintf(c("differentiation_global.%s.tsv",
                                           "differentiation_pairwise.%s.tsv"), stem))
    utils::write.table(gl, written[1], sep = "\t", quote = FALSE, row.names = FALSE)
    utils::write.table(pair_tab, written[2], sep = "\t", quote = FALSE, row.names = FALSE)
  }

  structure(
    list(global = gl, pairwise = pair_tab, pairwise_fst = pfst,
         pairwise_beta = beta_mat, pairwise_D = pD),
    class = "raddiv_differentiation",
    report = list(n_pop = r, n_rec = n_rec, n_loci = nL, nboot = nboot,
                  is_haplotype = .is_haplotype_H(H), files = written))
}

#' @rdname differentiation_stats
#' @param x A `raddiv_differentiation` object, as returned by
#'   `differentiation_stats()`.
#' @param ... Ignored.
#' @export
print.raddiv_differentiation <- function(x, ...) {
  rp <- attr(x, "report")
  old <- options(width = max(200L, getOption("width")))
  on.exit(options(old), add = TRUE)
  cat("\n=====================================================================\n")
  cat("  DIFFERENTIATION --", rp$n_pop, "populations,", format(rp$n_rec, big.mark = ","),
      "records on", format(rp$n_loci, big.mark = ","), "RAD loci\n")
  cat("  FST/FIS = Weir & Cockerham (1984); D = Jost (2008)\n")
  if (rp$nboot > 0) cat("  95% CI from ", format(rp$nboot, big.mark = ","),
                        " bootstrap replicates over RAD loci\n", sep = "")
  cat("=====================================================================\n\n")
  print(x$global, row.names = FALSE)
  cat("  _se columns: delete-one-block jackknife over RAD loci (same convention\n")
  cat("  as diversity_stats()); _lo/_hi: bootstrap 95% CI over the same RAD loci.\n")
  cat("\nPairwise FST / (beta) / D\n")
  print(x$pairwise, row.names = FALSE)
  if (!is.null(x$pairwise_beta))
    cat("  beta = Weir & Goudet (2017), via hierfstat::pairwise.betas() -- no bootstrap CI (hierfstat-only, point estimate).\n")
  else
    cat("  beta not available (hierfstat not installed).\n")
  cat("\n---------------------------------------------------------------------\n")
  if (rp$is_haplotype) {
    cat("  This is a HAPLOTYPE VCF: prefer D (or beta) over FST here -- FST's\n")
    cat("  ceiling shrinks as marker diversity grows, and haplotype loci are\n")
    cat("  exactly the highly-variable-marker case that distorts (see\n")
    cat("  ?differentiation_stats and Jost 2008).\n")
  }
  cat("---------------------------------------------------------------------\n")
  if (length(rp$files)) cat("\nWrote ", paste(rp$files, collapse = " and "), "\n", sep = "")
  cat("\n")
  invisible(x)
}
