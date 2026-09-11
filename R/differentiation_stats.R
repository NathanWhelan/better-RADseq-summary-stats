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
#  (functions .wc_components()/.fst_from_components()/.fis_from_components()
#  below), NOT by calling hierfstat, even when hierfstat is installed --
#  unlike diversity_stats(), which prefers hierfstat's own functions when
#  available. That is a deliberate choice, not an oversight: the bootstrap
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
#  hierfstat is also used at run time, when installed, for two things this
#  package does not reimplement: Weir & Goudet's beta, and a printed
#  cross-check of this run's own global FST against hierfstat::wc()'s.
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
## computes three numbers PER RECORD (PER LOCUS): `siga`, `sigb`, `sigw`.
## These are the three pieces of variance Weir & Cockerham's method splits
## every locus's genetic variation into:
##   siga  variance AMONG the populations listed in `pops`     ("a")
##   sigb  variance among INDIVIDUALS WITHIN those populations ("b")
##   sigw  variance WITHIN INDIVIDUALS, i.e. heterozygosity     ("c")
## FST is then just how big a share of the total (a+b+c) is the "among
## populations" share (a) -- see .fst_from_components() below, which is
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
  siga <- sigb <- sigw <- numeric(n_rec)
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
    if (r_eff < 1L || sum(n_i) < 2L) { siga[j] <- sigb[j] <- sigw[j] <- NA_real_; next }
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
    ## that fails this contributes exactly 0 here -- confirmed against
    ## hierfstat::wc()'s own source, which reaches the same "contributes
    ## nothing" outcome via a different route (a not-a-number value that its
    ## na.rm = TRUE summation silently treats as zero). `sigb`/`sigw` below
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
    siga[j] <- sum(a_comp); sigb[j] <- sum(b_comp); sigw[j] <- sum(c_comp)
  }
  list(siga = siga, sigb = sigb, sigw = sigw)
}

## Not exported. Turns a's/b's/c's (summed over whichever records a
## statistic needs -- every record for the point estimate, a resampled
## multiset for one bootstrap replicate, all-but-one-block for a jackknife
## replicate) into one FST value: a RATIO OF SUMS, never a mean of per-locus
## ratios -- the same principle this package already applies to FIS in
## R/estimators.R's fis_ratio_of_sums(), for the same reason (a locus with a
## tiny denominator would otherwise be able to dominate the average).
.fst_from_components <- function(siga, sigb, sigw) {
  ok <- is.finite(siga) & is.finite(sigb) & is.finite(sigw)
  if (!any(ok)) return(NA_real_)
  denom <- sum(siga[ok]) + sum(sigb[ok]) + sum(sigw[ok])
  if (!isTRUE(denom != 0)) return(NA_real_)
  sum(siga[ok]) / denom
}
## Not exported. Same idea, for the "metapopulation" FIS wc() reports: how
## much of the (within-population + within-individual) variance is within
## POPULATIONS rather than within individuals.
.fis_from_components <- function(sigb, sigw) {
  ok <- is.finite(sigb) & is.finite(sigw)
  if (!any(ok)) return(NA_real_)
  denom <- sum(sigb[ok]) + sum(sigw[ok])
  if (!isTRUE(denom != 0)) return(NA_real_)
  sum(sigb[ok]) / denom
}

## ---------------------------------------------------------------------------
## Jost's D (2008) building blocks
## ---------------------------------------------------------------------------

## Not exported. For one specific set of populations `pops`, computes Jost's
## bias-corrected within-population gene diversity (Hs_est) and total gene
## diversity (Ht_est) PER RECORD -- the two ingredients D is built from (see
## .jost_d_from_components() below). Deliberately a DIFFERENT weighting
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

## Not exported. Combines per-record Hs_est/Ht_est (over whichever records a
## statistic needs) into one D value, following Jost's (2008) own
## "average-Hs-and-Ht-first" convention (matching the `mmod` package's
## default -- see the file header comment above): `n_pop` is the FIXED
## number of populations in this particular comparison (2 for one pairwise
## D, or the full population count for the global D) -- NOT how many of them
## happened to have data at any one record, which is why it is passed in
## rather than read off the per-record data.
.jost_d_from_components <- function(Hs, Ht, n_pop) {
  ok <- is.finite(Hs) & is.finite(Ht)
  if (!any(ok)) return(NA_real_)
  hs_bar <- mean(Hs[ok]); ht_bar <- mean(Ht[ok])
  if (!isTRUE(hs_bar < 1)) return(NA_real_)
  (ht_bar - hs_bar) / (1 - hs_bar) * (n_pop / (n_pop - 1))
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
#' @return Invisibly, a list:
#'   \describe{
#'     \item{global}{One-row data frame: `FST`, `FST_se`, `FST_lo`, `FST_hi`,
#'       `FIS`, `FIS_se`, `FIS_lo`, `FIS_hi`, `D`, `D_se`, `D_lo`, `D_hi`.}
#'     \item{pairwise_fst, pairwise_beta, pairwise_D}{Symmetric matrices
#'       (`NA` on the diagonal), one row/column per population, of POINT
#'       ESTIMATES only. `pairwise_beta` is `NULL` (with an explanatory
#'       message) when hierfstat is not installed. The full detail
#'       (standard errors and bootstrap CIs per pair) is in the
#'       `differentiation_pairwise.*.tsv` file this function writes, in
#'       long format -- one row per pair.}
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
#' res <- differentiation_stats(H, popmap, nboot = 100, stem = "example",
#'                               outdir = tempdir(), verbose = FALSE)
#' res$global
#' @export
differentiation_stats <- function(vcf_file, popmap_f, nboot = 10000L,
                                   outdir = ".", seed = 2024, verbose = TRUE,
                                   stem = NULL) {

  nboot <- as.integer(nboot)
  if (is.na(nboot) || nboot < 0) stop("nboot must be a non-negative integer.")
  if (is.character(vcf_file) && !file.exists(vcf_file))
    stop("VCF file not found: ", vcf_file, "\n  Check the path and try again.")
  if (!is.character(vcf_file) && is.null(stem))
    stop("vcf_file is an already-parsed list rather than a file path, so its ",
         "filename can't be used to name the output files. Pass stem ",
         "explicitly, e.g. stem = \"haps\" or stem = \"snps\".")
  if (!file.exists(popmap_f))
    stop("Popmap file not found: ", popmap_f, "\n  Check the path and try again.")
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

  ## Same RNG-preservation convention as diversity_stats()/het_between_pops():
  ## restore the caller's own random-number state on exit, so calling this
  ## interactively never changes what random numbers the caller's own code
  ## gets afterwards.
  restore_rng <- .save_rng_state()
  on.exit(restore_rng(), add = TRUE)
  set.seed(seed)

  have_hf <- requireNamespace("hierfstat", quietly = TRUE)

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

  ## One statistic name per quantity this run tracks, so a single stat_from()
  ## can compute the whole named vector at once, exactly like
  ## diversity_stats()'s own stat_from() -- reused by the point estimate, the
  ## bootstrap, and the jackknife alike.
  pair_lab <- vapply(pair_names, paste, character(1), collapse = "__")
  stat_from <- function(sel) {
    g_fst <- .fst_from_components(wc_global$siga[sel], wc_global$sigb[sel], wc_global$sigw[sel])
    g_fis <- .fis_from_components(wc_global$sigb[sel], wc_global$sigw[sel])
    g_d   <- .jost_d_from_components(d_global$Hs[sel], d_global$Ht[sel], r)
    p_fst <- vapply(wc_pair, function(w) .fst_from_components(w$siga[sel], w$sigb[sel], w$sigw[sel]), numeric(1))
    p_d   <- vapply(d_pair, function(d) .jost_d_from_components(d$Hs[sel], d$Ht[sel], 2), numeric(1))
    stats::setNames(c(g_fst, g_fis, g_d, p_fst, p_d),
             c("FST", "FIS", "D", paste0("pFST_", pair_lab), paste0("pD_", pair_lab)))
  }
  point <- stat_from(seq_len(n_rec))

  ## Delete-one-BLOCK jackknife over RAD loci -- deterministic, needs no
  ## nboot, exactly mirroring diversity_stats()'s own jack_se().
  jk <- vapply(seq_len(nL), function(b) stat_from(which(li != b)), numeric(length(point)))
  theta_bar <- rowMeans(jk)
  se_jk <- sqrt(((nL - 1) / nL) * rowSums((jk - theta_bar)^2))
  names(se_jk) <- names(point)

  if (nboot > 0) {
    message(sprintf("Bootstrapping %s replicates over %s RAD loci ...",
                    format(nboot, big.mark = ","), format(nL, big.mark = ",")))
    rows_of <- split(seq_len(n_rec), li)
    resample_loci <- function() unlist(rows_of[sample.int(nL, nL, replace = TRUE)], use.names = FALSE)
    bt <- replicate(nboot, stat_from(resample_loci()))
    ci <- t(apply(bt, 1, stats::quantile, c(0.025, 0.975), na.rm = TRUE))
  } else {
    ci <- matrix(NA_real_, length(point), 2, dimnames = list(names(point), NULL))
  }

  ## Cross-check against hierfstat's own wc(), when installed -- exactly the
  ## same style of cross-check diversity_stats() already does for allelic
  ## richness (message the max difference; warn only if it is surprisingly
  ## large). This package's own formula is what actually gets reported (see
  ## the file header for why), hierfstat's copy is only a second opinion.
  if (have_hf) {
    ids <- unlist(pops, use.names = FALSE)
    popvec <- rep(seq_len(r), lengths(pops))
    Gm <- matrix(NA_integer_, length(ids), n_rec)
    for (j in seq_len(n_rec)) {
      a <- H$A1[j, ids]; b <- H$A2[j, ids]
      Gm[, j] <- pmin(a, b) * 1000L + pmax(a, b)
    }
    dat <- data.frame(pop = popvec, Gm); names(dat)[-1] <- paste0("L", seq_len(n_rec))
    hw <- try(hierfstat::wc(dat), silent = TRUE)
    if (!inherits(hw, "try-error")) {
      d_fst <- abs(hw$FST - point[["FST"]])
      message(sprintf("  cross-check vs hierfstat::wc(): FST diff %.2e", d_fst))
      if (d_fst > 1e-6) warning("This run's internal FST disagrees with hierfstat::wc().")
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
  ## Report
  ## ---------------------------------------------------------------------------
  cat("\n=====================================================================\n")
  cat("  DIFFERENTIATION --", r, "populations,", format(n_rec, big.mark = ","),
      "records on", format(nL, big.mark = ","), "RAD loci\n")
  cat("  FST/FIS = Weir & Cockerham (1984); D = Jost (2008)\n")
  if (nboot > 0) cat("  95% CI from ", format(nboot, big.mark = ","), " bootstrap replicates over RAD loci\n", sep = "")
  cat("=====================================================================\n\n")
  gl <- data.frame(
    FST = round(point[["FST"]], 4), FST_se = round(se_jk[["FST"]], 4),
    FST_lo = round(ci["FST", 1], 4), FST_hi = round(ci["FST", 2], 4),
    FIS = round(point[["FIS"]], 4), FIS_se = round(se_jk[["FIS"]], 4),
    FIS_lo = round(ci["FIS", 1], 4), FIS_hi = round(ci["FIS", 2], 4),
    D   = round(point[["D"]], 4),   D_se   = round(se_jk[["D"]], 4),
    D_lo = round(ci["D", 1], 4),    D_hi   = round(ci["D", 2], 4))
  print(gl, row.names = FALSE)
  cat("  _se columns: delete-one-block jackknife over RAD loci (same convention\n")
  cat("  as diversity_stats()); _lo/_hi: bootstrap 95% CI over the same RAD loci.\n")

  pfst <- matrix(NA_real_, r, r, dimnames = list(names(pops), names(pops)))
  pD    <- matrix(NA_real_, r, r, dimnames = list(names(pops), names(pops)))
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
  cat("\nPairwise FST / (beta) / D\n")
  print(pair_tab, row.names = FALSE)
  if (!is.null(beta_mat))
    cat("  beta = Weir & Goudet (2017), via hierfstat::pairwise.betas() -- no bootstrap CI (hierfstat-only, point estimate).\n")
  else
    cat("  beta not available (hierfstat not installed) -- see above.\n")

  cat("\n---------------------------------------------------------------------\n")
  is_hapvcf <- max(H$n_alleles) > 2L || any(nchar(unlist(H$alleles, use.names = FALSE)) > 1L)
  if (is_hapvcf) {
    cat("  This is a HAPLOTYPE VCF: prefer D (or beta) over FST here -- FST's\n")
    cat("  ceiling shrinks as marker diversity grows, and haplotype loci are\n")
    cat("  exactly the highly-variable-marker case that distorts (see @details,\n")
    cat("  ?differentiation_stats, and Jost 2008).\n")
  }
  cat("---------------------------------------------------------------------\n")

  stem <- .derive_stem(vcf_file, stem)
  f1 <- sprintf("differentiation_global.%s.tsv", stem)
  f2 <- sprintf("differentiation_pairwise.%s.tsv", stem)
  utils::write.table(gl, file.path(outdir, f1), sep = "\t", quote = FALSE, row.names = FALSE)
  utils::write.table(pair_tab, file.path(outdir, f2), sep = "\t", quote = FALSE, row.names = FALSE)
  cat(sprintf("\nWrote %s and %s\n\n", file.path(outdir, f1), file.path(outdir, f2)))

  invisible(list(global = gl, pairwise_fst = pfst, pairwise_beta = beta_mat, pairwise_D = pD))
}
