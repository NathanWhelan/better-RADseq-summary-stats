###############################################################################
#
#  R/filter_loci.R -- general-purpose locus and genotype filters.
#
#  Every function here takes `H`, the list returned by read_haps_vcf() (see
#  R/vcf_io.R), and returns a filtered version of it: same shape, same
#  elements (A1, A2, locus, locus_raw, alleles, n_alleles, samples, fields),
#  just fewer rows (loci) and/or some genotype cells set to missing. That
#  means these filters CHAIN: you can run several in a row, in any order,
#  e.g.
#
#    H <- read_haps_vcf("populations.snps.vcf")
#    H <- filter_call_rate(H, min_call = 0.8)
#    H <- filter_maf(H, min_maf = 0.05)
#    H <- filter_low_conf_alt(H, min_alt_reads = 2)$H
#
#  and the result can be handed straight to diversity_stats()/
#  het_between_pops() in place of a file path (see .resolve_H() in
#  R/vcf_io.R and the `stem` argument of diversity_stats()).
#
#  These are DIFFERENT from, and do not touch, the filtering logic already
#  built into diversity_stats() (the min_n/complete_case missing-data rule)
#  or het_between_pops() (the pooled min_call rule). Those two are each
#  tuned to one specific statistical need and are documented at length in
#  their own files -- the functions below are general-purpose, reusable
#  building blocks for cleaning up a dataset BEFORE analysis, export, or
#  both.
#
#  A NOTE FOR READERS NEW TO R: throughout this file, `NA` is R's way of
#  writing "no data here" (a missing value), and a "logical" vector is just
#  a vector of TRUE/FALSE (yes/no) values. `H$A1` and `H$A2` are matrices --
#  tables of numbers -- with one row per RAD locus and one column per
#  sample; `H$A1[j, i]` is individual i's FIRST allele at locus j and
#  `H$A2[j, i]` its SECOND allele, both stored as small integers (1 = the
#  REF allele, 2 = the first ALT allele, and so on), or `NA` if that
#  individual wasn't genotyped there.
#
###############################################################################

## Not exported. Every filter below ends by calling this to drop some rows
## (loci) from H and keep the rest. Written once here so every filter drops
## rows the same, correct way instead of repeating (and possibly
## mis-copying) the same seven lines of subsetting code.
##
## `keep` can be:
##   - a TRUE/FALSE (logical) vector, one value per locus, saying whether to
##     keep that locus, or
##   - a vector of row numbers to keep (e.g. c(1, 3, 4)).
## Either way, every piece of H that has "one entry per locus" is subset
## together, so they all stay lined up with each other afterwards. `samples`
## isn't touched here, since subsetting loci never changes which individuals
## are in the dataset.
.subset_H <- function(H, keep) {
  n_rec <- nrow(H$A1)
  if (is.logical(keep)) {
    if (length(keep) != n_rec || anyNA(keep))
      stop(".subset_H(): `keep` must be a TRUE/FALSE vector with exactly ",
           n_rec, " entries (one per locus in H) and no missing values.")
  } else {
    keep <- as.integer(keep)
    if (anyNA(keep) || any(keep < 1L | keep > n_rec))
      stop(".subset_H(): `keep` row numbers must all be between 1 and ", n_rec, ".")
  }
  ## `A1[keep, , drop = FALSE]` keeps only the rows named in `keep`, for every
  ## column (sample); `drop = FALSE` stops R from silently turning the result
  ## into a plain vector if only one locus (or one sample) is left.
  H$A1        <- H$A1[keep, , drop = FALSE]
  H$A2        <- H$A2[keep, , drop = FALSE]
  H$locus     <- H$locus[keep]
  H$locus_raw <- H$locus_raw[keep]
  H$alleles   <- H$alleles[keep]
  H$n_alleles <- H$n_alleles[keep]
  ## `fields` (the raw VCF text) isn't always present -- e.g. a small H built
  ## by hand for testing might skip it -- so only subset it if it's there.
  if (!is.null(H$fields)) H$fields <- H$fields[keep, , drop = FALSE]
  H
}

#' Per-locus allele-frequency statistics
#'
#' Computes, for every locus in `H`, how many allele copies were actually
#' genotyped, which allele was the most common ("major") one, the minor
#' allele frequency, and the minor allele count. This is the shared
#' groundwork behind [filter_maf()] and [filter_mac()] -- call it yourself
#' first if you want to look at the numbers before deciding on a threshold
#' (e.g. `hist(locus_allele_stats(H)$maf)`), or pass its result into
#' `filter_maf()`/`filter_mac()` via their `stats` argument to avoid
#' recomputing it twice when you're going to apply both filters.
#'
#' How this works, step by step: for every locus, every called allele copy
#' (from both `H$A1` and `H$A2`, ignoring anything missing) is counted up by
#' allele number. The allele with the highest count is the "major" allele;
#' everything else, added together, is the "minor" share. This is the usual
#' definition of minor allele frequency (MAF) when a locus has exactly two
#' possible alleles, and it generalizes sensibly to loci with more than two
#' (haplotype VCFs can have several alleles per locus): MAF there means "the
#' combined frequency of every allele except the single most common one".
#'
#' @param H A list as returned by [read_haps_vcf()] (or by another filter in
#'   this package, since they all return the same shape).
#' @return A data frame with one row per locus, in the same order as `H$A1`,
#'   and columns:
#'   \describe{
#'     \item{locus, locus_raw, n_alleles}{Copied straight from `H`, for
#'       convenience.}
#'     \item{n_observed_alleles}{How many DIFFERENT allele numbers were
#'       actually seen in the genotypes at this locus. Can be smaller than
#'       `n_alleles` if the VCF lists an ALT allele that nobody happened to
#'       carry.}
#'     \item{n_called}{Total number of allele copies genotyped at this
#'       locus (each diploid individual contributes up to 2).}
#'     \item{major_af}{Frequency of the single most common allele.}
#'     \item{maf}{Minor allele frequency, `1 - major_af`.}
#'     \item{mac}{Minor allele count: how many allele copies were NOT the
#'       major allele.}
#'   }
#'   A locus with zero genotyped individuals gets `NA` (not `0`) for
#'   `major_af`/`maf`/`mac`, since "zero data" and "definitely monomorphic"
#'   are different claims and must not be confused.
#' @examples
#' # A tiny made-up dataset: 2 loci, 3 samples, alleles coded 1 (REF) and
#' # 2 (ALT). read_haps_vcf() builds this same shape from a real VCF file.
#' # Each row of A1/A2 is one locus's alleles across all 3 samples.
#' H <- list(
#'   A1 = rbind(locus_1 = c(1, 1, 2), locus_2 = c(1, 2, NA)),
#'   A2 = rbind(locus_1 = c(1, 1, 2), locus_2 = c(2, 2, NA)),
#'   locus = c("locus_1", "locus_2"), locus_raw = c("locus_1", "locus_2"),
#'   n_alleles = c(2L, 2L), samples = c("ind1", "ind2", "ind3")
#' )
#' locus_allele_stats(H)  # locus_2 is missing in ind3, so n_called = 4 there
#' @export
locus_allele_stats <- function(H) {
  A1 <- H$A1; A2 <- H$A2
  n_rec  <- nrow(A1)
  n_samp <- ncol(A1)
  mx     <- max(H$n_alleles)  # the most alleles any single locus has

  ## Line up every genotyped allele copy (from A1 and from A2) next to the
  ## locus (row) it belongs to, so we can count them all in one pass instead
  ## of looping over loci one at a time.
  ##   c(A1, A2)          -- "unrolls" both matrices into one long list of
  ##                         numbers, column by column.
  ##   rep(rep(seq_len(n_rec), n_samp), 2)
  ##                      -- builds a matching list saying which locus (row
  ##                         number) each of those numbers came from.
  codes <- c(A1, A2)
  rows  <- rep(rep(seq_len(n_rec), n_samp), 2L)
  ok    <- !is.na(codes)  # drop the "no data here" entries before counting

  ## tabulate() is a fast base-R counting tool: give it a list of "bin"
  ## numbers and it tells you how many times each bin number showed up. By
  ## giving each (locus, allele) combination its own bin number, one call to
  ## tabulate() counts every allele at every locus at once.
  bin <- (codes[ok] - 1L) * n_rec + rows[ok]
  tab <- tabulate(bin, nbins = n_rec * mx)
  dim(tab) <- c(n_rec, mx)  # tab[j, a] = how many times allele `a` was seen at locus j

  n_called           <- rowSums(tab)
  n_observed_alleles <- rowSums(tab > 0L)

  ## The "major" allele count at each locus is just the largest count in
  ## that row of `tab`. pmax() compares two vectors position-by-position and
  ## keeps whichever value is bigger at each position; running it once per
  ## allele column finds the row-wise maximum without needing apply().
  major_count <- tab[, 1L]
  if (mx > 1L) for (a in 2:mx) major_count <- pmax(major_count, tab[, a])

  major_af <- major_count / n_called
  mac      <- n_called - major_count
  maf      <- mac / n_called

  ## A locus with n_called == 0 (nobody genotyped) divides zero by zero
  ## above, giving NaN ("not a number") -- replace that with a proper NA so
  ## it reads as "unknown", and make sure `mac` becomes NA too (not the 0 it
  ## would otherwise compute to) so nothing downstream mistakes "no data"
  ## for "definitely monomorphic".
  no_data <- n_called == 0L
  major_af[no_data] <- NA_real_
  maf[no_data]      <- NA_real_
  mac[no_data]       <- NA_real_

  data.frame(
    locus = H$locus, locus_raw = H$locus_raw, n_alleles = H$n_alleles,
    n_observed_alleles = n_observed_alleles, n_called = n_called,
    major_af = major_af, maf = maf, mac = as.integer(round(mac)),
    row.names = NULL
  )
}

## Shared by filter_maf() and filter_mac(): make sure a user-supplied `stats`
## data frame actually corresponds to this exact H (same loci, same order),
## and compute it fresh if none was supplied.
.check_or_make_stats <- function(H, stats) {
  if (is.null(stats)) return(locus_allele_stats(H))
  if (nrow(stats) != nrow(H$A1) || !identical(stats$locus, H$locus))
    stop("`stats` does not line up with `H` (different number of loci, or ",
         "the loci are in a different order). Pass the `stats` returned by ",
         "locus_allele_stats(H) computed on this exact H, or leave `stats` ",
         "at its default (NULL) to have it computed automatically.")
  stats
}

#' Filter loci by minor allele frequency
#'
#' Keeps only loci whose minor allele frequency (MAF) is at least `min_maf`.
#' Rare variants are the ones most likely to be sequencing/genotyping
#' errors rather than real biological variation, so a MAF filter is one of
#' the most common first QC steps for SNP data.
#'
#' @param H A list as returned by [read_haps_vcf()].
#' @param min_maf Minimum minor allele frequency a locus must have to be
#'   kept (a number between 0 and 1, e.g. `0.05` for 5%).
#' @param stats Optionally, the data frame already returned by
#'   [locus_allele_stats()] for this exact `H` -- pass it in to avoid
#'   recomputing it if you're also calling [filter_mac()] on the same
#'   dataset. Left as `NULL` (the default), it's computed automatically.
#' @param verbose Print how many loci were kept. Default `TRUE`.
#' @return `H`, with only the loci that passed the filter kept (a locus with
#'   no genotyped individuals at all has an undefined MAF and is always
#'   dropped, since it can't be shown to clear the threshold).
#' @examples
#' # Each row of A1/A2 is one locus's alleles across 3 samples; locus_1 is
#' # homozygous REF (allele 1) in everyone, locus_2 has 2 REF and 4 ALT
#' # copies (MAF = 2/6 = 0.33):
#' H <- list(
#'   A1 = rbind(locus_1 = c(1, 1, 1), locus_2 = c(1, 2, 2)),
#'   A2 = rbind(locus_1 = c(1, 1, 1), locus_2 = c(1, 2, 2)),
#'   locus = c("locus_1", "locus_2"), locus_raw = c("locus_1", "locus_2"),
#'   alleles = list(c("A", "C"), c("A", "C")), n_alleles = c(2L, 2L),
#'   samples = c("ind1", "ind2", "ind3")
#' )
#' # only locus_2 clears a 20% MAF threshold:
#' filter_maf(H, min_maf = 0.2, verbose = FALSE)$locus
#' @export
filter_maf <- function(H, min_maf, stats = NULL, verbose = TRUE) {
  stats <- .check_or_make_stats(H, stats)
  n_rec <- nrow(stats)
  ## `!is.na(...)` excludes loci with no data at all; `- 1e-9` guards against
  ## a locus landing exactly on the threshold but reading as just barely
  ## below it purely due to floating-point rounding.
  keep <- !is.na(stats$maf) & stats$maf >= min_maf - 1e-9
  if (verbose)
    message(sprintf(
      "MAF filter (keep loci with minor allele frequency >= %.4g): %s of %s loci kept (%.1f%%)",
      min_maf, format(sum(keep), big.mark = ","), format(n_rec, big.mark = ","), 100 * mean(keep)))
  if (sum(keep) < 0.05 * n_rec && verbose)
    message("  WARNING: fewer than 5% of loci passed this filter. Consider a lower min_maf.")
  if (!sum(keep))
    stop("No locus has a minor allele frequency >= ", min_maf, ". Lower min_maf, ",
         "or run locus_allele_stats(H) yourself to see the actual distribution of values.")
  .subset_H(H, keep)
}

#' Filter loci by minor allele count
#'
#' Keeps only loci with at least `min_mac` copies of the minor allele
#' (counting across every genotyped individual). This is the count-based
#' cousin of [filter_maf()] -- useful because a fixed count threshold (e.g.
#' "at least 3 copies") behaves more predictably than a frequency threshold
#' when sample sizes are small or uneven, where a single individual can
#' swing the frequency a lot.
#'
#' @param H A list as returned by [read_haps_vcf()].
#' @param min_mac Minimum minor allele count a locus must have to be kept
#'   (a whole number, e.g. `3`).
#' @param stats Optionally, a precomputed [locus_allele_stats()] result for
#'   this `H` (see [filter_maf()] for why you might pass this in).
#' @param verbose Print how many loci were kept. Default `TRUE`.
#' @return `H`, with only the loci that passed the filter kept.
#' @examples
#' # locus_2 has 2 minor (REF) allele copies out of 6; locus_1 has none:
#' H <- list(
#'   A1 = rbind(locus_1 = c(1, 1, 1), locus_2 = c(1, 2, 2)),
#'   A2 = rbind(locus_1 = c(1, 1, 1), locus_2 = c(1, 2, 2)),
#'   locus = c("locus_1", "locus_2"), locus_raw = c("locus_1", "locus_2"),
#'   alleles = list(c("A", "C"), c("A", "C")), n_alleles = c(2L, 2L),
#'   samples = c("ind1", "ind2", "ind3")
#' )
#' filter_mac(H, min_mac = 2, verbose = FALSE)$locus  # only locus_2 has >= 2 minor copies
#' @export
filter_mac <- function(H, min_mac, stats = NULL, verbose = TRUE) {
  stats <- .check_or_make_stats(H, stats)
  n_rec <- nrow(stats)
  keep <- !is.na(stats$mac) & stats$mac >= min_mac
  if (verbose)
    message(sprintf(
      "MAC filter (keep loci with minor allele count >= %d): %s of %s loci kept (%.1f%%)",
      as.integer(min_mac), format(sum(keep), big.mark = ","), format(n_rec, big.mark = ","), 100 * mean(keep)))
  if (sum(keep) < 0.05 * n_rec && verbose)
    message("  WARNING: fewer than 5% of loci passed this filter. Consider a lower min_mac.")
  if (!sum(keep))
    stop("No locus has a minor allele count >= ", min_mac, ". Lower min_mac, ",
         "or run locus_allele_stats(H) yourself to see the actual distribution of values.")
  .subset_H(H, keep)
}

#' Filter loci by genotyping (call) rate
#'
#' Keeps only loci that were genotyped often enough. Without `pops`, "often
#' enough" means across the whole dataset pooled together; with `pops`
#' (population assignments from [read_popmap()]), it's judged separately in
#' each population, which matters because a locus can look well-covered
#' overall while actually being poorly covered in just one population --
#' see `rule` below.
#'
#' This is a standalone, general-purpose version of a call-rate check --
#' it does not change or replace the call-rate logic already built into
#' [het_between_pops()] (which is tuned specifically for that function's own
#' diagnostics) or the `min_n`/`complete_case` logic in [diversity_stats()].
#'
#' @param H A list as returned by [read_haps_vcf()].
#' @param min_call Minimum fraction of individuals that must be genotyped at
#'   a locus for it to be kept (a number between 0 and 1, e.g. `0.8` for
#'   80%).
#' @param pops Optionally, a named list of sample-ID vectors (the output of
#'   [read_popmap()]) to judge the call rate separately per population
#'   instead of pooling everyone together. Default `NULL` (pooled).
#' @param rule Only used when `pops` is given. `"all"` (the default) keeps a
#'   locus only if EVERY population individually clears `min_call` there --
#'   the right choice if you're going to compare populations to each other
#'   afterwards, since a pooled or "any population" rule can let one
#'   well-covered population effectively borrow coverage from a poorly-
#'   covered one. `"any"` keeps a locus if AT LEAST ONE population clears
#'   it -- appropriate if each population's own results will be used
#'   independently rather than compared directly.
#' @param verbose Print how many loci were kept (and, with `pops`, a
#'   per-population breakdown). Default `TRUE`.
#' @return `H`, with only the loci that passed the filter kept.
#' @examples
#' H <- list(
#'   A1 = rbind(locus_1 = c(1, 1, 1, NA), locus_2 = c(1, NA, 1, NA)),
#'   A2 = rbind(locus_1 = c(1, 1, 1, NA), locus_2 = c(2, NA, 1, NA)),
#'   locus = c("locus_1", "locus_2"), locus_raw = c("locus_1", "locus_2"),
#'   n_alleles = c(2L, 2L), samples = c("a1", "a2", "a3", "a4")
#' )
#' # locus_1 is genotyped in 3 of 4 samples (75%); locus_2 only in 2 of 4 (50%).
#' filter_call_rate(H, min_call = 0.75, verbose = FALSE)$locus
#' @export
filter_call_rate <- function(H, min_call, pops = NULL, rule = "all", verbose = TRUE) {
  if (!(length(rule) == 1L && rule %in% c("all", "any")))
    stop("rule must be exactly \"all\" or \"any\" (got: ", paste(rule, collapse = ", "), ").")
  n_rec <- nrow(H$A1)

  if (is.null(pops)) {
    ## rowMeans(!is.na(H$A1)) is, for every locus, the fraction of samples
    ## that have a real (non-missing) genotype there.
    cr   <- rowMeans(!is.na(H$A1))
    keep <- cr >= min_call - 1e-9
    if (verbose)
      message(sprintf(
        "Call-rate filter (pooled across all %d samples, threshold %.0f%%): %s of %s loci kept (%.1f%%)",
        ncol(H$A1), 100 * min_call, format(sum(keep), big.mark = ","),
        format(n_rec, big.mark = ","), 100 * mean(keep)))
  } else {
    r <- length(pops)
    cr_pop <- sweep(.typed_by_pop(H, pops), 2L, lengths(pops), "/")
    ok_pop <- cr_pop >= min_call - 1e-9
    ## rowSums(ok_pop) counts how many populations clear the threshold at
    ## each locus; "all" needs that count to equal every population, "any"
    ## just needs it to be at least one.
    keep <- if (rule == "all") rowSums(ok_pop) == r else rowSums(ok_pop) > 0L
    if (verbose) {
      message(sprintf(
        "Call-rate filter (per population, threshold %.0f%%, rule = \"%s\"): %s of %s loci kept (%.1f%%)",
        100 * min_call, rule, format(sum(keep), big.mark = ","),
        format(n_rec, big.mark = ","), 100 * mean(keep)))
      for (p in names(pops))
        message(sprintf("    %-22s %s of %s loci clear the threshold",
                        p, format(sum(ok_pop[, p]), big.mark = ","), format(n_rec, big.mark = ",")))
    }
  }
  if (!sum(keep))
    stop("No locus passed the call-rate filter (min_call = ", min_call, "). Lower min_call, ",
         "or check for a widespread genotyping problem in this dataset.")
  .subset_H(H, keep)
}

#' Filter out loci with excessive heterozygosity
#'
#' Drops loci where an unusually large share of individuals are called
#' heterozygous. Real heterozygosity has a biological ceiling; a locus far
#' above it is the classic signature of an undetected paralog (two
#' similar-looking genomic regions being misread as one locus) or a
#' collapsed repeat, both of which produce fake "heterozygous" calls that
#' are really two different genes/copies being confused for two alleles of
#' the same one.
#'
#' @param H A list as returned by [read_haps_vcf()].
#' @param max_ho Maximum allowed observed heterozygosity at a locus (a
#'   number between 0 and 1, e.g. `0.5`). A locus above this is dropped.
#' @param verbose Print how many loci were kept. Default `TRUE`.
#' @return `H`, with only the loci that passed the filter kept. A locus with
#'   no genotyped individuals at all is kept (not dropped) by this filter,
#'   since "no data" is not the same claim as "too much heterozygosity" --
#'   use [filter_call_rate()] if you also want to remove those.
#' @examples
#' H <- list(
#'   A1 = rbind(locus_1 = c(1, 1, 1, 1), locus_2 = c(1, 2, 1, 2)),
#'   A2 = rbind(locus_1 = c(1, 1, 1, 1), locus_2 = c(2, 1, 2, 1)),
#'   locus = c("locus_1", "locus_2"), locus_raw = c("locus_1", "locus_2"),
#'   n_alleles = c(2L, 2L), samples = c("a1", "a2", "a3", "a4")
#' )
#' # locus_2 is heterozygous in all 4 samples (Ho = 1) -- dropped at max_ho = 0.6:
#' filter_max_het(H, max_ho = 0.6, verbose = FALSE)$locus
#' @export
filter_max_het <- function(H, max_ho, verbose = TRUE) {
  ## For each locus, the fraction of genotyped individuals where the two
  ## alleles differ (A1 != A2) is the observed heterozygosity, Ho.
  ## na.rm = TRUE skips missing genotypes; a locus with NO genotyped
  ## individuals at all then comes out as NaN ("not a number", from 0/0),
  ## which !is.finite() below recognizes as "no data" rather than "high Ho".
  ho   <- rowMeans(H$A1 != H$A2, na.rm = TRUE)
  keep <- !is.finite(ho) | ho <= max_ho + 1e-9
  n_rec <- length(keep)
  if (verbose)
    message(sprintf(
      "Excess-heterozygosity filter (drop loci with Ho > %.3f): %s of %s loci kept (%.1f%%)",
      max_ho, format(sum(keep), big.mark = ","), format(n_rec, big.mark = ","), 100 * mean(keep)))
  if (!sum(keep))
    stop("Every locus in this dataset exceeds max_ho = ", max_ho, ". Raise max_ho, ",
         "or double-check that A1/A2 were parsed correctly.")
  .subset_H(H, keep)
}

#' Thin a SNP dataset to one record per RAD locus
#'
#' Some analyses assume every locus is inherited independently of every
#' other one. Several SNPs sharing the same RAD tag do NOT meet that
#' assumption -- they're physically linked, so treating them as independent
#' overstates how much information the dataset really contains. This
#' function keeps exactly one record per RAD tag (`H$locus_raw` group),
#' dropping the rest, so that every remaining locus is a genuinely
#' independent RAD tag.
#'
#' @param H A list as returned by [read_haps_vcf()].
#' @param method `"first"` (the default) keeps the first record seen for
#'   each RAD tag, giving the same result every time. `"random"` picks one
#'   record per RAD tag at random.
#' @param seed Only used when `method = "random"`. An integer to make the
#'   random choice reproducible. Left `NULL` (the default), the choice
#'   depends on R's current random-number state, same as calling
#'   `sample()` yourself without first calling `set.seed()`. Either way,
#'   your own R session's random-number state is left exactly as it was
#'   before this function was called (see `?set.seed` if you're not
#'   familiar with this) -- calling this function will never change what
#'   random numbers you get afterwards in your own code.
#' @param verbose Print how many records were kept. Default `TRUE`.
#' @return `H`, thinned to one record per RAD locus, with the kept records
#'   left in their original order.
#' @examples
#' # Two SNPs on locus_1's RAD tag, one SNP alone on locus_2's:
#' H <- list(
#'   A1 = matrix(1L, nrow = 3, ncol = 2), A2 = matrix(1L, nrow = 3, ncol = 2),
#'   locus = c("locus_1_a", "locus_1_b", "locus_2"),
#'   locus_raw = c("locus_1", "locus_1", "locus_2"),
#'   n_alleles = c(2L, 2L, 2L), samples = c("a1", "a2")
#' )
#' filter_thin_one_snp(H, verbose = FALSE)$locus  # keeps "locus_1_a" and "locus_2"
#' @export
filter_thin_one_snp <- function(H, method = "first", seed = NULL, verbose = TRUE) {
  if (!(length(method) == 1L && method %in% c("first", "random")))
    stop("method must be exactly \"first\" or \"random\" (got: ", paste(method, collapse = ", "), ").")
  n_rec <- nrow(H$A1)
  n_loc <- length(unique(H$locus_raw))
  if (n_loc == n_rec) {
    if (verbose)
      message("Every record already belongs to its own RAD locus (locus_raw has no ",
              "repeats) -- there is nothing to thin.")
    return(H)
  }

  ## split() sorts the row numbers 1:n_rec into one group per distinct
  ## locus_raw value, in the order those values first appear.
  groups <- split(seq_len(n_rec), factor(H$locus_raw, levels = unique(H$locus_raw)))

  if (method == "first") {
    keep <- vapply(groups, function(i) i[1L], integer(1))
  } else {
    restore_rng <- .save_rng_state()
    on.exit(restore_rng(), add = TRUE)
    if (!is.null(seed)) set.seed(seed)
    ## sample.int(length(i), 1) draws one random POSITION out of `length(i)`
    ## and returns i AT that position. This is deliberately not
    ## `sample(i, 1)`: R's everyday sample() function has a surprising
    ## special case where, if you give it a single number x (not a vector),
    ## it treats that as shorthand for "pick randomly from 1:x" instead of
    ## "there's only one option, x, so return it". Many RAD loci have only
    ## ONE record on this dataset (a locus_raw group of size 1), so getting
    ## this right is not just a style choice -- sample.int() avoids that trap.
    keep <- vapply(groups, function(i) i[sample.int(length(i), 1L)], integer(1))
  }
  keep <- sort(keep)  # put the kept records back into their original file order

  if (verbose)
    message(sprintf(
      "One-SNP-per-locus thinning (method = \"%s\"): %s records -> %s RAD loci (%.1f%% of records kept)",
      method, format(n_rec, big.mark = ","), format(length(keep), big.mark = ","), 100 * length(keep) / n_rec))
  .subset_H(H, keep)
}

###############################################################################
#
#  Porting filter_low_conf_alt.py: flagging genotype calls whose ALT allele
#  is backed by very few sequencing reads.
#
#  A genotype call like "0/1" (heterozygous) or "1/1" (homozygous ALT) means
#  Stacks decided this individual carries at least one copy of a non-
#  reference allele. But if only 1 or 2 reads out of many actually showed
#  that ALT allele, the call could easily be a sequencing error rather than
#  a real allele -- especially for a rare variant. The functions below
#  identify exactly those low-confidence ALT calls (using the AD field --
#  "allele depth", i.e. how many reads supported each allele -- from the
#  original VCF) so they can be masked out (treated as missing) or used to
#  drop the whole locus if too many of its individuals look unreliable.
#
###############################################################################

## Not exported. Shared by filter_low_conf_alt() and low_conf_alt_sensitivity()
## so the (somewhat fiddly) work of pulling AD/DP back out of the raw VCF
## text is only written once.
##
## Returns one row per genotype call that has at least one ALT allele and is
## fully called (no missing alleles) -- i.e. every candidate that COULD be
## flagged, whether or not it actually gets flagged. Columns:
##   record        row number in H (which locus)
##   sample        column number in H (which individual)
##   GT             the genotype text, e.g. "0/1"
##   AD_ref, AD_alt reads supporting the REF allele, and reads supporting
##                  whichever ALT allele(s) this genotype carries
##   DP             total read depth Stacks used for this call
##   alt_fraction   AD_alt / DP
##   usable         TRUE if AD (and so AD_alt) could actually be read from
##                  the file. Some VCF rows legitimately drop trailing
##                  FORMAT fields when they're not needed (this is allowed by
##                  the VCF file format spec) -- such a call still clearly
##                  has an ALT allele, it just can't be checked for
##                  low-confidence support, so it is never flagged.
.parse_alt_ad <- function(H) {
  A1 <- H$A1; A2 <- H$A2
  n_rec <- nrow(A1)

  ## A cell "has an ALT allele" if both alleles are called (not missing) and
  ## at least one of them is bigger than 1 (allele 1 is always REF).
  alt_cell <- !is.na(A1) & !is.na(A2) & (A1 > 1L | A2 > 1L)
  idx <- which(alt_cell)  # a single "linear" index per flagged cell, R's usual
                          # way of numbering matrix cells column-by-column
  empty <- data.frame(record = integer(), sample = integer(), GT = character(),
                       AD_ref = integer(), AD_alt = integer(), DP = integer(),
                       alt_fraction = double(), usable = logical())
  if (!length(idx)) return(empty)

  ## Turn each linear index back into a (row, column) = (record, sample) pair.
  rec_i  <- ((idx - 1L) %% n_rec) + 1L
  samp_i <- ((idx - 1L) %/% n_rec) + 1L

  ## FORMAT (e.g. "GT:AD:DP") tells us WHERE in each sample's colon-separated
  ## genotype text the AD and DP values sit -- and that position is not
  ## always the same in every VCF. Stacks almost always uses one constant
  ## FORMAT for the whole file, so instead of re-parsing it for every single
  ## flagged cell, we parse each *distinct* FORMAT string once and then look
  ## up the right answer for each record from that small table.
  fmt_all   <- H$fields[, "FORMAT"]
  fmt_u     <- unique(fmt_all)
  fmt_split <- strsplit(fmt_u, ":", fixed = TRUE)
  find_pos  <- function(tag) vapply(fmt_split, function(f) {
    p <- which(f == tag); if (length(p)) p[1] else NA_integer_
  }, integer(1))
  ad_pos_u <- find_pos("AD")
  dp_pos_u <- find_pos("DP")
  fmt_map  <- match(fmt_all, fmt_u)     # for each record, which row of *_u applies
  ad_pos   <- ad_pos_u[fmt_map[rec_i]]  # for each FLAGGED CELL, its record's AD position
  dp_pos   <- dp_pos_u[fmt_map[rec_i]]

  ## Pull out the raw genotype text ("0/1:20:5,15:20:99" or similar) for
  ## each flagged cell, by record and by SAMPLE NAME (not a fixed column
  ## number), so this keeps working even if H$fields' column order ever
  ## differs from the usual layout.
  samp_col <- match(H$samples, colnames(H$fields))
  raw   <- H$fields[cbind(rec_i, samp_col[samp_i])]
  parts <- strsplit(raw, ":", fixed = TRUE)

  n <- length(idx)
  gt_str <- vapply(parts, `[`, character(1), 1L)  # GT is always the first subfield

  usable <- !is.na(ad_pos) & lengths(parts) >= ad_pos
  ad_str <- rep(NA_character_, n)
  ad_str[usable] <- vapply(which(usable), function(k) parts[[k]][ad_pos[k]], character(1))
  ## "." or a blank means the AD value itself is missing even though there
  ## was room for it -- also not usable.
  usable <- usable & !is.na(ad_str) & ad_str != "." & nzchar(ad_str)

  AD_ref <- rep(NA_integer_, n)
  AD_alt <- rep(NA_integer_, n)
  DP     <- rep(NA_integer_, n)
  a1 <- A1[idx]; a2 <- A2[idx]

  for (k in which(usable)) {
    ## AD is a comma-separated list, one count per possible allele in REF,ALT
    ## order (so ad[1] = REF reads, ad[2] = reads for the first ALT allele,
    ## and so on). A "." inside that list means "0 reads", per the VCF spec.
    ad_txt <- strsplit(ad_str[k], ",", fixed = TRUE)[[1]]
    ad <- suppressWarnings(as.integer(ifelse(ad_txt == ".", "0", ad_txt)))

    ## Add up AD for every DISTINCT ALT allele actually present in this
    ## genotype. "Distinct" matters for a homozygous ALT call like 1/1:
    ## both copies are the SAME allele, so its read support must only be
    ## counted once, not doubled.
    alts_here <- unique(c(a1[k], a2[k]))
    alts_here <- alts_here[alts_here > 1L]
    AD_alt[k] <- sum(vapply(alts_here, function(a) if (a <= length(ad)) ad[a] else 0L, integer(1)))
    AD_ref[k] <- if (length(ad) >= 1L) ad[1L] else 0L

    ## DP (total read depth) is usually given directly; if it's missing or
    ## ".", fall back to adding up the AD values instead.
    dp_val <- NA_integer_
    if (!is.na(dp_pos[k]) && length(parts[[k]]) >= dp_pos[k]) {
      s <- parts[[k]][dp_pos[k]]
      if (!is.na(s) && s != "." && nzchar(s)) dp_val <- suppressWarnings(as.integer(s))
    }
    DP[k] <- if (is.na(dp_val)) sum(ad) else dp_val
  }

  data.frame(
    record = rec_i, sample = samp_i, GT = gt_str,
    AD_ref = AD_ref, AD_alt = AD_alt, DP = DP,
    alt_fraction = ifelse(!is.na(DP) & DP > 0L, AD_alt / DP, NA_real_),
    usable = usable
  )
}

#' Flag and remove low-confidence ALT genotype calls
#'
#' Ports the standalone `filter_low_conf_alt.py` script's logic natively
#' into R. A genotype call carrying an ALT allele (e.g. `0/1`, `1/1`) is
#' flagged if the number of sequencing reads supporting that ALT allele is
#' `<= min_alt_reads` -- i.e. Stacks called an alternate allele on the
#' strength of very little read evidence, which is exactly the profile of a
#' sequencing-error false positive. You then choose what to do with flagged
#' calls: `mode = "mask"` blanks out just those individual genotype calls
#' (leaving the rest of that locus, and every other individual, untouched);
#' `mode = "drop"` instead removes whole loci where too large a share of
#' their ALT-containing calls were flagged.
#'
#' @param H A list as returned by [read_haps_vcf()].
#' @param min_alt_reads A genotype call with this many or fewer reads
#'   supporting its ALT allele(s) is flagged. Default `2`. See
#'   [low_conf_alt_sensitivity()] for a table to help pick this value.
#' @param mode `"mask"` (the default) sets only the flagged genotype calls
#'   to missing. `"drop"` removes whole loci instead (see `drop_frac`).
#'   Unlike the original Python script, this function can only apply ONE of
#'   these per call (it returns a single `H`, not two separate files) -- if
#'   you want both views, call this function twice on the SAME starting
#'   `H`, once with each `mode`, and keep the two results separately.
#'   Calling `mode = "drop"` on an `H` you already ran through
#'   `mode = "mask"` will not reproduce that: the flagged calls are already
#'   gone by then, so there is nothing left for a second pass to drop.
#' @param drop_frac Only used when `mode = "drop"`. A locus is dropped if
#'   the fraction of its ALT-containing calls that were flagged is greater
#'   than this (default `0`, meaning "drop a locus if even one of its calls
#'   is flagged"). A locus with no ALT-containing calls at all is never
#'   dropped by this rule.
#' @param calls Optionally, a precomputed result from the internal
#'   `.parse_alt_ad()` helper, to avoid re-parsing the same data if you're
#'   also calling [low_conf_alt_sensitivity()] on this `H`. Left `NULL`
#'   (the default), it's computed automatically.
#' @param verbose Print flagging counts. Default `TRUE`.
#' @return A list with three elements:
#'   \describe{
#'     \item{H}{The filtered/masked `H`.}
#'     \item{flagged_calls}{One row per flagged genotype call: `locus,
#'       sample, GT, AD_ref, AD_alt, DP, alt_fraction`.}
#'     \item{locus_summary}{One row per ORIGINAL locus, in original order:
#'       `locus, n_alt_calls, n_flagged, flagged_fraction, kept` -- `kept`
#'       records whether that locus survived (always `TRUE` under
#'       `mode = "mask"`, since masking never removes a locus).}
#'   }
#' @examples
#' # A tiny hand-built VCF: locus_1's second sample has an ALT call backed
#' # by only 1 read.
#' vcf_lines <- c(
#'   "##fileformat=VCFv4.2",
#'   "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tind1\tind2",
#'   "un\t1\tlocus_1\tA\tC\t.\tPASS\t.\tGT:AD:DP\t0/0:10,0:10\t0/1:9,1:10"
#' )
#' vcf_file <- tempfile(fileext = ".vcf")
#' writeLines(vcf_lines, vcf_file)
#' H <- read_haps_vcf(vcf_file, verbose = FALSE)
#' res <- filter_low_conf_alt(H, min_alt_reads = 2, verbose = FALSE)
#' res$flagged_calls  # ind2's call at locus_1 is flagged (1 <= 2 ALT reads)
#' @export
filter_low_conf_alt <- function(H, min_alt_reads = 2, mode = "mask",
                                 drop_frac = 0, calls = NULL, verbose = TRUE) {
  if (!(length(mode) == 1L && mode %in% c("mask", "drop")))
    stop("mode must be exactly \"mask\" or \"drop\" (got: ", paste(mode, collapse = ", "), ").")
  if (is.null(calls)) calls <- .parse_alt_ad(H)

  flagged <- calls$usable & calls$AD_alt <= min_alt_reads
  n_rec <- nrow(H$A1)

  ## tabulate() again: counting, per LOCUS, how many ALT-containing calls it
  ## has in total, and how many of those got flagged. nbins = n_rec makes
  ## sure every locus gets an entry (even 0) rather than only the ones that
  ## happen to appear in `calls`.
  n_alt_calls <- tabulate(calls$record, nbins = n_rec)
  n_flagged   <- tabulate(calls$record[flagged], nbins = n_rec)
  flagged_fraction <- ifelse(n_alt_calls > 0L, n_flagged / n_alt_calls, 0)

  if (verbose)
    message(sprintf(
      "Low-confidence ALT filter (AD_alt <= %d): %s of %s ALT-containing calls flagged, ",
      min_alt_reads, format(sum(flagged), big.mark = ","), format(nrow(calls), big.mark = ",")),
      sprintf("affecting %s of %s loci", format(sum(n_flagged > 0L), big.mark = ","),
              format(n_rec, big.mark = ",")))

  ## Build the two summary tables from the ORIGINAL H (before any masking or
  ## dropping happens below), so they always describe what was found, even
  ## for a locus that's about to be removed.
  flagged_calls <- data.frame(
    locus = H$locus[calls$record[flagged]], sample = H$samples[calls$sample[flagged]],
    GT = calls$GT[flagged], AD_ref = calls$AD_ref[flagged], AD_alt = calls$AD_alt[flagged],
    DP = calls$DP[flagged], alt_fraction = calls$alt_fraction[flagged], row.names = NULL
  )

  if (mode == "mask") {
    ## Blank out (set to NA, R's "no data" value) just the flagged cells --
    ## cbind(rows, columns) here builds exactly the (locus, sample) address
    ## of each flagged call, so only those cells change.
    if (any(flagged)) {
      addr <- cbind(calls$record[flagged], calls$sample[flagged])
      H$A1[addr] <- NA
      H$A2[addr] <- NA
    }
    kept <- rep(TRUE, n_rec)
    locus_summary <- data.frame(
      locus = H$locus, n_alt_calls = n_alt_calls, n_flagged = n_flagged,
      flagged_fraction = round(flagged_fraction, 4), kept = kept, row.names = NULL
    )
  } else {
    ## drop mode: a locus survives if the SHARE of its calls that were
    ## flagged is not more than drop_frac (a locus with 0 ALT calls always
    ## has flagged_fraction 0 already, so it's always kept by this rule).
    kept <- flagged_fraction <= drop_frac + 1e-9
    locus_summary <- data.frame(
      locus = H$locus, n_alt_calls = n_alt_calls, n_flagged = n_flagged,
      flagged_fraction = round(flagged_fraction, 4), kept = kept, row.names = NULL
    )
    if (verbose)
      message(sprintf("  drop_frac = %.3g: %s of %s loci kept", drop_frac,
                      format(sum(kept), big.mark = ","), format(n_rec, big.mark = ",")))
    if (!sum(kept))
      stop("Every locus would be dropped at drop_frac = ", drop_frac, ". Raise drop_frac, ",
           "or lower min_alt_reads so fewer calls are flagged in the first place.")
    H <- .subset_H(H, kept)
  }

  list(H = H, flagged_calls = flagged_calls, locus_summary = locus_summary)
}

#' Sensitivity table for the low-confidence-ALT threshold
#'
#' Ports the Python script's `--sensitivity-table` option: for a range of
#' possible `min_alt_reads` thresholds, shows how many ALT-containing
#' genotype calls would be flagged. Useful for picking a threshold for
#' [filter_low_conf_alt()] before committing to one.
#'
#' @param H A list as returned by [read_haps_vcf()].
#' @param thresholds Which `min_alt_reads` values to try. Default
#'   `c(1, 2, 3, 4, 5, 10)`.
#' @param calls Optionally, a precomputed result from the internal
#'   `.parse_alt_ad()` helper (see [filter_low_conf_alt()]'s matching
#'   argument). Left `NULL` (the default), it's computed automatically.
#' @return A data frame with one row per threshold: `threshold, n_flagged,
#'   pct_flagged, total_usable` (the last two out of every ALT-containing
#'   call whose AD value could actually be read).
#' @examples
#' vcf_lines <- c(
#'   "##fileformat=VCFv4.2",
#'   "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tind1\tind2",
#'   "un\t1\tlocus_1\tA\tC\t.\tPASS\t.\tGT:AD:DP\t0/0:10,0:10\t0/1:9,1:10"
#' )
#' vcf_file <- tempfile(fileext = ".vcf")
#' writeLines(vcf_lines, vcf_file)
#' H <- read_haps_vcf(vcf_file, verbose = FALSE)
#' low_conf_alt_sensitivity(H)
#' @export
low_conf_alt_sensitivity <- function(H, thresholds = c(1, 2, 3, 4, 5, 10), calls = NULL) {
  if (is.null(calls)) calls <- .parse_alt_ad(H)
  usable_ad <- calls$AD_alt[calls$usable]
  total <- length(usable_ad)
  n_flagged <- vapply(thresholds, function(t) sum(usable_ad <= t), integer(1))
  data.frame(
    threshold = thresholds, n_flagged = n_flagged,
    pct_flagged = if (total) round(100 * n_flagged / total, 2) else rep(NA_real_, length(thresholds)),
    total_usable = total
  )
}
