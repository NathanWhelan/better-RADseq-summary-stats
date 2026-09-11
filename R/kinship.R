###############################################################################
#
#  R/kinship.R -- screening every pair of individuals for close relatives,
#  BEFORE trusting a between-population comparison.
#
#  WHY THIS MATTERS. het_between_pops() (R/het_between_pops.R) already warns,
#  in its own printed output, to "check relatedness first" -- a population
#  sample that happens to include a parent-offspring pair or a set of full
#  siblings is not really "n independent individuals", it is fewer than that,
#  and every between-population statistic in this package (here and in
#  diversity_stats()/het_between_pops()) assumes independent individuals.
#  Nothing in this package actually performed that check until now.
#
#  TWO METHODS, NEITHER ONE "the" answer. This function offers two published,
#  cited kinship estimators -- KING-robust (default) and Goudet's beta -- not
#  because they are known to be adequate, but because the most directly
#  relevant benchmark for reduced-representation SNP data in non-model
#  species (McMaster et al. 2025) found BOTH have real weaknesses: KING-robust
#  has high precision but SENSITIVITY BELOW 0.5 on data like this (it misses
#  more real relatives than it catches), while Goudet's beta has better
#  sensitivity but degrades under real population structure or inbreeding.
#  That same benchmark's own top overall pick for default use was a THIRD
#  method (classic PLINK moment-of-moments IBD, `PI_HAT`) that this package
#  implements NEITHER of. Read this function's results as "two reasonable,
#  cited screening options", not as a validated, sufficient check -- and read
#  a CLEAN result (no flagged pairs) with real caution given KING's own
#  known under-sensitivity, not as proof the dataset has no relatives.
#
#  KING-ROBUST (Manichaikul et al. 2010) is the default because it is
#  closed-form and needs no allele-frequency estimate (so it isn't thrown off
#  by pooling structured populations to guess allele frequencies, unlike many
#  older kinship estimators) -- for a pair of individuals i,j, restricted to
#  loci where BOTH are genotyped:
#
#      phi_ij = (N11 - 2*N20) / (N1_i + N1_j)
#
#  N11 = loci where BOTH i and j are heterozygous; N20 = loci where they are
#  OPPOSITE homozygotes (share no allele -- only possible at a biallelic
#  locus if one is hom-REF and the other hom-ALT); N1_i / N1_j = i's / j's
#  own heterozygous-locus count, over that SAME shared both-called locus
#  set. Formula confirmed against its ORIGINAL source (Manichaikul et al.
#  2010) and against an independent re-derivation of the same statistic via
#  a 9-cell genotype contingency table (Waples, Albrechtsen & Moltke 2019,
#  more directly on-point since it targets low-depth/RADseq-style data
#  specifically -- full citations in @references below), and cross-checked
#  for internal consistency (a duplicate individual compared to itself
#  reads ~0.5; two individuals sharing no alleles at all read very
#  negative) in this package's own self-test. Restricted to BIALLELIC loci
#  for this v1 (KING's opposite-homozygote logic is defined per-diallelic-
#  site; PLINK2's own multiallelic handling via REF-vs-pooled-ALT is a real
#  option, deliberately deferred rather than guessed at here).
#
#  KING's own published scale: identical/duplicate individuals average
#  ~0.5, parent-offspring or full siblings ~0.25, half-siblings/
#  grandparent-grandchild/avuncular ~0.125, unrelated ~0. The CATEGORY
#  BOUNDARIES (not those expected values) are the standard cutpoints:
#  ~0.177 (first-degree), ~0.0884 (second-degree), ~0.0442 (third-degree).
#  `threshold` defaults to 0.0442, the loosest of these -- a first-degree
#  EXPECTED VALUE of 0.25 would already miss roughly half of true
#  first-degree pairs even before KING's own extra under-sensitivity on
#  this kind of data compounds it.
#
#  GOUDET'S BETA (`method = "beta"`, needs hierfstat) is the dosage-based
#  estimator behind `hierfstat::beta.dosage()`, restricted to the same
#  biallelic loci. hierfstat's own `matching()` (which `beta.dosage()` calls)
#  already tolerates missing genotypes internally, so no separate shared-
#  locus-count bookkeeping is needed for this method the way KING needs it.
#
###############################################################################

## Not exported. Restricts H to loci that are OBSERVED to be biallelic (same
## rule write_plink() already uses: locus_allele_stats()$n_observed_alleles,
## not the VCF's own declared n_alleles, since a locus can list an ALT allele
## nobody actually carries). Shared by both methods below.
.restrict_biallelic <- function(H, verbose) {
  stats <- locus_allele_stats(H)
  bad <- stats$n_observed_alleles > 2L
  if (any(bad)) {
    if (verbose)
      message(sprintf("kinship_check(): excluding %s of %s multiallelic loci (both methods need biallelic data)",
                      format(sum(bad), big.mark = ","), format(nrow(stats), big.mark = ",")))
    H <- .subset_H(H, !bad)
  }
  H
}

## Not exported. KING-robust (Manichaikul et al. 2010), vectorised over every
## pair of samples at once via matrix multiplication rather than an
## individual-by-individual loop (which would be needlessly slow once a
## dataset has more than a couple hundred individuals -- pairs grow with the
## SQUARE of sample count). Returns a list of n x n matrices (kinship,
## n_loci_used), `NA` on the diagonal (self-kinship isn't a meaningful
## question here).
##
## HOW THE MATRIX-MULTIPLICATION TRICK WORKS, IN PLAIN TERMS: for a 0/1
## indicator matrix X with one row per locus and one column per sample,
## `t(X) %*% X` is a sample-by-sample matrix whose (i, j) entry is exactly
## "how many loci have that indicator TRUE in BOTH sample i and sample j" --
## ordinary matrix multiplication sums exactly that product over loci. This
## is the same computation as a triple loop over (locus, sample_i, sample_j)
## would do, just handed to R's fast underlying matrix-multiply routine
## instead of looping row by row in R itself.
.king_kinship <- function(H) {
  A1 <- H$A1; A2 <- H$A2
  n_rec <- nrow(A1); n_samp <- ncol(A1)
  called <- !is.na(A1)                    # TRUE where this sample IS genotyped here
  het    <- called & (A1 != A2)           # TRUE where genotyped AND heterozygous
  ## Biallelic, so a homozygote is for allele 1 (REF) or allele 2 (ALT) --
  ## nothing else is possible once multiallelic loci have been excluded.
  hom1 <- called & (A1 == A2) & (A1 == 1L)
  hom2 <- called & (A1 == A2) & (A1 == 2L)

  ## Coerce logical -> 0/1 numeric for matrix multiplication (R's %*% needs a
  ## numeric matrix, not TRUE/FALSE directly).
  Called <- matrix(as.numeric(called), n_rec, n_samp)
  Het    <- matrix(as.numeric(het),    n_rec, n_samp)
  Hom1   <- matrix(as.numeric(hom1),   n_rec, n_samp)
  Hom2   <- matrix(as.numeric(hom2),   n_rec, n_samp)

  n_loci_used <- crossprod(Called)                    # [i,j] = both-called locus count
  N11 <- crossprod(Het)                                # [i,j] = both-heterozygous locus count
  N20 <- crossprod(Hom1, Hom2) + crossprod(Hom2, Hom1) # [i,j] = opposite-homozygote locus count
  ## HetOtherCalled[i, j] = loci where i is heterozygous AND j is (separately)
  ## called there -- i's own heterozygous-locus count, restricted to the
  ## shared both-called set with j. Its TRANSPOSE gives the same thing from
  ## j's side, so adding a matrix to its own transpose gives the symmetric
  ## "N1_i + N1_j" KING's denominator needs.
  HetOtherCalled <- crossprod(Het, Called)
  denom <- HetOtherCalled + t(HetOtherCalled)

  kinship <- (N11 - 2 * N20) / denom
  dimnames(kinship) <- dimnames(n_loci_used) <- list(H$samples, H$samples)
  diag(kinship) <- NA_real_
  list(kinship = kinship, n_loci_used = n_loci_used)
}

#' Screen every pair of individuals for close relatives
#'
#' Estimates pairwise kinship for every pair of individuals in `vcf_file`,
#' using either KING-robust (Manichaikul et al. 2010, the default -- needs no
#' extra package) or Goudet's beta (`method = "beta"`, needs `hierfstat`).
#' Both methods are screening tools, not a validated final answer -- see
#' `@details` and the header comment in `R/kinship.R` for why. Intended to be
#' run BEFORE [het_between_pops()] or [diversity_stats()], excluding one
#' individual from each flagged pair so those functions' "independent
#' individuals" assumption actually holds.
#'
#' @details
#' **This is a screening tool, not a validated final answer.** The most
#' directly relevant benchmark for reduced-representation SNP data in
#' non-model species (McMaster et al. 2025) found KING-robust has high
#' precision but sensitivity below 0.5 on this kind of data (it misses more
#' real relatives than it flags), while Goudet's beta has better sensitivity
#' but degrades under real population structure or inbreeding. That
#' benchmark's own overall top pick for default use was a THIRD method
#' (classic PLINK moment-of-moments IBD, `PI_HAT`) that this function
#' implements neither of. Treat a clean result (nothing flagged) with real
#' caution, not as proof the dataset has no close relatives.
#'
#' **The KING cutpoint table**, for setting `threshold` deliberately rather
#' than accepting the loosest default: ~0.354 (duplicate/identical twin),
#' ~0.177 (first-degree: parent-offspring or full sibling), ~0.0884
#' (second-degree: half-sibling, grandparent-grandchild, avuncular), ~0.0442
#' (third-degree). `threshold` defaults to 0.0442, the loosest (most
#' inclusive) of these.
#'
#' **`n_loci_used` matters as much as the kinship value.** Two pairs at the
#' same kinship value can rest on very different amounts of evidence once
#' RADseq missingness is uneven -- pairs with fewer than 30 shared
#' both-called loci get `NA` instead of a number (too little data to trust
#' at all), and pairs under 200 are still reported but noted as
#' low-confidence, rather than silently trusted the same as a
#' well-supported pair.
#'
#' **References.** Manichaikul, A., Mychaleckyj, J.C., Rich, S.S., Daly, K.,
#' Sale, M. & Chen, W.-M. (2010) Robust relationship inference in
#' genome-wide association studies. *Bioinformatics* 26:2867-2873.
#' <https://doi.org/10.1093/bioinformatics/btq559> (KING-robust, primary
#' source) -- Waples, R.K., Albrechtsen, A. & Moltke, I. (2019) Allele
#' frequency-free inference of close familial relationships from genotypes
#' or low-depth sequencing data. *Molecular Ecology* 28:35-48.
#' <https://doi.org/10.1111/mec.14954> (independent re-derivation of the
#' same statistic, more directly on-point for RADseq/low-depth data) --
#' McMaster, E.S. et al. (2025) Evaluating kinship estimation methods for
#' reduced-representation SNP data in non-model species. *Molecular
#' Ecology Resources*. <https://doi.org/10.1111/1755-0998.70038> (the
#' sensitivity/precision benchmark this function's `@details` cites
#' above).
#'
#' @param vcf_file Path to a Stacks VCF, or an already-parsed `H` list (see
#'   [diversity_stats()] for the same option there).
#' @param method `"king"` (default, no extra package needed) or `"beta"`
#'   (Goudet's beta via `hierfstat::beta.dosage()`, needs hierfstat
#'   installed).
#' @param threshold Kinship value above which a pair is reported in
#'   `flagged_pairs`. Default `0.0442` (the third-degree-or-closer KING
#'   cutpoint -- see `@details`).
#' @param outdir If given, write `kinship_pairwise.tsv` there. Default `NULL`
#'   (no file written).
#' @param verbose Print progress and the locus-count summary. Default `TRUE`.
#' @return Invisibly, a list:
#'   \describe{
#'     \item{pairwise}{Data frame, one row per pair: `sample1`, `sample2`,
#'       `kinship` (`NA` if `n_loci_used` is below 30), `n_loci_used`,
#'       `low_confidence` (`TRUE` below 200 shared loci).}
#'     \item{flagged_pairs}{The subset of `pairwise` with `kinship >
#'       threshold` (and therefore never the `NA`-kinship rows).}
#'   }
#' @examples
#' # 3 samples: s1 and s2 are constructed as identical (every genotype
#' # matches -- KING should read this pair at exactly 0.5), s3 differs.
#' # The 6-locus pattern below is repeated 6x (36 loci total) purely so this
#' # toy example clears the "at least 30 shared loci" floor described above
#' # and actually prints a number instead of NA -- a real dataset would not
#' # need this repetition trick.
#' samp <- c("s1", "s2", "s3")
#' base1 <- matrix(c(1,1,1, 1,1,2, 2,2,1, 1,1,2, 2,2,1, 1,1,2), nrow = 6, byrow = TRUE)
#' base2 <- matrix(c(1,1,2, 2,2,1, 2,2,2, 2,2,1, 2,2,2, 2,2,1), nrow = 6, byrow = TRUE)
#' n_loc <- 36
#' H <- list(
#'   A1 = `dimnames<-`(do.call(rbind, replicate(6, base1, simplify = FALSE)), list(NULL, samp)),
#'   A2 = `dimnames<-`(do.call(rbind, replicate(6, base2, simplify = FALSE)), list(NULL, samp)),
#'   locus = paste0("locus_", seq_len(n_loc)), locus_raw = paste0("locus_", seq_len(n_loc)),
#'   n_alleles = rep(2L, n_loc), alleles = replicate(n_loc, c("A", "C"), simplify = FALSE),
#'   samples = samp
#' )
#' res <- kinship_check(H, verbose = FALSE)
#' res$pairwise  # s1-s2 (identical) reads 0.5; s1-s3/s2-s3 read noticeably lower
#' @export
kinship_check <- function(vcf_file, method = "king", threshold = 0.0442,
                           outdir = NULL, verbose = TRUE) {
  if (!(length(method) == 1L && method %in% c("king", "beta")))
    stop("method must be exactly \"king\" or \"beta\" (got: ", paste(method, collapse = ", "), ").")
  if (!is.null(threshold) && (length(threshold) != 1L || is.na(threshold)))
    stop("threshold must be a single number (or NULL to report every pair).")

  if (is.character(vcf_file)) message("Reading ", vcf_file, " ...")
  H <- .resolve_H(vcf_file, verbose = verbose)
  n_samp <- length(H$samples)
  if (n_samp < 2) stop("Need at least 2 samples to compute any kinship.")
  H <- .restrict_biallelic(H, verbose)
  if (!nrow(H$A1))
    stop("No biallelic locus remains after excluding multiallelic loci -- kinship can't be computed.")

  if (method == "king") {
    res <- .king_kinship(H)
    kin_mat <- res$kinship; n_mat <- res$n_loci_used
  } else {
    ## .hierfstat_available() (defined in R/write_formats.R) is a thin
    ## wrapper around requireNamespace() reused here for the same reason it
    ## exists there: it lives in THIS package's own namespace, so tests can
    ## mock it via testthat::with_mocked_bindings() to exercise this path
    ## without actually needing hierfstat uninstalled -- mocking base R's
    ## own requireNamespace() directly isn't possible.
    if (!.hierfstat_available())
      stop("method = \"beta\" needs the hierfstat package: install.packages(\"hierfstat\")")
    ## Dosage matrix: samples x loci (hierfstat::matching()'s own required
    ## orientation, confirmed against its source this session), values 0/1/2
    ## = copies of allele 2 (ALT), NA for a missing genotype -- hierfstat's
    ## matching() already tolerates NA internally, so no extra bookkeeping
    ## is needed here the way KING's shared-locus-count logic needs.
    dos <- t(H$A1 - 1L) + t(H$A2 - 1L)
    beta <- hierfstat::beta.dosage(dos, inb = FALSE)
    dimnames(beta) <- list(H$samples, H$samples)
    diag(beta) <- NA_real_
    kin_mat <- beta
    ## beta.dosage() has no notion of "shared locus count" the way KING's
    ## formula needs one (hierfstat::matching() folds missingness into its
    ## own weighting instead) -- report every locus in H as available to
    ## every pair, which is conservative (never UNDER-states how little
    ## data a thin pair actually rests on, only fails to flag it as such).
    n_mat <- matrix(nrow(H$A1), n_samp, n_samp, dimnames = list(H$samples, H$samples))
    diag(n_mat) <- NA_integer_
  }

  ## Long-format table: one row per pair, upper triangle only (kinship is
  ## symmetric, so the lower triangle would just repeat the same pairs).
  ut <- which(upper.tri(kin_mat), arr.ind = TRUE)
  pairwise <- data.frame(
    sample1 = H$samples[ut[, 1]], sample2 = H$samples[ut[, 2]],
    kinship = round(kin_mat[ut], 4),
    n_loci_used = as.integer(n_mat[ut]),
    row.names = NULL)
  ## Below 30 shared loci: too little data to report a number at all (an
  ## explicit NA, not a suspicious-looking value that happens to be small).
  ## Below 200: still reported, but flagged so a low_confidence result isn't
  ## silently treated the same as a well-supported one -- see @details.
  too_thin <- pairwise$n_loci_used < 30L
  pairwise$kinship[too_thin] <- NA_real_
  pairwise$low_confidence <- pairwise$n_loci_used < 200L

  if (verbose) {
    message(sprintf("kinship_check() [method = \"%s\"]: %s pairs, %s with < 30 shared loci (reported as NA), %s low-confidence (< 200 shared loci)",
                    method, format(nrow(pairwise), big.mark = ","), format(sum(too_thin), big.mark = ","),
                    format(sum(pairwise$low_confidence & !too_thin), big.mark = ",")))
  }

  flagged_pairs <- if (is.null(threshold)) pairwise[0, ]
    else pairwise[!is.na(pairwise$kinship) & pairwise$kinship > threshold, , drop = FALSE]
  if (verbose && !is.null(threshold)) {
    if (nrow(flagged_pairs))
      message(sprintf("  %d pair(s) exceed threshold = %.4f -- see flagged_pairs. Consider excluding\n  one individual from each flagged pair before het_between_pops()/diversity_stats().",
                      nrow(flagged_pairs), threshold))
    else
      message(sprintf("  No pair exceeds threshold = %.4f. See @details before treating this as \"no relatives\" --\n  read the sensitivity caveat for whichever method you used.", threshold))
  }

  if (!is.null(outdir)) {
    dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
    f <- file.path(outdir, "kinship_pairwise.tsv")
    utils::write.table(pairwise, f, sep = "\t", quote = FALSE, row.names = FALSE)
    if (verbose) message("  Wrote ", f)
  }

  invisible(list(pairwise = pairwise, flagged_pairs = flagged_pairs))
}
