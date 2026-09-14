###############################################################################
#
#  R/kinship.R -- screening every pair of individuals for close relatives,
#  BEFORE trusting a between-population comparison.
#
#  WHY THIS MATTERS. A population sample that includes a parent-offspring
#  pair or full siblings is not "n independent individuals"; it is fewer.
#  Every between-population statistic in this package (diversity_stats(),
#  het_between_pops()) assumes independent individuals, and the
#  het_between_pops() report tells the user to check relatedness first. This
#  file provides that check.
#
#  TWO METHODS, NEITHER ONE "the" answer. This function offers two published,
#  cited kinship estimators -- KING-robust (default) and Goudet's beta -- not
#  because they are known to be adequate, but because the most directly
#  relevant benchmark for reduced-representation SNP data in non-model
#  species (McMaster et al. 2025) found BOTH have real weaknesses: KING-robust
#  had high precision but SENSITIVITY BELOW 0.5 (it misses more real
#  relatives than it catches), while Goudet's beta had high sensitivity
#  (> 0.8) but limited precision, which got worse with population structure
#  and inbreeding. For general use that benchmark recommended a THIRD method,
#  PLINK's method-of-moments IBD (`PI_HAT`), as a balance of the two -- while
#  noting it was sensitive to filtering and often underestimated relatedness.
#  This package does not implement it. Read this function's results as "two
#  reasonable, cited screening options", not as a validated, sufficient check
#  -- and read a CLEAN result (no flagged pairs) with real caution given
#  KING's own known under-sensitivity, not as proof the dataset has no
#  relatives.
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
#  set. The formula is from the original source (Manichaikul et al. 2010)
#  and agrees with the re-derivation from a 9-cell genotype table by Waples,
#  Albrechtsen & Moltke (2019), which targets low-depth, RAD-like data. The
#  tests check that a duplicated individual reads 0.5. Only BIALLELIC loci
#  are used: KING's opposite-homozygote count is defined for two alleles.
#  (PLINK 2 handles multi-allelic sites as REF versus pooled ALT; that is not
#  implemented here.)
#
#  WHICH KING FORMULA. Manichaikul et al. give two: the one above (their
#  eq. 9, denominator N1_i + N1_j) and a "between-family" version whose
#  denominator uses the smaller of N1_i and N1_j. PLINK 2's --make-king and
#  Hail use the between-family version for every pair, so their numbers can
#  differ from this function's when two individuals differ in
#  heterozygosity. The sum form is kept deliberately: when one individual's
#  heterozygotes are under-called (allele dropout at low depth, common in
#  RAD data), it stays closer to the truth. In a simulation with 20% of one
#  individual's heterozygotes called homozygous, a parent-offspring pair
#  (true 0.25) read 0.226 with this form and 0.187 with the between-family
#  form, and an unrelated pair read -0.135 and -0.240.
#
#  KING's own published scale: identical/duplicate individuals average
#  ~0.5, parent-offspring or full siblings ~0.25, half-siblings/
#  grandparent-grandchild/avuncular ~0.125, unrelated ~0. The CATEGORY
#  BOUNDARIES (not those expected values) are the standard cutpoints:
#  ~0.177 (first-degree), ~0.0884 (second-degree), ~0.0442 (third-degree).
#  `threshold` defaults to 0.0442, the loosest of these. Using a first-degree
#  EXPECTED VALUE (0.25) as the threshold would miss about half of true
#  first-degree pairs, before KING's under-sensitivity on this kind of data
#  is even considered.
#
#  GOUDET'S BETA (`method = "beta"`, needs hierfstat) is the dosage-based
#  estimator behind `hierfstat::beta.dosage()`, restricted to the same
#  biallelic loci. hierfstat handles missing genotypes itself, so no
#  shared-locus count is kept for this method.
#
#  BETA NEEDS THE POPMAP. Beta measures how much more alike two individuals
#  are than the AVERAGE pair of individuals it was given. Given every
#  individual of several differentiated populations at once, two unrelated
#  members of the same population are more alike than that pooled average,
#  and read as relatives. In a simulation of two populations at FST = 0.1
#  (20 individuals each), 277 of 378 unrelated same-population pairs
#  exceeded 0.0442 when beta was computed on the pooled sample, and none did
#  when each population was its own reference -- while a planted full-sib
#  pair (0.259) and half-sib pair (0.138) were still found. So with `popmap`,
#  beta is computed within each population, which is also where relatives
#  are looked for; pairs from different populations get NA. KING does not
#  use allele frequencies, so it is computed on every pair either way.
#
#  Checked by tests/testthat/test-kinship.R.
#
###############################################################################

## Not exported. Restricts H to records OBSERVED to be biallelic (the same
## rule as write_plink(): the alleles actually seen, not the VCF's declared
## count, since a record can list an ALT allele nobody carries). Used by both
## methods below.
.restrict_biallelic <- function(H, verbose) {
  allele_stats <- locus_allele_stats(H)
  multiallelic <- allele_stats$n_observed_alleles > 2L
  if (any(multiallelic)) {
    .inform(verbose, sprintf("kinship_check(): excluding %s of %s multiallelic records (both methods need biallelic data)",
                             .big(sum(multiallelic)), .big(nrow(allele_stats))))
    H <- .subset_H(H, !multiallelic)
  }
  H
}

## Not exported. Renumbers each record's observed alleles as 1 and 2. Both
## methods below assume a biallelic record uses allele numbers 1 and 2 (KING
## looks for homozygotes of allele 1 and of allele 2; beta counts copies of
## allele 2). A record kept by .restrict_biallelic() can still use other
## numbers: in a haplotype VCF a record may declare three alleles of which
## only the 2nd and 3rd are carried by the screened samples, so its genotypes
## are 2/2, 2/3 and 3/3. Without renumbering, a 3/3 homozygote is neither
## "allele 1" nor "allele 2", opposite homozygotes go uncounted, and every
## pair reads as related. The smaller observed allele number becomes 1 and
## the larger becomes 2; a monomorphic record becomes all 1.
.recode_biallelic <- function(H) {
  if (!nrow(H$A1)) return(H)
  counts <- .allele_counts(H$A1, H$A2)                             # records x alleles
  observed <- counts > 0L
  ## Records with no genotype at all keep their (all-NA) genotypes.
  none <- rowSums(observed) == 0L
  observed[none, 1L] <- TRUE
  first <- max.col(observed, ties.method = "first")
  recode <- function(A) {
    out <- ifelse(A == first, 1L, 2L)      # `first` recycles down each column
    storage.mode(out) <- "integer"
    dimnames(out) <- dimnames(A)
    out
  }
  H$A1 <- recode(H$A1)
  H$A2 <- recode(H$A2)
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
.king_kinship <- function(H, chunk_rows = 20000L) {
  n_rec <- nrow(H$A1)
  n_samp <- ncol(H$A1)
  zero <- matrix(0, n_samp, n_samp)
  n_loci_used <- N11 <- N20 <- HetOtherCalled <- zero

  ## The counts are sums over loci, so they are accumulated a block of
  ## `chunk_rows` records at a time: the 0/1 matrices below then never exceed
  ## that many rows, which keeps memory bounded on large datasets.
  for (start in seq.int(1L, n_rec, by = chunk_rows)) {
    rows <- start:min(n_rec, start + chunk_rows - 1L)
    A1 <- H$A1[rows, , drop = FALSE]
    A2 <- H$A2[rows, , drop = FALSE]
    called <- !is.na(A1)                  # TRUE where this sample IS genotyped here
    het    <- called & (A1 != A2)         # TRUE where genotyped AND heterozygous
    ## Biallelic, renumbered 1/2 (.recode_biallelic()), so a homozygote is
    ## for allele 1 or allele 2 -- nothing else is possible.
    hom1 <- called & (A1 == A2) & (A1 == 1L)
    hom2 <- called & (A1 == A2) & (A1 == 2L)

    ## Coerce logical -> 0/1 numeric for matrix multiplication (R's %*% needs
    ## a numeric matrix, not TRUE/FALSE directly).
    as_01 <- function(x) matrix(as.numeric(x), length(rows), n_samp)
    Called <- as_01(called)
    Het    <- as_01(het)
    Hom1   <- as_01(hom1)
    Hom2   <- as_01(hom2)

    n_loci_used <- n_loci_used + crossprod(Called)                  # both-called loci
    N11 <- N11 + crossprod(Het)                                      # both heterozygous
    N20 <- N20 + crossprod(Hom1, Hom2) + crossprod(Hom2, Hom1)       # opposite homozygotes
    ## HetOtherCalled[i, j] = loci where i is heterozygous AND j is
    ## (separately) called there -- i's own heterozygous-locus count,
    ## restricted to the shared both-called set with j.
    HetOtherCalled <- HetOtherCalled + crossprod(Het, Called)
  }
  ## The TRANSPOSE of HetOtherCalled gives the same count from j's side, so
  ## adding the matrix to its own transpose gives the symmetric "N1_i + N1_j"
  ## KING's denominator needs.
  denom <- HetOtherCalled + t(HetOtherCalled)

  kinship <- (N11 - 2 * N20) / denom
  ## A pair where NEITHER individual is heterozygous at any shared locus has
  ## denom = 0, giving 0/0 = NaN or -x/0 = -Inf: no information, not a
  ## kinship value. Report NA, the package's usual "no data" code.
  kinship[!is.finite(kinship)] <- NA_real_
  dimnames(kinship) <- dimnames(n_loci_used) <- list(H$samples, H$samples)
  diag(kinship) <- NA_real_
  list(kinship = kinship, n_loci_used = n_loci_used)
}

## Not exported. Goudet's beta kinship matrix (hierfstat::beta.dosage()).
## Without `pop_of`, every sample is compared against the whole sample's
## average pair. With `pop_of` (population of each sample, named by sample),
## each population is its own reference and pairs from different populations
## are NA -- see "BETA NEEDS THE POPMAP" in the header for why.
.beta_kinship <- function(H, pop_of, verbose) {
  ## Dosage matrix, individuals in rows and loci in columns (as
  ## hierfstat::beta.dosage() expects): copies of allele 2 (ALT), 0/1/2, NA
  ## for a missing genotype.
  dosage <- t(H$A1 - 1L) + t(H$A2 - 1L)
  rownames(dosage) <- H$samples
  kinship <- matrix(NA_real_, length(H$samples), length(H$samples),
                    dimnames = list(H$samples, H$samples))
  if (is.null(pop_of)) {
    .inform(verbose, "kinship_check(): method = \"beta\" without `popmap` compares every pair with ",
            "the average pair of the WHOLE sample. If the sample holds more than one ",
            "population, unrelated members of the same population read as relatives; ",
            "pass `popmap` to compute beta within each population.")
    kinship[] <- hierfstat::beta.dosage(dosage, inb = FALSE)
  } else {
    for (p in unique(pop_of)) {
      members <- names(pop_of)[pop_of == p]
      if (length(members) < 3L) {
        .inform(verbose, sprintf("kinship_check(): population %s has %d individual(s); beta needs at least 3 as its own reference (kinship NA).",
                                 p, length(members)))
        next
      }
      kinship[members, members] <- hierfstat::beta.dosage(dosage[members, , drop = FALSE], inb = FALSE)
    }
  }
  diag(kinship) <- NA_real_
  kinship
}

#' Screen every pair of individuals for close relatives
#'
#' Estimates kinship for every pair of individuals, with KING-robust
#' (Manichaikul et al. 2010, the default, no extra package needed) or
#' Goudet's beta (`method = "beta"`, needs `hierfstat`). Both are screening
#' tools, not a validated final answer (see Details). Run it BEFORE
#' [het_between_pops()] or [diversity_stats()] and exclude one individual from
#' each flagged pair, so that those functions' assumption of independent
#' individuals holds.
#'
#' @details
#' **A screening tool, not a final answer.** The most relevant benchmark for
#' reduced-representation SNP data in non-model species (McMaster et al. 2025)
#' found that KING-robust has high precision but sensitivity below 0.5 on
#' such data (it misses more real relatives than it finds), while Goudet's
#' beta has high sensitivity but limited precision, worse under population
#' structure or inbreeding. For general use that benchmark recommended a
#' third method (PLINK's method-of-moments IBD, `PI_HAT`), as a balance of
#' the two, while noting it was sensitive to filtering and often
#' underestimated relatedness; this function does not implement it. Treat a
#' clean result (nothing flagged) with caution, not as proof that there are
#' no close relatives.
#'
#' **Relatives are looked for within populations.** Goudet's beta measures
#' how much more alike two individuals are than the average pair it is
#' given. Computed on several differentiated populations at once, unrelated
#' members of the same population look related (in a simulation at FST =
#' 0.1, about 70% of unrelated same-population pairs were flagged). With
#' `popmap`, beta is computed within each population, each population being
#' its own reference, and pairs from different populations are `NA`. In a
#' small population, or one with many relatives, that reference itself
#' contains relatedness, so beta then understates kinship. KING uses no
#' allele frequencies and is computed on every pair either way; a pair from
#' different populations flagged by KING is worth checking as a possible
#' sample mix-up.
#'
#' **Which KING formula.** This is Manichaikul et al.'s (2010) estimator
#' with the SUM of the two individuals' heterozygote counts in the
#' denominator. PLINK 2 (`--make-king`) uses their between-family version,
#' with the SMALLER count, so values can differ when two individuals differ
#' in heterozygosity. The sum form is less affected by heterozygotes that
#' were called homozygous in a low-coverage library (see `R/kinship.R`).
#'
#' **KING cutpoints**, for choosing `threshold`: ~0.354 (duplicate/identical
#' twin), ~0.177 (first degree: parent-offspring or full sibling), ~0.0884
#' (second degree: half-sibling, grandparent-grandchild, avuncular), ~0.0442
#' (third degree). The default, 0.0442, is the most inclusive.
#'
#' **`n_loci_used` matters as much as the kinship value.** With uneven RAD-seq
#' missing data, two pairs with the same kinship can rest on very different
#' amounts of evidence. Pairs with fewer than `min_shared_loci` loci genotyped
#' in both individuals get `NA`, and pairs with fewer than
#' `low_confidence_loci` are marked `low_confidence`.
#'
#' **References.** Manichaikul, A., Mychaleckyj, J.C., Rich, S.S., Daly, K.,
#' Sale, M. & Chen, W.-M. (2010) Robust relationship inference in
#' genome-wide association studies. *Bioinformatics* 26:2867-2873.
#' \doi{10.1093/bioinformatics/btq559} -- Waples, R.K., Albrechtsen, A. &
#' Moltke, I. (2019) Allele frequency-free inference of close familial
#' relationships from genotypes or low-depth sequencing data. *Molecular
#' Ecology* 28:35-48. \doi{10.1111/mec.14954} -- McMaster, E.S. et al. (2025)
#' Evaluating kinship estimation methods for reduced-representation SNP data
#' in non-model species. *Molecular Ecology Resources*.
#' \doi{10.1111/1755-0998.70038} --
#' Goudet, J., Kay, T. & Weir, B.S. (2018) How to estimate kinship.
#' *Molecular Ecology* 27:4121-4135. (`method = "beta"`,
#' `hierfstat::beta.dosage()`.) --
#' Weir, B.S. & Goudet, J. (2017) A unified characterization of population
#' structure and relatedness. *Genetics* 206:2085-2103.
#' \doi{10.1534/genetics.116.198424} --
#' Chang, C.C., Chow, C.C., Tellier, L.C.A.M., Vattikuti, S., Purcell, S.M. &
#' Lee, J.J. (2015) Second-generation PLINK: rising to the challenge of larger
#' and richer datasets. *GigaScience* 4:7. (PLINK 2's `--make-king`.)
#'
#' @param vcf Path to a VCF file, or the object returned by [read_stacks_vcf()]
#'   (optionally filtered).
#' @param popmap Optional: a popmap file path or the list returned by
#'   [read_popmap()]. Only its samples are screened, and each pair is
#'   labeled with the samples' populations. Needed for `method = "beta"` on
#'   more than one population (see Details). Default `NULL`.
#' @param method `"king"` (default, no extra package needed) or `"beta"`
#'   (Goudet's beta via `hierfstat::beta.dosage()`; needs hierfstat).
#' @param threshold Kinship above which a pair is listed in `flagged_pairs`.
#'   Default `0.0442` (third degree or closer; see Details). `NULL` lists every
#'   pair that has a kinship value.
#' @param min_shared_loci Pairs with fewer loci genotyped in both individuals
#'   get kinship `NA`: too little data to report a number. Default `30`.
#' @param low_confidence_loci Pairs with fewer shared loci than this are
#'   marked `low_confidence = TRUE`. Default `200`.
#' @param outdir If given, write `kinship_pairwise.tsv` there. Default `NULL`.
#' @param verbose Print progress and a summary. Default `TRUE`.
#' @return A list:
#'   \describe{
#'     \item{pairwise}{One row per pair: `sample1`, `sample2`, (with `popmap`)
#'       `population1`, `population2`, then `kinship` (`NA` below
#'       `min_shared_loci`), `n_loci_used`, `low_confidence`.}
#'     \item{flagged_pairs}{The rows of `pairwise` with `kinship > threshold`
#'       (never an `NA` kinship).}
#'   }
#' @examples
#' H <- read_stacks_vcf(system.file("extdata", "small.snps.vcf",
#'                                  package = "RADdiversity"), verbose = FALSE)
#' # The toy data have few loci, so lower the evidence floors for this example.
#' popmap <- system.file("extdata", "small_popmap.tsv", package = "RADdiversity")
#' kin <- kinship_check(H, popmap, min_shared_loci = 5, low_confidence_loci = 10)
#' kin$pairwise
#' kin$flagged_pairs
#' @export
kinship_check <- function(vcf, popmap = NULL, method = "king", threshold = 0.0442,
                          min_shared_loci = 30L, low_confidence_loci = 200L, outdir = NULL,
                          verbose = TRUE) {
  .check_choice(method, "method", c("king", "beta"))
  if (!is.null(threshold)) .check_number(threshold, "threshold")
  min_shared_loci <- .check_count(min_shared_loci, "min_shared_loci")
  low_confidence_loci <- .check_count(low_confidence_loci, "low_confidence_loci")
  if (!is.null(outdir)) .check_string(outdir, "outdir")
  .check_flag(verbose, "verbose")
  if (is.character(popmap) && length(popmap) == 1L && !file.exists(popmap))
    stop("Popmap file not found: ", popmap, "\n  Check the path and try again.", call. = FALSE)
  if (method == "beta" && !.hierfstat_available())
    stop("method = \"beta\" needs the hierfstat package: install.packages(\"hierfstat\")",
         call. = FALSE)

  H <- .resolve_H(vcf, verbose = verbose)
  ## With a popmap, only its samples are screened, and each is labeled with
  ## its population.
  pop_of <- NULL
  if (!is.null(popmap)) {
    pops <- .resolve_pops(popmap, H$samples, verbose = verbose)
    ids <- unlist(pops, use.names = FALSE)
    H$A1 <- H$A1[, ids, drop = FALSE]
    H$A2 <- H$A2[, ids, drop = FALSE]
    H$samples <- ids
    pop_of <- stats::setNames(rep(names(pops), lengths(pops)), ids)
  }
  n_samp <- length(H$samples)
  if (n_samp < 2) stop("Need at least 2 samples to compute any kinship.", call. = FALSE)
  H <- .restrict_biallelic(H, verbose)
  if (!nrow(H$A1))
    stop("No biallelic record remains after excluding multiallelic records -- kinship ",
         "can't be computed.", call. = FALSE)
  H <- .recode_biallelic(H)

  if (method == "king") {
    king <- .king_kinship(H)
    kinship <- king$kinship
    shared_loci <- king$n_loci_used
  } else {
    kinship <- .beta_kinship(H, pop_of, verbose)
    ## beta.dosage() has no per-pair count of shared loci, so every pair is
    ## given the total number of records: this never UNDERstates the data a
    ## pair rests on, but cannot mark thin pairs either.
    shared_loci <- matrix(nrow(H$A1), n_samp, n_samp, dimnames = list(H$samples, H$samples))
    diag(shared_loci) <- NA_integer_
  }

  ## One row per pair (upper triangle: kinship is symmetric).
  upper <- which(upper.tri(kinship), arr.ind = TRUE)
  pairwise <- data.frame(sample1 = H$samples[upper[, 1]], sample2 = H$samples[upper[, 2]],
                         kinship = kinship[upper], n_loci_used = as.integer(shared_loci[upper]),
                         row.names = NULL)
  if (!is.null(pop_of)) {
    pairwise <- data.frame(pairwise[c("sample1", "sample2")],
                           population1 = unname(pop_of[pairwise$sample1]),
                           population2 = unname(pop_of[pairwise$sample2]),
                           pairwise[c("kinship", "n_loci_used")], row.names = NULL)
  }
  too_thin <- pairwise$n_loci_used < min_shared_loci
  pairwise$kinship[too_thin] <- NA_real_
  pairwise$low_confidence <- pairwise$n_loci_used < low_confidence_loci

  .inform(verbose, sprintf(
    "kinship_check() [method = \"%s\"]: %s pairs, %s with < %d shared loci (kinship NA), %s low-confidence (< %d shared loci)",
    method, .big(nrow(pairwise)), .big(sum(too_thin)), min_shared_loci,
    .big(sum(pairwise$low_confidence & !too_thin)), low_confidence_loci))

  has_value <- !is.na(pairwise$kinship)
  flagged_pairs <- if (is.null(threshold)) pairwise[has_value, , drop = FALSE]
                   else pairwise[has_value & pairwise$kinship > threshold, , drop = FALSE]
  if (!is.null(threshold)) {
    if (nrow(flagged_pairs))
      .inform(verbose, sprintf("  %d pair(s) exceed threshold = %.4f -- see flagged_pairs. Consider excluding\n  one individual from each flagged pair before het_between_pops()/diversity_stats().",
                               nrow(flagged_pairs), threshold))
    else
      .inform(verbose, sprintf("  No pair exceeds threshold = %.4f. See ?kinship_check before treating this as \"no relatives\":\n  both methods miss real relatives in reduced-representation data.", threshold))
  }

  if (!is.null(outdir)) {
    written <- .write_tables(list(pairwise = .round_table(pairwise)), outdir,
                             c(pairwise = "kinship_pairwise.tsv"))
    .inform(verbose, "  Wrote ", written)
  }
  list(pairwise = pairwise, flagged_pairs = flagged_pairs)
}
