###############################################################################
#
#  R/adegenet.R -- the data object as an adegenet `genind` or `genlight`, so
#  a filtered dataset can go straight into adegenet, dartR, poppr and the
#  other packages that read those classes.
#
#  adegenet is only suggested, not required: these two functions check for it
#  and say how to install it. Both follow the write_*() functions' rules
#  (R/write_formats.R): `popmap` is optional, and when it is given every
#  sample in H must be in it (filter_samples() removes the others).
#
#    as_genind()    any marker, including multi-allelic RAD haplotypes. The
#                   individuals x alleles count table is built here and handed
#                   to adegenet::genind(), which is much faster than
#                   adegenet::df2genind() on large data.
#    as_genlight()  biallelic records only (SNPs), stored as the count of
#                   the second allele (0, 1, 2), which is what dartR uses.
#
###############################################################################

## Not exported. TRUE when adegenet is installed (a wrapper so tests can
## pretend it is not; see .hierfstat_available() in R/utils.R).
.adegenet_available <- function() requireNamespace("adegenet", quietly = TRUE)

## Not exported. Stops with install instructions when adegenet is missing.
.need_adegenet <- function(what) {
  if (!.adegenet_available())
    stop(what, " needs the adegenet package: install.packages(\"adegenet\")", call. = FALSE)
  invisible(NULL)
}

## Not exported. The population factor for the writers' rule: NULL without a
## popmap; otherwise each sample's population, levels in popmap order. Stops if
## a sample of H is not in the popmap (see .pop_id_vector()).
.adegenet_pop <- function(H, popmap, what) {
  pops <- .writer_pops(H, popmap, required = FALSE, what)
  if (is.null(pops)) return(NULL)
  factor(names(pops)[.pop_id_vector(H, pops)], levels = names(pops))
}

## Not exported. Keeps only records where at least one individual is
## genotyped (a record with no genotype has no allele to count), with a
## message saying how many were dropped.
.drop_empty_records <- function(H, verbose, what) {
  typed <- rowSums(!is.na(H$A1)) > 0L
  if (!all(typed)) {
    .inform(verbose, sprintf("%s: leaving out %s record(s) with no genotyped individual.",
                             what, .big(sum(!typed))))
    H <- .subset_H(H, typed)
  }
  if (!nrow(H$A1)) stop(what, ": no record has a genotyped individual.", call. = FALSE)
  H
}

#' Convert to an adegenet genind or genlight object
#'
#' Turns the data (a VCF path, or the object from [read_stacks_vcf()],
#' optionally filtered) into an object of the adegenet package, so it can be
#' used with adegenet, dartR, poppr and other packages that read these classes.
#'
#' * `as_genind()` works for any marker, including multi-allelic RAD
#'   haplotypes from `populations.haps.vcf`. Each allele is named by its
#'   sequence. It stores one number per individual and allele, so for a large
#'   SNP dataset `as_genlight()` uses much less memory.
#' * `as_genlight()` is for biallelic records (SNPs). Each genotype is stored
#'   as the number of copies of the record's second allele (the ALT allele of a
#'   normal SNP): 0, 1 or 2. For a normal SNP the second allele is always the
#'   ALT, even if only one allele is seen, so objects made from different
#'   subsets of samples agree. A record with more than two observed alleles
#'   stops the function unless `drop_multiallelic = TRUE`, as in
#'   [write_plink()].
#'
#' Both leave out records with no genotyped individual, and write a missing
#' genotype as `NA`.
#'
#' **Locus names.** adegenet does not allow a `.` in a genind locus name, and
#' this package names the SNPs of one RAD locus `1`, `1.1`, `1.2`, ... So in
#' `as_genind()` every `.` becomes `_` (`1_1`, `1_2`). The original names, and
#' each record's CHROM and POS, are kept in `other(x)`. `as_genlight()` keeps
#' the names as they are, and also sets `chromosome()` and `position()`, and
#' `alleles()` when every allele is a single base.
#'
#' **Which samples.** Every sample in the data is converted. With `popmap`,
#' each sample gets its population (`pop()`), and every sample must be in the
#' popmap; remove the others first with [filter_samples()].
#'
#' @references
#' Jombart, T. (2008) adegenet: a R package for the multivariate analysis of
#' genetic markers. *Bioinformatics* 24:1403-1405.
#'
#' Jombart, T. & Ahmed, I. (2011) adegenet 1.3-1: new tools for the analysis of
#' genome-wide SNP data. *Bioinformatics* 27:3070-3071.
#'
#' @param H A VCF path, or the object returned by [read_stacks_vcf()]
#'   (optionally filtered).
#' @param popmap Optional: a popmap file path or the list returned by
#'   [read_popmap()], to set each individual's population.
#' @param verbose Print what was left out. Default `TRUE`.
#' @return `as_genind()`: an adegenet `genind` object. `as_genlight()`: an
#'   adegenet `genlight` object.
#' @examples
#' if (requireNamespace("adegenet", quietly = TRUE)) {
#'   popmap <- system.file("extdata", "small_popmap.tsv", package = "RADdiversity")
#'   haps <- read_stacks_vcf(system.file("extdata", "small.haps.vcf",
#'                                       package = "RADdiversity"), verbose = FALSE)
#'   gi <- as_genind(haps, popmap)          # haplotype loci, many alleles
#'   gi
#'   snps <- read_stacks_vcf(system.file("extdata", "small.snps.vcf",
#'                                       package = "RADdiversity"), verbose = FALSE)
#'   gl <- as_genlight(snps, popmap)        # SNPs as 0/1/2
#'   gl
#' }
#' @export
as_genind <- function(H, popmap = NULL, verbose = TRUE) {
  .need_adegenet("as_genind()")
  .check_flag(verbose, "verbose")
  H <- .resolve_H(H, verbose = FALSE)
  pop <- .adegenet_pop(H, popmap, "as_genind()")
  H <- .drop_empty_records(H, verbose, "as_genind()")
  n_rec <- nrow(H$A1)
  k <- max(1L, H$n_alleles)

  ## Which alleles each record actually shows, and the column each one gets:
  ## records in order, and within a record its alleles in VCF order.
  observed <- .allele_counts(H$A1, H$A2, k) > 0L              # records x alleles
  ## Number the observed (record, allele) cells record by record: which() on
  ## the transpose walks the records in order.
  by_record <- matrix(0L, k, n_rec)
  by_record[which(t(observed))] <- seq_len(sum(observed))
  column <- t(by_record)                                       # records x alleles

  ## Allele copies of each individual: one column per (record, allele).
  tab <- matrix(NA_real_, ncol(H$A1), sum(observed), dimnames = list(H$samples, NULL))
  for (a in seq_len(k)) {
    rows <- which(observed[, a])
    if (!length(rows)) next
    copies <- (H$A1[rows, , drop = FALSE] == a) + (H$A2[rows, , drop = FALSE] == a)
    tab[, column[rows, a]] <- t(copies)                        # NA where missing
  }

  ## genind column names are LOCUS.ALLELE, so neither part may contain a ".".
  locus_name <- .names_without_dots(H$locus)
  allele_name <- character(sum(observed))
  locus_of_column <- character(sum(observed))
  for (a in seq_len(k)) {
    rows <- which(observed[, a])
    allele_name[column[rows, a]] <- vapply(H$alleles[rows], `[`, character(1), a)
    locus_of_column[column[rows, a]] <- locus_name[rows]
  }
  allele_name <- gsub(".", "_", allele_name, fixed = TRUE)
  colnames(tab) <- paste0(locus_of_column, ".", allele_name)

  gi <- adegenet::genind(tab, pop = pop, ploidy = 2L, type = "codom")
  adegenet::`other<-`(gi, value = list(
    locus = stats::setNames(H$locus, locus_name),
    chromosome = if (!is.null(H$fields)) H$fields[, "CHROM"],
    position = if (!is.null(H$fields)) H$fields[, "POS"]))
}

#' @rdname as_genind
#' @param drop_multiallelic If `FALSE` (default), a record with more than two
#'   observed alleles stops `as_genlight()`. If `TRUE`, such records are left
#'   out, with a message.
#' @export
as_genlight <- function(H, popmap = NULL, drop_multiallelic = FALSE, verbose = TRUE) {
  .need_adegenet("as_genlight()")
  .check_flag(drop_multiallelic, "drop_multiallelic")
  .check_flag(verbose, "verbose")
  H <- .resolve_H(H, verbose = FALSE)
  pop <- .adegenet_pop(H, popmap, "as_genlight()")
  H <- .drop_empty_records(H, verbose, "as_genlight()")

  multiallelic <- locus_allele_stats(H)$n_observed_alleles > 2L
  if (any(multiallelic)) {
    if (!drop_multiallelic)
      stop(sum(multiallelic), " record(s) have more than 2 observed alleles, which a genlight ",
           "object can't hold (it is strictly biallelic). Set drop_multiallelic = TRUE to leave ",
           "them out, use the SNP VCF, or use as_genind().", call. = FALSE)
    .inform(verbose, sprintf("as_genlight(): leaving out %s of %s multiallelic records (drop_multiallelic = TRUE)",
                             .big(sum(multiallelic)), .big(length(multiallelic))))
    H <- .subset_H(H, !multiallelic)
  }
  if (!nrow(H$A1)) stop("as_genlight(): no biallelic record is left.", call. = FALSE)

  ## Each record's two alleles. A normal SNP is always REF/ALT (alleles 1 and
  ## 2), even when only one of them is seen in these samples, so a SNP fixed
  ## for ALT gets dosage 2 and objects made from different subsets agree. Only
  ## a record where a third declared allele is seen uses its two observed
  ## alleles, in VCF order. A record with one declared allele gets dosage 0.
  observed <- .allele_counts(H$A1, H$A2, max(1L, H$n_alleles)) > 0L
  ref_alt <- H$n_alleles >= 2L &
    rowSums(observed[, -(1:2), drop = FALSE]) == 0L
  first <- ifelse(ref_alt, 1L, max.col(observed, ties.method = "first"))
  second <- vapply(seq_len(nrow(observed)), function(j) {
    if (ref_alt[j]) return(2L)
    seen <- which(observed[j, ])
    if (length(seen) == 2L) seen[2L]
    else if (H$n_alleles[j] >= 2L) setdiff(seq_len(H$n_alleles[j]), first[j])[1L]
    else NA_integer_
  }, integer(1))
  dosage <- (H$A1 == second) + (H$A2 == second)                 # NA where missing
  dosage[is.na(dosage) & !is.na(H$A1)] <- 0L                     # only one allele declared
  dimnames(dosage) <- list(H$locus, H$samples)

  gl <- adegenet::as.genlight(t(dosage))
  if (!is.null(pop)) gl <- adegenet::`pop<-`(gl, value = pop)
  if (!is.null(H$fields)) {
    gl <- adegenet::`chromosome<-`(gl, value = H$fields[, "CHROM"])
    gl <- adegenet::`position<-`(gl, value = suppressWarnings(as.integer(H$fields[, "POS"])))
  }
  first_allele <- vapply(seq_along(first), function(j) H$alleles[[j]][first[j]], character(1))
  second_allele <- vapply(seq_along(second), function(j)
    if (is.na(second[j])) first_allele[j] else H$alleles[[j]][second[j]], character(1))
  ## genlight's alleles() takes single bases only ("a/c"); haplotype alleles
  ## go in other() instead.
  pairs <- paste0(first_allele, "/", second_allele)
  if (all(nchar(first_allele) == 1L & nchar(second_allele) == 1L))
    gl <- adegenet::`alleles<-`(gl, value = pairs)
  else
    gl <- adegenet::`other<-`(gl, value = list(alleles = pairs))
  gl
}
