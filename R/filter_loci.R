###############################################################################
#
#  R/filter_loci.R -- general-purpose locus and genotype filters.
#
#  Every filter takes `H`, the object returned by read_stacks_vcf() (see
#  R/vcf_io.R), and returns the same kind of object with fewer records
#  and/or some genotypes set to missing. So filters CHAIN, in any order:
#
#    H <- read_stacks_vcf("populations.snps.vcf")
#    H <- filter_call_rate(H, min_call = 0.8)
#    H <- filter_maf(H, min_maf = 0.05)
#    H <- filter_low_conf_alt(H, min_alt_reads = 2)
#
#  (or, with R's pipe, H |> filter_call_rate(0.8) |> filter_maf(0.05)), and
#  the result goes straight into any analysis function in place of a path.
#  Each filter also adds a line to H$filter_log, so print(H) lists what was
#  done to the data (see "The filter log" below).
#
#  These filters are separate from the record rules built into
#  diversity_stats() (min_n / complete_case) and het_between_pops()
#  (min_call), which are documented with those functions. The filters here
#  clean a dataset BEFORE analysis or export.
#
#  WHICH SAMPLES. The filters use EVERY sample in H, including samples that
#  are not in the popmap (an outgroup, a failed library), whereas every
#  analysis uses only the popmap's samples. filter_samples() removes the
#  others from H, so that filters and exported files see only the analyzed
#  individuals.
#
#  A NOTE FOR READERS NEW TO R: `NA` is R's code for a missing value, and a
#  "logical" vector is a vector of TRUE/FALSE values. `H$A1` and `H$A2` are
#  matrices with one row per VCF record and one column per sample;
#  `H$A1[j, i]` is individual i's FIRST allele at record j and `H$A2[j, i]`
#  its SECOND allele, stored as small integers (1 = REF, 2 = first ALT, ...),
#  or NA if the individual was not genotyped there.
#
#  Checked by tests/testthat/test-filter-loci.R.
#
###############################################################################

## Not exported. Every filter that removes records ends with this. It keeps
## the records `keep` (a TRUE/FALSE vector with one value per record, or a
## vector of record numbers) in every per-record element of H, so they all
## stay lined up. `samples` is untouched: removing records never removes
## individuals.
.subset_H <- function(H, keep) {
  n_rec <- nrow(H$A1)
  if (is.logical(keep)) {
    if (length(keep) != n_rec || anyNA(keep))
      stop(".subset_H(): `keep` must be a TRUE/FALSE vector with exactly ",
           n_rec, " entries (one per record in H) and no missing values.", call. = FALSE)
  } else {
    keep <- as.integer(keep)
    if (anyNA(keep) || any(keep < 1L | keep > n_rec))
      stop(".subset_H(): `keep` record numbers must all be between 1 and ", n_rec, ".",
           call. = FALSE)
  }
  ## drop = FALSE keeps a matrix a matrix even when one record is left.
  H$A1        <- H$A1[keep, , drop = FALSE]
  H$A2        <- H$A2[keep, , drop = FALSE]
  H$locus     <- H$locus[keep]
  H$locus_raw <- H$locus_raw[keep]
  H$alleles   <- H$alleles[keep]
  H$n_alleles <- H$n_alleles[keep]
  ## `fields` and the depth matrices are absent from a hand-built H.
  for (name in c("fields", "depth", "ad_ref", "ad_alt"))
    if (!is.null(H[[name]])) H[[name]] <- H[[name]][keep, , drop = FALSE]
  H
}

## Not exported. The sample-wise twin of .subset_H(): keeps the samples
## (columns) `keep` in every per-sample element of H. `keep` is either a
## TRUE/FALSE vector with one value per sample (H's column order is kept) or
## a vector of sample names (the samples come out in that order). Records are
## untouched, so `n_records_read` stays as it is.
.subset_samples <- function(H, keep) {
  if (is.character(keep)) {
    if (anyNA(keep) || !all(keep %in% H$samples))
      stop(".subset_samples(): every name in `keep` must be a sample of H.", call. = FALSE)
    keep <- match(keep, H$samples)
  } else if (!is.logical(keep) || length(keep) != length(H$samples) || anyNA(keep)) {
    stop(".subset_samples(): `keep` must be a TRUE/FALSE vector with one entry per ",
         "sample in H, or a vector of sample names.", call. = FALSE)
  }
  H$A1 <- H$A1[, keep, drop = FALSE]
  H$A2 <- H$A2[, keep, drop = FALSE]
  ## The depth matrices are absent from a hand-built H, or from a VCF with no
  ## DP or AD.
  for (name in c("depth", "ad_ref", "ad_alt"))
    if (!is.null(H[[name]])) H[[name]] <- H[[name]][, keep, drop = FALSE]
  H$samples <- H$samples[keep]
  H
}

## ---------------------------------------------------------------------------
## The filter log
##
## Every exported filter adds one row to H$filter_log (except filter_samples()
## and filter_thin_one_snp() when they have nothing to do, which return H
## unchanged), so the data object itself records what was done to it since
## read_stacks_vcf(): print(H) lists the steps (useful for a methods
## section), and diversity_stats() uses the
## log to tell apart RAD loci emptied by an allele-frequency filter (their
## sequenced sites still count) from loci removed as unreliable (their sites
## should not). Internal subsetting (.subset_H() inside kinship_check() or
## write_plink()) is not a filter and is not logged.
## ---------------------------------------------------------------------------

## Not exported. The log of an object nothing has filtered yet: a data frame
## with no rows and these columns:
##   filter           the function, e.g. "filter_mac"
##   setting          its settings as text, e.g. "min_mac = 3"
##   records_removed  records it removed
##   loci_removed     RAD loci it removed entirely (every record of the locus)
##   calls_masked     genotype calls it set to missing
##   samples_removed  samples it removed
.empty_filter_log <- function() {
  data.frame(filter = character(0), setting = character(0), records_removed = integer(0),
             loci_removed = integer(0), calls_masked = integer(0),
             samples_removed = integer(0), stringsAsFactors = FALSE)
}

## Not exported. The counts a filter's log row is measured against, taken
## when the filter starts.
.filter_start <- function(H) {
  list(records = nrow(H$A1), loci = length(unique(H$locus_raw)), samples = length(H$samples))
}

## Not exported. Adds one row to H$filter_log for `filter` (its name) run with
## `setting` (text), comparing H now with `before` (from .filter_start()).
## An object built by hand has no log yet; it gets one here.
.log_filter <- function(H, before, filter, setting, calls_masked = 0L) {
  log <- if (is.null(H$filter_log)) .empty_filter_log() else H$filter_log
  row <- data.frame(filter = filter, setting = setting,
                    records_removed = as.integer(before$records - nrow(H$A1)),
                    loci_removed = as.integer(before$loci - length(unique(H$locus_raw))),
                    calls_masked = as.integer(calls_masked),
                    samples_removed = as.integer(before$samples - length(H$samples)),
                    stringsAsFactors = FALSE)
  H$filter_log <- rbind(log, row)
  H
}

## Not exported. The log as numbered lines of text, one per filter, e.g.
## "1. filter_mac (min_mac = 3): 127 records removed (48 RAD loci emptied)".
## character(0) for an empty or missing log.
.format_filter_log <- function(log) {
  if (is.null(log) || !nrow(log)) return(character(0))
  effect <- function(i) {
    parts <- c(
      if (log$samples_removed[i] > 0) sprintf("%s samples removed", .big(log$samples_removed[i])),
      if (log$records_removed[i] > 0)
        sprintf("%s records removed%s", .big(log$records_removed[i]),
                if (log$loci_removed[i] > 0)
                  sprintf(" (%s RAD loci emptied)", .big(log$loci_removed[i])) else ""),
      if (log$calls_masked[i] > 0) sprintf("%s calls masked", .big(log$calls_masked[i])))
    if (length(parts)) paste(parts, collapse = ", ") else "nothing changed"
  }
  vapply(seq_len(nrow(log)), function(i)
    sprintf("%d. %s (%s): %s", i, log$filter[i], log$setting[i], effect(i)), character(1))
}

#' Keep only the individuals in a popmap
#'
#' Removes from `H` every sample that is not in `popmap`: its genotypes, read
#' depths and allele depths. Records are kept as they are.
#'
#' Every analysis function uses only the popmap's individuals, but the
#' `filter_*()` functions and [locus_allele_stats()] work on EVERY sample in
#' `H`, and the `write_*()` functions write every sample. A sample left out
#' of the popmap on purpose -- an outgroup, a failed library, a replicate --
#' still counts toward minor allele frequencies, call rates and
#' heterozygosity there, and still appears in exported files. Run
#' `filter_samples()` first to filter, and export, only the individuals you
#' analyze.
#'
#' Removing samples can leave records with no genotyped individual, or with
#' only one allele among the remaining individuals; the message counts them.
#' Remove them with [filter_call_rate()] or `filter_mac(H, min_mac = 1)` if
#' needed.
#'
#' @param H The object returned by [read_stacks_vcf()] (optionally filtered).
#' @param popmap A popmap file path or the list returned by [read_popmap()].
#' @param verbose Print which samples were removed. Default `TRUE`.
#' @return `H` with only the samples in `popmap`, in their original order.
#' @examples
#' H <- read_stacks_vcf(system.file("extdata", "small.snps.vcf",
#'                                  package = "RADdiversity"), verbose = FALSE)
#' # A popmap that leaves out one individual:
#' pops <- list(popA = c("popA_1", "popA_2", "popA_3"),
#'              popB = c("popB_1", "popB_2", "popB_3"))
#' H_pop <- filter_samples(H, pops)
#' H_pop$samples
#' # Filters now see only these individuals:
#' H_pop <- filter_mac(H_pop, min_mac = 1)
#' @export
filter_samples <- function(H, popmap, verbose = TRUE) {
  H <- .resolve_H(H, verbose = FALSE)
  .check_flag(verbose, "verbose")
  before <- .filter_start(H)
  pops <- .resolve_pops(popmap, H$samples, verbose = FALSE)
  ids <- unlist(pops, use.names = FALSE)
  keep <- H$samples %in% ids
  if (all(keep)) {
    ## Nothing to do, so H comes back unchanged (and unlogged).
    .inform(verbose, "filter_samples(): every sample in H is in the popmap; nothing removed.")
    return(H)
  }
  removed <- H$samples[!keep]
  H <- .subset_samples(H, keep)
  H <- .log_filter(H, before, "filter_samples", sprintf("popmap with %d populations", length(pops)))
  typed <- rowSums(!is.na(H$A1))
  n_observed <- rowSums(.allele_counts(H$A1, H$A2, k = max(1L, H$n_alleles)) > 0L)
  .inform(verbose, sprintf("filter_samples(): kept %d of %d samples; removed %d not in the popmap: %s%s",
                           sum(keep), length(keep), length(removed),
                           paste(utils::head(removed, 10), collapse = ", "),
                           if (length(removed) > 10) ", ..." else ""))
  if (any(typed == 0L | n_observed < 2L))
    .inform(verbose, sprintf(paste0("  Among the remaining samples, %s of %s records have no ",
                                    "genotype and %s have only one allele; filter_call_rate() ",
                                    "or filter_mac(H, min_mac = 1) removes them."),
                             .big(sum(typed == 0L)), .big(length(typed)),
                             .big(sum(typed > 0L & n_observed < 2L))))
  H
}

## Not exported. The "N of M records kept (x%)" message every filter prints.
.report_kept <- function(verbose, what, kept, n_rec) {
  .inform(verbose, sprintf("%s: %s of %s records kept (%.1f%%)", what, .big(sum(kept)),
                           .big(n_rec), 100 * mean(kept)))
}

#' Per-record allele-frequency statistics
#'
#' For every record in `H`: how many allele copies were genotyped, the
#' frequency of the most common ("major") allele, the minor allele frequency
#' and the minor allele count. This is what [filter_maf()] and [filter_mac()]
#' filter on; look at it before choosing a threshold (for example
#' `hist(locus_allele_stats(H)$maf)`).
#'
#' Every genotyped allele copy is counted by allele. The most common allele is
#' the "major" allele; all others together are the "minor" share. At a
#' biallelic record this is the usual minor allele frequency (MAF); at a
#' multi-allelic (haplotype) record it is the combined frequency of every
#' allele except the most common one.
#'
#' **Which samples.** This uses every sample in `H`, including any that are
#' not in your popmap (an outgroup, say). Run [filter_samples()] first to base
#' it on the individuals you analyze.
#'
#' @param H The object returned by [read_stacks_vcf()] (optionally filtered).
#' @return A data frame with one row per record, in the order of `H`:
#'   \describe{
#'     \item{locus, locus_raw, n_alleles}{Copied from `H`.}
#'     \item{n_observed_alleles}{How many different alleles were actually
#'       seen. Can be smaller than `n_alleles` when the VCF lists an ALT
#'       allele that nobody carries.}
#'     \item{n_called}{Allele copies genotyped (2 per genotyped individual).}
#'     \item{major_af}{Frequency of the most common allele.}
#'     \item{maf}{Minor allele frequency, `1 - major_af`.}
#'     \item{mac}{Minor allele count: copies that are not the major allele.}
#'   }
#'   A record with no genotyped individual gets `NA` (not `0`) for
#'   `major_af`, `maf` and `mac`: "no data" is not "monomorphic".
#' @examples
#' H <- read_stacks_vcf(system.file("extdata", "small.snps.vcf",
#'                                  package = "RADdiversity"), verbose = FALSE)
#' allele_stats <- locus_allele_stats(H)
#' head(allele_stats)
#' summary(allele_stats$maf)
#' @export
locus_allele_stats <- function(H) {
  H <- .resolve_H(H, verbose = FALSE)
  n_rec <- nrow(H$A1)
  ## tab[j, a] = copies of allele a at record j (one tabulate() call; see
  ## .allele_counts() in R/vcf_io.R).
  tab <- .allele_counts(H$A1, H$A2, k = max(1L, H$n_alleles))

  n_called <- rowSums(tab)
  n_observed_alleles <- rowSums(tab > 0L)
  ## The largest count in each row: max.col() finds its column.
  major_count <- if (n_rec) tab[cbind(seq_len(n_rec), max.col(tab, ties.method = "first"))]
                 else integer(0)
  major_af <- major_count / n_called
  mac <- n_called - major_count
  maf <- mac / n_called

  ## A record with no genotyped individual would give 0/0 = NaN here. Mark it
  ## as NA, and its minor allele count as NA rather than 0.
  no_data <- n_called == 0L
  major_af[no_data] <- NA_real_
  maf[no_data] <- NA_real_
  mac[no_data] <- NA_integer_

  data.frame(locus = H$locus, locus_raw = H$locus_raw, n_alleles = H$n_alleles,
             n_observed_alleles = as.integer(n_observed_alleles),
             n_called = as.integer(n_called), major_af = major_af, maf = maf,
             mac = as.integer(mac), row.names = NULL)
}

## Not exported. What the data look filtered at, for data that may have been
## filtered before they reached this package (a Stacks run with --min-mac,
## -R, --max-obs-het; vcftools; dDocent). Computed over the popmap's
## individuals only (`pops`), the same individuals every statistic of
## diversity_stats() uses: a VCF sample left out of the popmap, such as an
## outgroup, could otherwise supply or remove the rare alleles this looks
## for. Returns a list:
##   n_samples         individuals it was computed over
##   n_variable        records with at least 2 observed alleles
##   min_allele_count  over variable records, copies of the rarest observed
##                     allele (for a SNP, the minor allele count)
##   rare_share        share of variable records whose rarest allele has 1 or
##                     2 copies
##   looks_mac_filtered  TRUE when there are at least 200 variable records and
##                     none has an allele with 1 or 2 copies
##   min_call_pooled   lowest call rate of any record, all samples pooled
##   min_call_by_pop   lowest call rate of any record within each population
##   max_ho            highest observed heterozygosity of any record (pooled)
## WHY 200 AND "NONE". Under a neutral site frequency spectrum a site whose
## minor allele is a singleton is the most common kind of SNP: about
## (1 + 1/(N - 1)) / sum_{i<N} 1/i of SNPs among N gene copies, 23% at N = 50.
## Bottlenecks and structure lower that share, but a dataset with hundreds of
## SNPs and not one allele seen once or twice has almost certainly had rare
## alleles removed. Nothing here can detect an HWE filter.
.prior_filter_signals <- function(H, pops) {
  ids <- unlist(pops, use.names = FALSE)
  A1 <- H$A1[, ids, drop = FALSE]
  A2 <- H$A2[, ids, drop = FALSE]
  counts <- .allele_counts(A1, A2)
  observed <- counts > 0L
  n_observed <- rowSums(observed)
  variable <- n_observed >= 2L
  ## Copies of the rarest observed allele of each variable record: a loop over
  ## the (few) allele columns rather than over the (many) records.
  rarest <- rep(.Machine$integer.max, sum(variable))
  for (a in seq_len(ncol(counts))) {
    column <- counts[variable, a]
    rarest <- ifelse(column > 0L, pmin(rarest, column), rarest)
  }
  n_variable <- sum(variable)
  call_rate_pooled <- rowMeans(!is.na(A1))
  typed <- .typed_by_pop(H, pops)
  call_by_pop <- vapply(names(pops), function(p) min(typed[, p]) / length(pops[[p]]), numeric(1))
  ho <- rowMeans(A1 != A2, na.rm = TRUE)
  list(n_samples = length(ids), n_variable = n_variable,
       min_allele_count = if (n_variable) min(rarest) else NA_integer_,
       rare_share = if (n_variable) mean(rarest <= 2L) else NA_real_,
       looks_mac_filtered = n_variable >= 200L && all(rarest > 2L),
       min_call_pooled = min(call_rate_pooled),
       min_call_by_pop = call_by_pop,
       max_ho = suppressWarnings(max(ho, na.rm = TRUE)))
}

## Not exported. Used by filter_maf() and filter_mac(): checks that a
## user-supplied `allele_stats` belongs to this exact H (same records, same
## order), or computes it when not supplied.
.check_or_make_allele_stats <- function(H, allele_stats) {
  if (is.null(allele_stats)) return(locus_allele_stats(H))
  if (!is.data.frame(allele_stats) || nrow(allele_stats) != nrow(H$A1) ||
      !identical(allele_stats$locus, H$locus))
    stop("`allele_stats` does not line up with `H` (different number of records, or ",
         "the records are in a different order). Pass locus_allele_stats(H) computed ",
         "on this exact H, or leave `allele_stats` as NULL to compute it automatically.",
         call. = FALSE)
  allele_stats
}

#' Filter records by minor allele frequency
#'
#' Keeps only records whose minor allele frequency (MAF) is at least
#' `min_maf`. Rare variants are the ones most likely to be sequencing or
#' genotyping errors, so a MAF filter is a common first quality-control step
#' for SNP data.
#'
#' **What it changes.** MAF and MAC filters are routine for RAD-seq data, but
#' rare variants are also real diversity. Removing them lowers He per
#' sequenced site, pct_poly, allelic richness and especially private allelic
#' richness, by an amount that depends on sample size. The threshold also
#' affects inferences of population structure (Linck & Battey 2019). Report
#' the threshold used, and compare diversity values only between datasets
#' filtered the same way. `summary(details = TRUE)` of a [diversity_stats()]
#' result shows the lowest allele count remaining in the data (its short
#' `summary()` says whether the data look filtered by minor allele count).
#'
#' **Which samples.** This uses every sample in `H`, including any that are
#' not in your popmap (an outgroup, say). Run [filter_samples()] first to base
#' it on the individuals you analyze.
#'
#' @references
#' Linck, E. & Battey, C.J. (2019) Minor allele frequency thresholds strongly
#' affect population structure inference with genomic data sets. *Molecular
#' Ecology Resources* 19:639-647.
#'
#' @param H The object returned by [read_stacks_vcf()] (optionally filtered).
#' @param min_maf Minimum minor allele frequency to keep a record, a number
#'   from 0 to 1 (e.g. `0.05` for 5%).
#' @param allele_stats Optional: [locus_allele_stats()] already computed for
#'   this `H`, to save computing it again when you also run [filter_mac()].
#'   Default `NULL`: computed automatically.
#' @param verbose Print how many records were kept. Default `TRUE`.
#' @return `H` with only the records that pass. A record with no genotyped
#'   individual has no MAF and is always dropped.
#' @examples
#' H <- read_stacks_vcf(system.file("extdata", "small.snps.vcf",
#'                                  package = "RADdiversity"), verbose = FALSE)
#' H_common <- filter_maf(H, min_maf = 0.1)
#' H_common
#' @export
filter_maf <- function(H, min_maf, allele_stats = NULL, verbose = TRUE) {
  H <- .resolve_H(H, verbose = FALSE)
  .check_number(min_maf, "min_maf", min = 0, max = 1)
  .check_flag(verbose, "verbose")
  allele_stats <- .check_or_make_allele_stats(H, allele_stats)
  before <- .filter_start(H)
  n_rec <- nrow(allele_stats)
  keep <- !is.na(allele_stats$maf) & allele_stats$maf >= min_maf - .threshold_tol
  .report_kept(verbose, sprintf("MAF filter (minor allele frequency >= %.4g)", min_maf), keep, n_rec)
  if (sum(keep) < 0.05 * n_rec)
    .inform(verbose, "  WARNING: fewer than 5% of records passed this filter. Consider a lower min_maf.")
  if (!any(keep))
    stop("No record has a minor allele frequency >= ", min_maf, ". Lower min_maf, ",
         "or look at locus_allele_stats(H)$maf to see the actual values.", call. = FALSE)
  .log_filter(.subset_H(H, keep), before, "filter_maf", sprintf("min_maf = %g", min_maf))
}

#' Filter records by minor allele count
#'
#' Keeps only records with at least `min_mac` copies of the minor allele,
#' counted over every genotyped individual. The count-based version of
#' [filter_maf()]: a fixed count (e.g. "at least 3 copies") behaves more
#' predictably than a frequency when sample sizes are small or uneven.
#'
#' See [filter_maf()] for what the filter changes. For scale: under a neutral
#' site frequency spectrum, sites with a minor allele count of 2 or less carry
#' about 4/(N - 1) of nucleotide diversity for N gene copies (21% at 10
#' diploids, 10% at 20), so `min_mac = 3` lowers per-site diversity by about
#' that much.
#'
#' **Which samples.** This uses every sample in `H`, including any that are
#' not in your popmap (an outgroup, say). Run [filter_samples()] first to base
#' it on the individuals you analyze.
#'
#' @references
#' Linck, E. & Battey, C.J. (2019) Minor allele frequency thresholds strongly
#' affect population structure inference with genomic data sets. *Molecular
#' Ecology Resources* 19:639-647.
#'
#' @inheritParams filter_maf
#' @param min_mac Minimum minor allele count to keep a record (e.g. `3`).
#' @return `H` with only the records that pass.
#' @examples
#' H <- read_stacks_vcf(system.file("extdata", "small.snps.vcf",
#'                                  package = "RADdiversity"), verbose = FALSE)
#' filter_mac(H, min_mac = 3)
#' @export
filter_mac <- function(H, min_mac, allele_stats = NULL, verbose = TRUE) {
  H <- .resolve_H(H, verbose = FALSE)
  .check_number(min_mac, "min_mac", min = 0)
  .check_flag(verbose, "verbose")
  allele_stats <- .check_or_make_allele_stats(H, allele_stats)
  before <- .filter_start(H)
  n_rec <- nrow(allele_stats)
  keep <- !is.na(allele_stats$mac) & allele_stats$mac >= min_mac
  .report_kept(verbose, sprintf("MAC filter (minor allele count >= %g)", min_mac), keep, n_rec)
  if (sum(keep) < 0.05 * n_rec)
    .inform(verbose, "  WARNING: fewer than 5% of records passed this filter. Consider a lower min_mac.")
  if (!any(keep))
    stop("No record has a minor allele count >= ", min_mac, ". Lower min_mac, ",
         "or look at locus_allele_stats(H)$mac to see the actual values.", call. = FALSE)
  .log_filter(.subset_H(H, keep), before, "filter_mac", sprintf("min_mac = %g", min_mac))
}

#' Filter records by genotyping (call) rate
#'
#' Keeps only records genotyped in enough individuals. Without `popmap`,
#' "enough" is judged over all individuals pooled; with `popmap`, it is
#' judged within each population, because a record can look well covered
#' overall while being poorly covered in one population (see `rule`).
#'
#' This is independent of the call-rate rule inside [het_between_pops()] and
#' the `min_n` / `complete_case` rules of [diversity_stats()].
#'
#' **Which samples.** Without `popmap`, the pooled call rate counts every
#' sample in `H`, including any that are not in your popmap; with `popmap`,
#' only the popmap's samples are used. Run [filter_samples()] first to remove
#' the other samples from `H` altogether.
#'
#' @inheritParams filter_maf
#' @param min_call Minimum fraction of individuals genotyped to keep a record,
#'   a number from 0 to 1 (e.g. `0.8` for 80%).
#' @param popmap Optional: a popmap file path or the list returned by
#'   [read_popmap()], to judge the call rate within each population. Default
#'   `NULL`: pooled over all individuals.
#' @param rule Used only with `popmap`. `"all"` (default) keeps a record only
#'   if EVERY population clears `min_call` -- the right choice before
#'   comparing populations, because otherwise a well-covered population can
#'   carry a poorly covered one. `"any"` keeps a record if AT LEAST ONE
#'   population clears it -- appropriate when each population will be analyzed
#'   on its own.
#' @param verbose Print how many records were kept (and, with `popmap`, a
#'   breakdown by population). Default `TRUE`.
#' @return `H` with only the records that pass.
#' @examples
#' H <- read_stacks_vcf(system.file("extdata", "small.haps.vcf",
#'                                  package = "RADdiversity"), verbose = FALSE)
#' popmap <- system.file("extdata", "small_popmap.tsv", package = "RADdiversity")
#' filter_call_rate(H, min_call = 0.75)                    # pooled
#' filter_call_rate(H, min_call = 0.75, popmap = popmap)   # within each population
#' @export
filter_call_rate <- function(H, min_call, popmap = NULL, rule = "all", verbose = TRUE) {
  H <- .resolve_H(H, verbose = FALSE)
  .check_number(min_call, "min_call", min = 0, max = 1)
  .check_choice(rule, "rule", c("all", "any"))
  .check_flag(verbose, "verbose")
  before <- .filter_start(H)
  n_rec <- nrow(H$A1)

  if (is.null(popmap)) {
    ## rowMeans(!is.na(H$A1)): for each record, the fraction of samples with
    ## a genotype.
    keep <- rowMeans(!is.na(H$A1)) >= min_call - .threshold_tol
    .report_kept(verbose, sprintf("Call-rate filter (pooled across %d samples, >= %.0f%%)",
                                  ncol(H$A1), 100 * min_call), keep, n_rec)
  } else {
    pops <- .resolve_pops(popmap, H$samples, verbose = FALSE)
    passes <- .population_locus_sets(H, pops, min_call)     # records x populations
    ## rowSums(passes) = how many populations clear the threshold.
    keep <- if (rule == "all") rowSums(passes) == length(pops) else rowSums(passes) > 0L
    .report_kept(verbose, sprintf("Call-rate filter (within each population, >= %.0f%%, rule = \"%s\")",
                                  100 * min_call, rule), keep, n_rec)
    for (p in names(pops))
      .inform(verbose, sprintf("    %-22s %s of %s records clear the threshold",
                               p, .big(sum(passes[, p])), .big(n_rec)))
  }
  if (!any(keep))
    stop("No record passed the call-rate filter (min_call = ", min_call, "). Lower min_call, ",
         "or check for a widespread genotyping problem in this dataset.", call. = FALSE)
  setting <- if (is.null(popmap)) sprintf("min_call = %g, pooled", min_call)
             else sprintf("min_call = %g within populations, rule = \"%s\"", min_call, rule)
  .log_filter(.subset_H(H, keep), before, "filter_call_rate", setting)
}

#' Filter out records with excessive heterozygosity
#'
#' Drops records where an unusually large share of individuals are
#' heterozygous. Real heterozygosity has a biological ceiling; a record far
#' above it is the classic sign of an undetected paralog or collapsed repeat,
#' where two similar genomic regions are read as one locus and their
#' differences look like two alleles.
#'
#' This is a crude screen. A threshold on Ho alone also removes genuinely
#' variable loci, and it misses paralogs at low frequency. HDplot (McKinney et
#' al. 2017), which combines heterozygosity with read-ratio deviation, is the
#' more specific tool.
#'
#' **Which samples.** This uses every sample in `H`, including any that are
#' not in your popmap (an outgroup, say). Run [filter_samples()] first to base
#' it on the individuals you analyze.
#'
#' @references
#' McKinney, G.J., Waples, R.K., Seeb, L.W. & Seeb, J.E. (2017) Paralogs are
#' revealed by proportion of heterozygotes and deviations in read ratios in
#' genotyping-by-sequencing data from natural populations. *Molecular Ecology
#' Resources* 17:656-669.
#'
#' @inheritParams filter_maf
#' @param max_ho Maximum observed heterozygosity at a record, a number from 0
#'   to 1 (e.g. `0.5`). Records above it are dropped.
#' @return `H` with only the records that pass. A record with no genotyped
#'   individual is kept ("no data" is not "too heterozygous"); use
#'   [filter_call_rate()] to remove those.
#' @examples
#' H <- read_stacks_vcf(system.file("extdata", "small.haps.vcf",
#'                                  package = "RADdiversity"), verbose = FALSE)
#' filter_max_het(H, max_ho = 0.7)
#' @export
filter_max_het <- function(H, max_ho, verbose = TRUE) {
  H <- .resolve_H(H, verbose = FALSE)
  .check_number(max_ho, "max_ho", min = 0, max = 1)
  .check_flag(verbose, "verbose")
  before <- .filter_start(H)
  ## Observed heterozygosity per record: the fraction of genotyped individuals
  ## whose two alleles differ. With nobody genotyped this is 0/0 = NaN, which
  ## !is.finite() below treats as "no data".
  ho <- rowMeans(H$A1 != H$A2, na.rm = TRUE)
  keep <- !is.finite(ho) | ho <= max_ho + .threshold_tol
  .report_kept(verbose, sprintf("Excess-heterozygosity filter (Ho <= %.3f)", max_ho),
               keep, length(keep))
  if (!any(keep))
    stop("Every record exceeds max_ho = ", max_ho, ". Raise max_ho, ",
         "or check that the genotypes were read correctly.", call. = FALSE)
  .log_filter(.subset_H(H, keep), before, "filter_max_het", sprintf("max_ho = %g", max_ho))
}

#' Thin a SNP dataset to one record per RAD locus
#'
#' Some analyses assume every locus is inherited independently. SNPs on the
#' same RAD tag are physically linked, so treating them as independent
#' overstates how much information the data hold. This keeps exactly one
#' record per RAD locus (`H$locus_raw` group).
#'
#' @inheritParams filter_maf
#' @param method `"first"` (default) keeps the first record of each RAD locus,
#'   the same every time. `"random"` picks one record per RAD locus at random.
#' @param seed Used only when `method = "random"`. Default `NULL`: use R's
#'   current random-number stream (call `set.seed()` first for a reproducible
#'   choice). A number makes the choice reproducible on its own and leaves your
#'   session's random-number stream as it was.
#' @return `H` with one record per RAD locus, in the original order.
#' @examples
#' H <- read_stacks_vcf(system.file("extdata", "small.snps.vcf",
#'                                  package = "RADdiversity"), verbose = FALSE)
#' filter_thin_one_snp(H)
#' @export
filter_thin_one_snp <- function(H, method = "first", seed = NULL, verbose = TRUE) {
  H <- .resolve_H(H, verbose = FALSE)
  .check_choice(method, "method", c("first", "random"))
  .check_seed(seed)
  .check_flag(verbose, "verbose")
  before <- .filter_start(H)
  n_rec <- nrow(H$A1)
  if (length(unique(H$locus_raw)) == n_rec) {
    ## Nothing to do, so H comes back unchanged (and unlogged).
    .inform(verbose, "Every record already belongs to its own RAD locus (locus_raw has no ",
            "repeats) -- there is nothing to thin.")
    return(H)
  }

  ## The record numbers of each RAD locus, in order of first appearance.
  records_of_locus <- split(seq_len(n_rec), factor(H$locus_raw, levels = unique(H$locus_raw)))
  if (method == "first") {
    keep <- vapply(records_of_locus, function(i) i[1L], integer(1))
  } else {
    if (!is.null(seed)) {
      restore_rng <- .save_rng_state()
      on.exit(restore_rng(), add = TRUE)
      set.seed(seed)
    }
    ## i[sample.int(length(i), 1)], not sample(i, 1): given a single number x,
    ## sample(x, 1) draws from 1:x instead of returning x, which would be
    ## wrong for every RAD locus with only one record.
    keep <- vapply(records_of_locus, function(i) i[sample.int(length(i), 1L)], integer(1))
  }
  keep <- sort(keep)

  .inform(verbose, sprintf(
    "One-record-per-locus thinning (method = \"%s\"): %s records -> %s RAD loci (%.1f%% of records kept)",
    method, .big(n_rec), .big(length(keep)), 100 * length(keep) / n_rec))
  .log_filter(.subset_H(H, keep), before, "filter_thin_one_snp", sprintf("method = \"%s\"", method))
}

###############################################################################
#
#  Genotype depth: set calls backed by too few (or too many) reads to missing.
#
#  The rule looks only at a call's read depth, never at the genotype itself,
#  so heterozygous and homozygous calls at the same depth are treated the
#  same. That matters: a filter that removes some genotypes more readily than
#  others changes Ho, He and FIS, not just the amount of missing data.
#
###############################################################################

## Not exported. Read depth of every genotype call of H: the DP subfield, or
## the sum of AD where DP is absent, as read_stacks_vcf() stored it in
## H$depth (see .parse_vcf_chunk() in R/vcf_io.R). A records x samples
## matrix, NA where neither could be read.
.genotype_depth <- function(H) {
  if (is.null(H$depth))
    stop("No genotype call in H has a readable depth (DP or AD): H$depth is empty. ",
         "A Stacks haplotype VCF holds genotypes only; filter populations.snps.vcf ",
         "instead, and read it with read_stacks_vcf().", call. = FALSE)
  H$depth
}

#' Set genotypes with too few or too many reads to missing
#'
#' Masks (sets to missing) every genotype call whose read depth is below
#' `min_dp` or above `max_dp`, whatever the genotype. Depth is the call's
#' `DP` value, or the sum of its `AD` (allele depths) where `DP` is absent.
#'
#' Low depth is where heterozygotes are mistaken for homozygotes: with `d`
#' reads, a true heterozygote shows only one allele with probability
#' `2^(1 - d)` (25% at 3 reads, 3% at 6). Very high depth, far above the
#' typical depth in the dataset, is a sign of a collapsed paralog or repeat.
#'
#' Because the rule ignores the genotype, heterozygous and homozygous calls at
#' the same depth are masked alike, and the message reports the share of each
#' that was masked. Prefer this to [filter_low_conf_alt()] before computing
#' diversity statistics: that function masks only calls carrying an ALT
#' allele, which removes heterozygotes preferentially and biases Ho down and
#' FIS up.
#'
#' Stacks' `populations --min-gt-depth` (Stacks 2.67 and later) applies the
#' same minimum-depth rule when the VCF is written. The haplotype VCF
#' (`populations.haps.vcf`) holds genotypes only, without depths, so this
#' filter needs the SNP VCF.
#'
#' @references
#' Rochette, N.C., Rivera-Colon, A.G. & Catchen, J.M. (2019) Stacks 2:
#' analytical methods for paired-end sequencing improve RADseq-based
#' population genomics. *Molecular Ecology* 28:4737-4754.
#'
#' @inheritParams filter_maf
#' @param H The object returned by [read_stacks_vcf()] (optionally filtered).
#'   The read depths come from `H$depth`, which [read_stacks_vcf()] fills from
#'   each call's `DP` (or the sum of its `AD`); it is `NULL` for a VCF with
#'   neither, such as a Stacks haplotype VCF.
#' @param min_dp Calls with fewer reads than this are masked (e.g. `6`).
#' @param max_dp Calls with more reads than this are masked. Default `Inf`:
#'   no upper limit.
#' @return `H`, with the masked calls set to missing. Calls whose depth cannot
#'   be read are left as they are (the message counts them).
#' @examples
#' vcf_lines <- c(
#'   "##fileformat=VCFv4.2",
#'   "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tind1\tind2\tind3",
#'   "un\t1\tlocus_1\tA\tC\t.\tPASS\t.\tGT:DP:AD\t0/0:3:3,0\t0/1:12:6,6\t1/1:40:0,40"
#' )
#' vcf_file <- tempfile(fileext = ".vcf")
#' writeLines(vcf_lines, vcf_file)
#' H <- read_stacks_vcf(vcf_file, verbose = FALSE)
#' H_depth <- filter_genotype_depth(H, min_dp = 6, max_dp = 30)
#' H_depth$A1          # ind1 (3 reads) and ind3 (40 reads) are now missing
#' @export
filter_genotype_depth <- function(H, min_dp, max_dp = Inf, verbose = TRUE) {
  H <- .resolve_H(H, verbose = FALSE)
  .check_number(min_dp, "min_dp", min = 0)
  .check_number(max_dp, "max_dp", min = 0)
  if (max_dp < min_dp)
    stop("`max_dp` (", max_dp, ") must be at least `min_dp` (", min_dp, ").", call. = FALSE)
  .check_flag(verbose, "verbose")
  ## Taken before anything is masked, as in every other filter. Masking removes
  ## no record, RAD locus or sample, so these counts do not change here; taking
  ## them first keeps the log honest if that ever ceases to be true.
  before <- .filter_start(H)
  depth <- .genotype_depth(H)
  typed <- !is.na(H$A1)
  known <- typed & !is.na(depth)
  if (any(typed) && !any(known))
    stop("No genotype call in H has a readable depth (DP or AD). A Stacks haplotype VCF ",
         "holds genotypes only; filter populations.snps.vcf instead.", call. = FALSE)
  masked <- known & (depth < min_dp | depth > max_dp)
  heterozygous <- typed & H$A1 != H$A2
  share <- function(of) if (any(of)) 100 * sum(masked & of) / sum(of) else NA_real_
  .inform(verbose, sprintf(
    "Genotype-depth filter (%g <= reads <= %g): %s of %s calls masked (%.1f%% of heterozygous, %.1f%% of homozygous calls)",
    min_dp, max_dp, .big(sum(masked)), .big(sum(typed)), share(heterozygous), share(typed & !heterozygous)))
  if (any(typed & !known))
    .inform(verbose, sprintf("  %s calls have no readable depth and were kept.", .big(sum(typed & !known))))
  H$A1[masked] <- NA_integer_
  H$A2[masked] <- NA_integer_
  .log_filter(H, before, "filter_genotype_depth",
              sprintf("min_dp = %g, max_dp = %g", min_dp, max_dp), calls_masked = sum(masked))
}

###############################################################################
#
#  Low-confidence ALT calls: genotype calls whose ALT allele is supported by
#  very few sequencing reads.
#
#  A call such as "0/1" or "1/1" means the genotyper decided the individual
#  carries a non-reference allele. If only 1 or 2 reads showed that allele,
#  the call may be a sequencing error. The functions below find such calls
#  from the AD ("allele depth": reads supporting each allele) field of the
#  original VCF, so they can be set to missing, or used to drop records
#  where too many calls look unreliable.
#
###############################################################################

## Not exported. Shared by filter_low_conf_alt(), low_conf_alt_calls() and
## low_conf_alt_sensitivity(): the allele depths read_stacks_vcf() stored in
## H$ad_ref, H$ad_alt and H$depth (see .parse_vcf_chunk() in R/vcf_io.R).
## Returns one row per genotype call that is fully called and carries at
## least one ALT allele -- every call that COULD be flagged:
##   record, sample   row and column of the call in H
##   GT               the genotype, e.g. "0/1" (VCF allele numbers)
##   AD_ref, AD_alt   reads supporting REF, and the ALT allele(s) of the call
##                    (each distinct ALT allele counted once)
##   DP               read depth (sum of AD when DP is absent)
##   alt_fraction     AD_alt / DP
##   usable           TRUE if AD could be read. The VCF specification allows
##                    trailing FORMAT fields to be dropped, so a call can have
##                    no AD; such a call is never flagged.
.parse_alt_ad <- function(H) {
  ## A VCF with no AD field (such as a Stacks haplotype VCF) has no allele
  ## depths to read: every call is then "not usable" and nothing is flagged.
  ## Only an H built by hand, without read_stacks_vcf()'s columns, is an error.
  if ((is.null(H$ad_ref) || is.null(H$ad_alt)) && is.null(H$fields))
    stop("This needs the AD (allele depth) values that read_stacks_vcf() stores in ",
         "H$ad_ref and H$ad_alt, but this H has neither those nor H$fields -- for ",
         "example because it was built by hand. Read the data with read_stacks_vcf().",
         call. = FALSE)
  A1 <- H$A1
  A2 <- H$A2
  n_rec <- nrow(A1)

  ## A call "has an ALT allele" if both alleles are called and at least one is
  ## above 1 (allele 1 is REF).
  cell <- which(!is.na(A1) & !is.na(A2) & (A1 > 1L | A2 > 1L))
  if (!length(cell))
    return(data.frame(record = integer(), sample = integer(), GT = character(),
                      AD_ref = integer(), AD_alt = integer(), DP = integer(),
                      alt_fraction = double(), usable = logical()))
  ## Matrix cells are numbered down each column; turn the numbers back into
  ## (record, sample) pairs.
  record <- ((cell - 1L) %% n_rec) + 1L
  sample <- ((cell - 1L) %/% n_rec) + 1L

  AD_ref <- if (is.null(H$ad_ref)) rep(NA_integer_, length(cell)) else H$ad_ref[cell]
  AD_alt <- if (is.null(H$ad_alt)) rep(NA_integer_, length(cell)) else H$ad_alt[cell]
  usable <- !is.na(AD_ref)
  DP <- if (is.null(H$depth)) rep(NA_integer_, length(cell)) else H$depth[cell]
  DP[!usable] <- NA_integer_
  gt_text <- paste0(A1[cell] - 1L, "/", A2[cell] - 1L)

  data.frame(record = record, sample = sample, GT = gt_text,
             AD_ref = AD_ref, AD_alt = AD_alt, DP = DP,
             alt_fraction = ifelse(!is.na(DP) & DP > 0L, AD_alt / DP, NA_real_),
             usable = usable)
}

## Not exported. Per record: ALT-containing calls, flagged calls, and the
## flagged fraction, from the .parse_alt_ad() table.
.flag_summary <- function(H, calls, min_alt_reads) {
  n_rec <- nrow(H$A1)
  flagged <- calls$usable & calls$AD_alt <= min_alt_reads
  n_alt_calls <- tabulate(calls$record, nbins = n_rec)
  n_flagged <- tabulate(calls$record[flagged], nbins = n_rec)
  list(flagged = flagged,
       per_record = data.frame(locus = H$locus, n_alt_calls = n_alt_calls, n_flagged = n_flagged,
                               flagged_fraction = ifelse(n_alt_calls > 0L, n_flagged / n_alt_calls, 0),
                               row.names = NULL))
}

#' Remove low-confidence ALT genotype calls
#'
#' A genotype call carrying an ALT allele (e.g. `0/1`, `1/1`) is flagged when
#' the reads supporting that ALT allele number `min_alt_reads` or fewer: the
#' genotyper called an alternate allele on very little evidence, which is the
#' profile of a sequencing error. `mode = "mask"` sets just the flagged calls
#' to missing; `mode = "drop"` removes whole records where too large a share of
#' the ALT-containing calls were flagged. Use [low_conf_alt_calls()] to see
#' which calls are flagged, and [low_conf_alt_sensitivity()] to choose
#' `min_alt_reads`.
#'
#' **Not before diversity statistics.** Only calls that carry an ALT allele
#' can be flagged; a REF homozygote on equally few reads never is. Masking
#' therefore removes heterozygotes more readily than homozygotes, which biases
#' Ho and allele frequencies down and FIS up, most in low-coverage libraries.
#' A true heterozygote read 6 times shows 2 or fewer ALT reads 34% of the
#' time (15% at 8 reads, 5.5% at 10), before the genotype caller's own
#' threshold. Before [diversity_stats()], [het_between_pops()] or
#' [individual_inbreeding()], use [filter_genotype_depth()], which masks by
#' depth whatever the genotype. The message reports the share of heterozygous
#' calls masked, so the imbalance can be seen.
#'
#' **Which samples.** This uses every sample in `H`, including any that are
#' not in your popmap (an outgroup, say). Run [filter_samples()] first to base
#' it on the individuals you analyze.
#'
#' @param H The object returned by [read_stacks_vcf()] (optionally filtered).
#'   The allele depths come from `H$ad_ref` and `H$ad_alt`, which
#'   [read_stacks_vcf()] fills from each call's `AD` (allele depth) field; they
#'   are `NULL` for a VCF without `AD`, in which case nothing is flagged.
#' @param min_alt_reads A call with this many or fewer reads supporting its
#'   ALT allele(s) is flagged. Default `2`.
#' @param mode `"mask"` (default): set only the flagged calls to missing.
#'   `"drop"`: remove whole records instead (see `drop_frac`). To compare both,
#'   run the function twice on the same starting `H`: after masking, the
#'   flagged calls are gone, so a second, dropping pass finds nothing to drop.
#' @param drop_frac Used only with `mode = "drop"`. A record is removed when
#'   the fraction of its ALT-containing calls that were flagged is greater
#'   than this. Default `0`: remove a record if any call is flagged. A record
#'   with no ALT-containing call is never removed.
#' @param verbose Print how many calls were flagged. Default `TRUE`.
#' @return `H`, with flagged calls set to missing (`mode = "mask"`) or
#'   records removed (`mode = "drop"`).
#' @examples
#' # A tiny VCF: ind2's ALT call at locus_1 is supported by only 1 read.
#' vcf_lines <- c(
#'   "##fileformat=VCFv4.2",
#'   "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tind1\tind2",
#'   "un\t1\tlocus_1\tA\tC\t.\tPASS\t.\tGT:AD:DP\t0/0:10,0:10\t0/1:9,1:10"
#' )
#' vcf_file <- tempfile(fileext = ".vcf")
#' writeLines(vcf_lines, vcf_file)
#' H <- read_stacks_vcf(vcf_file, verbose = FALSE)
#' low_conf_alt_calls(H, min_alt_reads = 2)$flagged_calls
#' H_masked <- filter_low_conf_alt(H, min_alt_reads = 2)
#' H_masked$A1          # ind2's genotype is now missing
#' @export
filter_low_conf_alt <- function(H, min_alt_reads = 2, mode = "mask", drop_frac = 0,
                                verbose = TRUE) {
  H <- .resolve_H(H, verbose = FALSE)
  .check_number(min_alt_reads, "min_alt_reads", min = 0)
  .check_choice(mode, "mode", c("mask", "drop"))
  .check_number(drop_frac, "drop_frac", min = 0, max = 1)
  .check_flag(verbose, "verbose")
  before <- .filter_start(H)
  calls <- .parse_alt_ad(H)
  flags <- .flag_summary(H, calls, min_alt_reads)
  n_rec <- nrow(H$A1)
  .inform(verbose, sprintf(
    "Low-confidence ALT filter (AD_alt <= %g): %s of %s ALT-containing calls flagged, ",
    min_alt_reads, .big(sum(flags$flagged)), .big(nrow(calls))),
    sprintf("in %s of %s records", .big(sum(flags$per_record$n_flagged > 0L)), .big(n_rec)))

  if (mode == "mask") {
    ## Masking is one-sided (see @details), so say how unevenly it falls.
    typed <- !is.na(H$A1)
    n_het <- sum(typed & H$A1 != H$A2)
    cell <- cbind(calls$record, calls$sample)
    het_flagged <- sum(flags$flagged & H$A1[cell] != H$A2[cell])
    if (n_het > 0)
      .inform(verbose, sprintf(
        "  %.1f%% of heterozygous calls masked; REF-homozygous calls are never masked, so this lowers Ho (see ?filter_low_conf_alt).",
        100 * het_flagged / n_het))
    ## cbind(record, sample) addresses exactly the flagged cells.
    if (any(flags$flagged)) {
      address <- cbind(calls$record[flags$flagged], calls$sample[flags$flagged])
      H$A1[address] <- NA_integer_
      H$A2[address] <- NA_integer_
    }
    return(.log_filter(H, before, "filter_low_conf_alt",
                       sprintf("min_alt_reads = %g, mode = \"mask\"", min_alt_reads),
                       calls_masked = sum(flags$flagged)))
  }
  keep <- flags$per_record$flagged_fraction <= drop_frac + .threshold_tol
  .inform(verbose, sprintf("  drop_frac = %.3g: %s of %s records kept", drop_frac,
                           .big(sum(keep)), .big(n_rec)))
  if (!any(keep))
    stop("Every record would be dropped at drop_frac = ", drop_frac, ". Raise drop_frac, ",
         "or lower min_alt_reads so fewer calls are flagged.", call. = FALSE)
  .log_filter(.subset_H(H, keep), before, "filter_low_conf_alt",
              sprintf("min_alt_reads = %g, mode = \"drop\", drop_frac = %g", min_alt_reads, drop_frac))
}

#' List low-confidence ALT genotype calls
#'
#' The calls [filter_low_conf_alt()] would flag, and a per-record summary, so
#' they can be inspected before (or instead of) filtering.
#'
#' @inheritParams filter_low_conf_alt
#' @return A list:
#'   \describe{
#'     \item{flagged_calls}{One row per flagged call: `locus`, `sample`, `GT`,
#'       `AD_ref`, `AD_alt`, `DP`, `alt_fraction`.}
#'     \item{locus_summary}{One row per record: `locus`, `n_alt_calls`,
#'       `n_flagged`, `flagged_fraction`.}
#'   }
#' @examples
#' H <- read_stacks_vcf(system.file("extdata", "small.haps.vcf",
#'                                  package = "RADdiversity"), verbose = FALSE)
#' flagged <- low_conf_alt_calls(H, min_alt_reads = 2)
#' head(flagged$locus_summary)
#' @export
low_conf_alt_calls <- function(H, min_alt_reads = 2) {
  H <- .resolve_H(H, verbose = FALSE)
  .check_number(min_alt_reads, "min_alt_reads", min = 0)
  calls <- .parse_alt_ad(H)
  flags <- .flag_summary(H, calls, min_alt_reads)
  f <- flags$flagged
  list(flagged_calls = data.frame(locus = H$locus[calls$record[f]],
                                  sample = H$samples[calls$sample[f]],
                                  GT = calls$GT[f], AD_ref = calls$AD_ref[f],
                                  AD_alt = calls$AD_alt[f], DP = calls$DP[f],
                                  alt_fraction = calls$alt_fraction[f], row.names = NULL),
       locus_summary = flags$per_record)
}

#' Sensitivity table for the low-confidence-ALT threshold
#'
#' For a range of `min_alt_reads` thresholds, how many ALT-containing genotype
#' calls [filter_low_conf_alt()] would flag. Use it to choose a threshold.
#'
#' @inheritParams filter_low_conf_alt
#' @param thresholds The `min_alt_reads` values to try. Default
#'   `c(1, 2, 3, 4, 5, 10)`.
#' @return A data frame with one row per threshold: `threshold`, `n_flagged`,
#'   `pct_flagged` and `total_usable` (ALT-containing calls whose AD could be
#'   read).
#' @examples
#' H <- read_stacks_vcf(system.file("extdata", "small.haps.vcf",
#'                                  package = "RADdiversity"), verbose = FALSE)
#' low_conf_alt_sensitivity(H)
#' @export
low_conf_alt_sensitivity <- function(H, thresholds = c(1, 2, 3, 4, 5, 10)) {
  H <- .resolve_H(H, verbose = FALSE)
  if (!(is.numeric(thresholds) && length(thresholds) && !anyNA(thresholds)))
    stop("`thresholds` must be a vector of numbers.", call. = FALSE)
  calls <- .parse_alt_ad(H)
  usable_ad <- calls$AD_alt[calls$usable]
  total <- length(usable_ad)
  n_flagged <- vapply(thresholds, function(t) sum(usable_ad <= t), integer(1))
  data.frame(threshold = thresholds, n_flagged = n_flagged,
             pct_flagged = if (total) 100 * n_flagged / total else rep(NA_real_, length(thresholds)),
             total_usable = total)
}
