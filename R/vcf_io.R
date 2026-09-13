###############################################################################
#
#  R/vcf_io.R -- reading VCF and popmap files, and the small helpers every
#  analysis uses to turn them into counts.
#
#  THE DATA OBJECT. read_stacks_vcf() returns an object of class "raddiv_vcf"
#  (called `H` throughout the code). It is a list:
#
#    A1, A2     integer matrices, one row per VCF record, one column per
#               sample: the two alleles of each genotype, numbered 1 = REF,
#               2 = first ALT, 3 = second ALT, ...; NA = missing genotype
#    locus      a unique name for each record (used as a marker name)
#    locus_raw  the RAD locus each record belongs to; SNPs on one RAD tag
#               share it. Standard errors resample these groups.
#    alleles    for each record, the allele sequences (REF first)
#    n_alleles  for each record, how many alleles it declares
#    samples    sample names, in column order
#    fields     the VCF's text columns (CHROM ... FORMAT and genotypes),
#               kept for filter_low_conf_alt() and write_vcf()
#    n_records_read  how many records the file had; filters leave it alone,
#               so it records whether records were removed after reading
#
#  Every filter_*() function returns the same object with fewer records or
#  more missing genotypes, and every analysis function accepts it in place of
#  a file path.
#
#  HAPLOTYPE VS SNP VCF. `populations.haps.vcf` has one record per RAD locus:
#  REF is one haplotype (the variable sites of the tag, concatenated) and ALT
#  lists the others, so records are often multi-allelic. Haplotype loci are
#  one per RAD tag, so they are independent by construction, and allelic
#  richness is a real quantity rather than "1 or 2". On the other hand, one
#  sequencing error anywhere in the tag creates a spurious haplotype, and a
#  paralog adds haplotypes rather than just heterozygotes.
#  `populations.snps.vcf` has one record per SNP.
#
#  Requires only base R.
#
###############################################################################

#' Read a Stacks haplotype or SNP VCF
#'
#' Reads `populations.haps.vcf` or `populations.snps.vcf` (gzip-compressed or
#' not) into genotype matrices plus locus and allele information. The result
#' can be passed to every analysis function in place of the file path, and
#' through any `filter_*()` function first.
#'
#' VCFs from other pipelines (ipyrad, dDocent/freebayes, GATK, ...) are read
#' the same way, provided GT is the first FORMAT field; `locus_from` says how
#' their records are grouped into RAD loci.
#'
#' @param path Path to the VCF file (`.vcf` or `.vcf.gz`).
#' @param locus_from How to tell which RAD locus each record belongs to. RAD
#'   loci are the unit that standard errors and bootstrap intervals resample,
#'   so that linked SNPs on one RAD tag move together.
#'   * `"auto"` (default): the ID column when every record has one (Stacks
#'     writes the locus number, or `locus:column:strand` in a SNP VCF), and
#'     `"window"` otherwise.
#'   * `"ID"`: the ID column, up to its first `:`.
#'   * `"CHROM"`: one locus per CHROM value, for de novo formats that put the
#'     locus there (some Stacks 2 versions, ipyrad, dDocent).
#'   * `"window"`: records on the same CHROM within `window_bp` of the
#'     previous record form one locus. This gives one locus per CHROM for de
#'     novo formats and one per RAD locus (or run of adjacent loci) for
#'     reference-aligned data. Choose it explicitly for a reference-aligned
#'     VCF whose IDs are unique per SNP (e.g. dbSNP rs numbers).
#' @param window_bp Largest gap, in base pairs, between consecutive records of
#'   one locus under `locus_from = "window"`. Default `1000`.
#' @param verbose Print a short summary of what was read. Default `TRUE`.
#' @return An object of class `raddiv_vcf`: a list with elements `A1`, `A2`
#'   (allele-number matrices, one row per record, one column per sample;
#'   1 = REF, 2 = first ALT, ...; `NA` = missing), `locus` (unique record
#'   names), `locus_raw` (the RAD locus of each record), `alleles` (allele
#'   sequences per record), `n_alleles`, `samples`, `fields` (the raw VCF
#'   columns) and `n_records_read` (records in the file, unchanged by the
#'   `filter_*()` functions). Printing it shows a short summary.
#' @examples
#' H <- read_stacks_vcf(system.file("extdata", "small.haps.vcf",
#'                                  package = "RADdiversity"))
#' H                  # a short summary
#' dim(H$A1)          # records x samples
#' H$alleles[[1]]     # the haplotype alleles of the first RAD locus
#' H$A1[1:3, 1:4]     # first allele of each genotype (1 = REF, 2 = first ALT, ...)
#' @export
read_stacks_vcf <- function(path, locus_from = "auto", window_bp = 1000, verbose = TRUE) {
  .check_string(path, "path")
  .check_choice(locus_from, "locus_from", c("auto", "ID", "CHROM", "window"))
  .check_number(window_bp, "window_bp", min = 0)
  .check_flag(verbose, "verbose")
  if (!file.exists(path))
    stop("VCF file not found: ", path, "\n  Check the path and try again.", call. = FALSE)

  ## readLines() decompresses .gz files by itself.
  lines <- readLines(path)
  header_line <- grep("^#CHROM", lines)
  if (length(header_line) != 1L)
    stop("Expected exactly one '#CHROM' header line; found ", length(header_line), ".",
         call. = FALSE)
  header <- strsplit(sub("^#", "", lines[header_line]), "\t", fixed = TRUE)[[1]]
  if (length(header) < 10L)
    stop("The #CHROM header line has ", length(header), " columns, but a VCF with ",
         "genotypes needs at least 10 (CHROM ... FORMAT, then one per sample). ",
         "This file has no sample columns.", call. = FALSE)

  ## Everything after the header line, without blank lines. (Indexing with
  ## -seq_len() rather than (header_line + 1):length(lines), which counts
  ## DOWN when the header is the last line.)
  body <- lines[-seq_len(header_line)]
  body <- body[nzchar(body)]
  if (!length(body)) stop("No variant records found.", call. = FALSE)

  ## Every record must have exactly as many tab-separated fields as the
  ## header. rbind() of unequal-length vectors silently RECYCLES the short
  ## ones, which would shift every later column without any warning.
  split_body <- strsplit(body, "\t", fixed = TRUE)
  n_fields <- lengths(split_body)
  bad_len <- which(n_fields != length(header))
  if (length(bad_len))
    stop(sprintf(
      "Malformed VCF: record %d of %d has %d tab-separated field(s), expected %d ",
      bad_len[1], length(body), n_fields[bad_len[1]], length(header)),
      "to match the #CHROM header line. The file is truncated, corrupted, or ",
      "not a VCF. Fix or regenerate the file and re-run.", call. = FALSE)
  fields <- do.call(rbind, split_body)
  colnames(fields) <- header
  samples <- header[-(1:9)]

  .check_gt_first(fields[, "FORMAT"])

  ## Which RAD locus each record belongs to. See `locus_from` above. Getting
  ## this wrong is silent (linked SNPs would be resampled as if independent),
  ## so the rule used is always reported.
  rule <- .choose_locus_rule(locus_from, fields[, "ID"])
  locus_raw <- .locus_labels(rule, fields[, "CHROM"], fields[, "POS"], fields[, "ID"], window_bp)
  ## `locus` is unique per record and only names markers. Anything that
  ## groups records must use `locus_raw`.
  has_ids <- all(nzchar(fields[, "ID"]) & fields[, "ID"] != ".")
  locus <- if (has_ids) sub(":.*$", "", fields[, "ID"])
           else paste0(fields[, "CHROM"], "_", fields[, "POS"])
  locus <- make.unique(locus)
  n_loci <- length(unique(locus_raw))
  .inform(verbose, sprintf("  RAD loci from %s: %s records on %s loci (%.2f per locus)",
                           .describe_locus_rule(rule, window_bp), .big(nrow(fields)),
                           .big(n_loci), nrow(fields) / n_loci))

  ## Alleles per record: REF, then the comma-separated ALT list ("." = none).
  alt <- fields[, "ALT"]
  alleles <- Map(function(ref, alt_text) {
    if (alt_text == "." || !nzchar(alt_text)) ref
    else c(ref, strsplit(alt_text, ",", fixed = TRUE)[[1]])
  }, fields[, "REF"], alt, USE.NAMES = FALSE)
  n_alleles <- lengths(alleles)

  ## Genotypes -> two allele-number matrices.
  gt <- .parse_gt(fields[, -(1:9), drop = FALSE])
  dimnames(gt$A1) <- dimnames(gt$A2) <- list(NULL, samples)

  ## An allele number above the number of declared alleles means the record
  ## and its genotypes disagree: stop rather than compute nonsense. (Matrix
  ## rows are records, so comparing with the per-record `n_alleles` vector
  ## lines up row by row.)
  over <- which((!is.na(gt$A1) & gt$A1 > n_alleles) | (!is.na(gt$A2) & gt$A2 > n_alleles))
  if (length(over))
    stop(length(over), " genotype calls reference an allele number beyond REF+ALT. ",
         "The VCF is malformed.", call. = FALSE)

  if (verbose) {
    message(sprintf("  %s records, %d samples", .big(nrow(fields)), length(samples)))
    message(sprintf("  alleles per record: median %d, max %d; %.1f%% of records have more than 2",
                    as.integer(stats::median(n_alleles)), max(n_alleles),
                    100 * mean(n_alleles > 2)))
    message(sprintf("  missing genotype rate: %.2f%%", 100 * mean(is.na(gt$A1))))
    if (length(gt$unrecognised))
      message("  treated as missing: ", paste(utils::head(gt$unrecognised, 8), collapse = ", "))
  }

  ## `n_records_read` is never changed by the filter_*() functions, so an
  ## analysis can tell that records were removed after reading (see the
  ## per-site check in diversity_stats()).
  structure(list(A1 = gt$A1, A2 = gt$A2, locus = locus, locus_raw = locus_raw,
                 alleles = alleles, n_alleles = n_alleles, samples = samples,
                 fields = fields, n_records_read = nrow(fields)),
            class = "raddiv_vcf")
}

#' @rdname read_stacks_vcf
#' @param x A `raddiv_vcf` object, as returned by `read_stacks_vcf()`.
#' @param ... Ignored.
#' @export
print.raddiv_vcf <- function(x, ...) {
  n_rec <- nrow(x$A1)
  cat("<raddiv_vcf> genotype data read by read_stacks_vcf()\n")
  cat(sprintf("  %s records on %s RAD loci; %d samples\n", .big(n_rec),
              .big(length(unique(x$locus_raw))), ncol(x$A1)))
  if (n_rec > 0) {
    cat(sprintf("  record type: %s\n", if (.is_haplotype_H(x))
      "haplotype (one multi-allelic record per RAD locus)" else "SNP (one site per record)"))
    cat(sprintf("  missing genotypes: %.2f%%\n", 100 * mean(is.na(x$A1))))
  }
  cat("  samples:", paste(utils::head(x$samples, 6), collapse = ", "),
      if (length(x$samples) > 6) "..." else "", "\n")
  cat("  elements: $A1 $A2 (allele matrices), $locus, $locus_raw, $alleles,",
      "$n_alleles, $samples, $fields\n")
  invisible(x)
}

## ---------------------------------------------------------------------------
## VCF parsing helpers (shared with pi_allsites(), which reads in chunks)
## ---------------------------------------------------------------------------

## Not exported. Stops unless every FORMAT value starts with GT. The VCF
## specification puts GT first whenever it is present, and the genotype
## parser relies on it; a FORMAT such as "DP:GT" would otherwise be read as
## plausible-looking but wrong allele numbers.
.check_gt_first <- function(format_values, first_record = 1L) {
  bad <- which(!grepl("^GT(:|$)", format_values))
  if (length(bad))
    stop(sprintf("Malformed VCF: FORMAT field at record %d is \"%s\", not GT-first. ",
                 first_record - 1L + bad[1], format_values[bad[1]]),
         "This package requires GT to be the first FORMAT subfield, as the VCF ",
         "specification and Stacks output have it.", call. = FALSE)
  invisible(NULL)
}

## Not exported. Genotype text ("0/1:12:5,7", "1|1", "./.", ...) -> two
## integer allele-number matrices (VCF allele 0 becomes 1, so the numbers
## index `alleles` directly). Only a clean diploid call "<number>/<number>"
## or "<number>|<number>" is kept; anything else (./., ., a haploid call, a
## half-missing call) becomes NA in both matrices. `unrecognised` lists the
## distinct rejected values other than the usual missing codes.
.parse_gt <- function(genotype_text) {
  gt <- sub(":.*$", "", genotype_text)          # GT is the first subfield
  ok <- grepl("^[0-9]+[/|][0-9]+$", gt)
  A1 <- A2 <- matrix(NA_integer_, nrow(genotype_text), ncol(genotype_text))
  A1[ok] <- as.integer(sub("[/|].*$", "", gt[ok])) + 1L
  A2[ok] <- as.integer(sub("^.*[/|]", "", gt[ok])) + 1L
  unrecognised <- setdiff(unique(gt[!ok]), c("./.", ".|.", ".", "", "./", "/."))
  list(A1 = A1, A2 = A2, unrecognised = unrecognised)
}

## Not exported. Resolves `locus_from = "auto"` to a concrete rule, and stops
## if "ID" was requested but some records have no ID.
.choose_locus_rule <- function(locus_from, ids) {
  has_ids <- all(nzchar(ids) & ids != ".")
  rule <- if (locus_from == "auto") (if (has_ids) "ID" else "window") else locus_from
  if (rule == "ID" && !has_ids)
    stop("locus_from = \"ID\", but some records have no ID (\".\"). Use ",
         "locus_from = \"CHROM\" or \"window\" (see ?read_stacks_vcf).", call. = FALSE)
  rule
}

## Not exported. The RAD locus label of each record under `rule`.
.locus_labels <- function(rule, chrom, pos, ids, window_bp) {
  switch(rule,
    ID     = sub(":.*$", "", ids),
    CHROM  = chrom,
    window = .window_blocks(chrom, pos, window_bp))
}

## Not exported. How a locus rule is described in progress messages.
.describe_locus_rule <- function(rule, window_bp) {
  switch(rule,
    ID     = "the ID column",
    CHROM  = "CHROM",
    window = sprintf("CHROM + POS (records within %s bp)", .big(window_bp)))
}

## Not exported. Groups records into loci by position: sorted within each
## CHROM, a record starts a new locus when it is on a different CHROM from the
## previous record or more than `window_bp` beyond it. Each locus is labelled
## "CHROM:POS" of its first record; labels come back in the input order.
.window_blocks <- function(chrom, pos, window_bp) {
  pos <- suppressWarnings(as.numeric(pos))
  if (anyNA(pos))
    stop("Some POS values are not numbers, so records cannot be grouped by ",
         "position. Use locus_from = \"ID\" or \"CHROM\".", call. = FALSE)
  o <- order(chrom, pos)
  chrom_sorted <- chrom[o]
  pos_sorted <- pos[o]
  starts_locus <- c(TRUE, chrom_sorted[-1] != chrom_sorted[-length(chrom_sorted)] |
                          diff(pos_sorted) > window_bp)
  first_label <- paste0(chrom_sorted, ":",
                        format(pos_sorted, scientific = FALSE, trim = TRUE))[starts_locus]
  out <- character(length(chrom))
  out[o] <- first_label[cumsum(starts_locus)]
  out
}

## ---------------------------------------------------------------------------
## Popmaps
## ---------------------------------------------------------------------------

#' Read a Stacks-style popmap
#'
#' Reads a two-column, no-header, tab-separated popmap (`sample_id <TAB>
#' population`) into a list of sample IDs per population. Every function that
#' takes a `popmap` argument accepts either the file path or this list.
#'
#' @param path Path to the popmap file.
#' @param samples Optional character vector of sample IDs to keep (typically
#'   `H$samples`, the VCF's sample columns). Samples not in it are dropped,
#'   and so is a population left with no samples. Default `NULL`: keep every
#'   sample in the file.
#' @param verbose Print a short summary. Default `TRUE`.
#' @return A named list of character vectors, one per population (in the
#'   order populations first appear in the file), each holding that
#'   population's sample IDs.
#' @examples
#' popmap <- system.file("extdata", "small_popmap.tsv", package = "RADdiversity")
#' pops <- read_popmap(popmap)
#' pops
#' lengths(pops)    # individuals per population
#' @export
read_popmap <- function(path, samples = NULL, verbose = TRUE) {
  .check_string(path, "path")
  .check_flag(verbose, "verbose")
  if (!file.exists(path))
    stop("Popmap file not found: ", path, "\n  Check the path and try again.", call. = FALSE)
  ## colClasses = "character": otherwise read.delim() guesses column types,
  ## so an all-numeric sample ID such as "001" becomes the number 1 (and no
  ## longer matches the VCF), and a population named "T" becomes TRUE.
  pm <- utils::read.delim(path, header = FALSE, colClasses = "character",
                          col.names = c("sample", "pop"), strip.white = TRUE,
                          comment.char = "#", quote = "")
  ## Keep only samples in the VCF BEFORE grouping, so populations come out in
  ## the order they first appear among the retained samples.
  if (!is.null(samples)) pm <- pm[pm$sample %in% samples, , drop = FALSE]
  if (!nrow(pm))
    stop(if (is.null(samples)) "The popmap file has no samples."
         else "No popmap sample names match the VCF.", call. = FALSE)
  pops <- split(pm$sample, factor(pm$pop, levels = unique(pm$pop)))
  .clean_pops(pops, samples, verbose)
}

## Not exported. Accepts `popmap` as either a path to a popmap file or a list
## already returned by read_popmap() (or written by hand), and returns the
## list restricted to `samples`. This is what lets every analysis function
## take either form.
.resolve_pops <- function(popmap, samples, verbose = TRUE) {
  if (is.character(popmap) && length(popmap) == 1L)
    return(read_popmap(popmap, samples, verbose = verbose))
  .clean_pops(popmap, samples, verbose)
}

## Not exported. Checks a population list and restricts it to `samples`
## (when given): drops unknown samples and emptied populations, and stops on
## a sample assigned to two populations.
.clean_pops <- function(pops, samples = NULL, verbose = TRUE) {
  if (!is.list(pops) || is.null(names(pops)) || any(!nzchar(names(pops))) ||
      !all(vapply(pops, is.character, logical(1))))
    stop("`popmap` must be the path to a popmap file, or a named list of sample ",
         "IDs per population as returned by read_popmap() (got an object of class ",
         paste(class(pops), collapse = "/"), ").", call. = FALSE)
  if (anyDuplicated(names(pops)))
    stop("Population names must be unique; repeated: ",
         paste(unique(names(pops)[duplicated(names(pops))]), collapse = ", "), call. = FALSE)
  if (!is.null(samples)) {
    pops <- lapply(pops, function(ids) ids[ids %in% samples])
    pops <- pops[lengths(pops) > 0]
    if (!length(pops)) stop("No popmap sample names match the VCF.", call. = FALSE)
  }
  ## A sample listed twice would be counted in both populations and silently
  ## corrupt every between-population statistic.
  ids <- unlist(pops, use.names = FALSE)
  repeated <- unique(ids[duplicated(ids)])
  if (length(repeated))
    stop(length(repeated), " sample(s) appear more than once in the popmap: ",
         paste(utils::head(repeated, 10), collapse = ", "),
         "\n  Each individual must be assigned to exactly one population.", call. = FALSE)
  if (!is.null(samples)) {
    n_absent <- length(setdiff(samples, ids))
    if (n_absent)
      .inform(verbose, "  ", n_absent, " VCF samples absent from the popmap; excluded.")
  }
  .inform(verbose, "  populations: ",
          paste(sprintf("%s (n=%d)", names(pops), lengths(pops)), collapse = ", "))
  pops
}

## ---------------------------------------------------------------------------
## Accepting a path or an already-read object
## ---------------------------------------------------------------------------

## Not exported. Every analysis function accepts `vcf` as EITHER a path to a
## VCF file (read here with read_stacks_vcf()) OR the object read_stacks_vcf()
## returns, optionally passed through filter_*() functions. A list is checked
## for the elements the analyses need, so that a hand-built or damaged object
## fails here with a clear message rather than deep inside a calculation.
.resolve_H <- function(vcf, verbose = TRUE) {
  if (is.character(vcf) && length(vcf) == 1L) {
    .inform(verbose, "Reading ", vcf, " ...")
    return(read_stacks_vcf(vcf, verbose = verbose))
  }
  if (!is.list(vcf))
    stop("`vcf` must be either the path to a VCF file, or the object returned by ",
         "read_stacks_vcf() (optionally passed through filter_*() functions). Got an ",
         "object of class: ", paste(class(vcf), collapse = "/"), call. = FALSE)
  required <- c("A1", "A2", "locus", "locus_raw", "alleles", "n_alleles", "samples")
  missing_elements <- setdiff(required, names(vcf))
  if (length(missing_elements))
    stop("`vcf` is a list, but it is missing element(s) that read_stacks_vcf() ",
         "always includes: ", paste(missing_elements, collapse = ", "),
         ". Pass the object returned by read_stacks_vcf() (optionally filtered).",
         call. = FALSE)
  if (!is.matrix(vcf$A1) || !is.matrix(vcf$A2) || !identical(dim(vcf$A1), dim(vcf$A2)))
    stop("`vcf$A1` and `vcf$A2` must both be matrices of the same size ",
         "(as read_stacks_vcf() always produces).", call. = FALSE)
  vcf
}

## Not exported. Checks shared by the functions that take `vcf` and `popmap`,
## run before any (possibly slow) reading.
.check_run_inputs <- function(vcf, popmap, stem, outdir) {
  if (is.character(vcf) && length(vcf) == 1L && !file.exists(vcf))
    stop("VCF file not found: ", vcf, "\n  Check the path and try again.", call. = FALSE)
  if (is.character(popmap) && length(popmap) == 1L && !file.exists(popmap))
    stop("Popmap file not found: ", popmap, "\n  Check the path and try again.", call. = FALSE)
  if (!is.null(stem)) .check_string(stem, "stem")
  if (!is.null(outdir)) .check_string(outdir, "outdir")
  ## `stem` names output files. A path supplies one (see .derive_stem()); an
  ## already-read object does not, so it needs `stem` when files are written.
  if (!is.character(vcf) && is.null(stem) && !is.null(outdir))
    stop("`vcf` is an already-read object rather than a file path, so its file ",
         "name can't be used to name the output files. Pass `stem` as well, e.g. ",
         "stem = \"haps\" or stem = \"snps\", matching which VCF the data came from.",
         call. = FALSE)
  invisible(NULL)
}

## Not exported. The text used in output file names, e.g.
## diversity_per_population.<stem>.tsv. A given `stem` always wins; otherwise
## it comes from the VCF's file name (populations.haps.vcf.gz -> "haps"), so
## runs on the SNP and haplotype VCFs can share an output directory without
## overwriting each other.
.derive_stem <- function(vcf, stem = NULL) {
  if (!is.null(stem)) return(stem)
  s <- sub("\\.gz$", "", basename(vcf))
  s <- sub("\\.vcf$", "", s)
  s <- sub("^.*\\.", "", s)                            # populations.haps -> haps
  if (!nzchar(s) || grepl("[^A-Za-z0-9_-]", s)) "out" else s
}

## ---------------------------------------------------------------------------
## Counting helpers used by the analyses
## ---------------------------------------------------------------------------

## Not exported. Is H a haplotype VCF (one multi-allelic record per RAD tag)
## rather than a SNP VCF (one site per record)? Either signal is decisive:
## Stacks writes multi-nucleotide alleles ("AC", "CA") for haplotype records
## and single bases for SNP records, and haplotype records are often
## multi-allelic.
.is_haplotype_H <- function(H) {
  max(H$n_alleles) > 2L || any(nchar(unlist(H$alleles, use.names = FALSE)) > 1L)
}

## Not exported. Typed (non-missing) individuals per record and population,
## for the records `rows`: an integer matrix, one row per record, one column
## per population.
.typed_by_pop <- function(H, pops, rows = seq_len(nrow(H$A1))) {
  out <- vapply(pops, function(ids) rowSums(!is.na(H$A1[rows, ids, drop = FALSE])),
                numeric(length(rows)))
  out <- matrix(out, nrow = length(rows), dimnames = list(NULL, names(pops)))
  storage.mode(out) <- "integer"
  out
}

## Not exported. records x populations TRUE/FALSE matrix: does population p
## ITSELF genotype at least `min_call` of its individuals at the record?
## Decided for each population independently of every other population. Used
## by het_between_pops(), individual_inbreeding(), identity_disequilibrium()
## and filter_call_rate().
.population_locus_sets <- function(H, pops, min_call) {
  call_rate <- sweep(.typed_by_pop(H, pops), 2L, lengths(pops), "/")
  call_rate >= min_call - .threshold_tol
}

## Not exported. Allele counts per record: a records x alleles integer matrix
## built from the records x individuals allele matrices `a1`/`a2` (allele
## numbers 1, 2, ...; NA = missing). Column a holds the number of copies of
## allele a. `k` is the number of allele columns (default: the largest allele
## number present).
##
## How: every allele copy gets a "bin" number that encodes both its record
## (row) and its allele (column), and one tabulate() call counts all bins at
## once -- much faster than looping over records.
.allele_counts <- function(a1, a2, k = max(c(1L, a1, a2), na.rm = TRUE)) {
  n_rec <- nrow(a1)
  codes <- c(a1, a2)
  record_of_code <- rep(rep(seq_len(n_rec), ncol(a1)), 2L)
  typed <- !is.na(codes)
  bins <- (codes[typed] - 1L) * n_rec + record_of_code[typed]
  matrix(tabulate(bins, nbins = n_rec * k), n_rec, k)
}

## Not exported. One allele-count matrix per population for the records
## `rows` of H: a named list of records x `k` integer matrices (see
## .allele_counts()). All matrices have the same `k` columns, so alleles
## absent from a record simply have count 0.
.pop_counts <- function(H, pops, rows = seq_len(nrow(H$A1)),
                        k = max(1L, H$n_alleles[rows])) {
  lapply(pops, function(ids)
    .allele_counts(H$A1[rows, ids, drop = FALSE], H$A2[rows, ids, drop = FALSE], k))
}

## Not exported. H (records `rows`) as the data frame hierfstat expects: a
## population number, then one column per locus with each genotype as a
## 3-digits-per-allele integer (alleles 1 and 12 -> 1012). Used only for the
## optional hierfstat cross-checks and Weir & Goudet's beta.
.to_hierfstat_df <- function(H, pops, rows = seq_len(nrow(H$A1))) {
  ids <- unlist(pops, use.names = FALSE)
  a <- H$A1[rows, ids, drop = FALSE]
  b <- H$A2[rows, ids, drop = FALSE]
  genotype_codes <- t(pmin(a, b) * 1000L + pmax(a, b))      # individuals x loci
  dat <- data.frame(pop = rep(seq_along(pops), lengths(pops)), genotype_codes,
                    row.names = NULL)
  names(dat)[-1] <- paste0("L", seq_along(rows))
  dat
}
