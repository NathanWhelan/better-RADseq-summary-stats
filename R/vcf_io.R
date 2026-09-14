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
#    fields     the VCF's nine fixed text columns (CHROM ... FORMAT), kept for
#               write_vcf() and write_plink()
#    depth      integer matrix like A1: read depth of each call (DP, or the
#               sum of AD); NULL when the VCF has neither
#    ad_ref, ad_alt  integer matrices like A1: reads for REF and for the
#               call's own ALT allele(s), from AD; NULL without AD. Used by
#               filter_genotype_depth() and filter_low_conf_alt(), so the
#               genotype text itself is never kept (it dominated memory)
#    n_records_read, n_loci_read  how many records and RAD loci the file had;
#               filters leave them alone, so they show whether records or whole
#               loci were removed after reading
#    filter_log one row per filter_*() call since reading: which filter, its
#               settings, and the records, RAD loci, calls and samples it
#               removed (see "The filter log" in R/filter_loci.R)
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
#  WHAT STACKS 2.68 ACTUALLY WRITES (export_formats.cc: VcfExport::write_site
#  and VcfHapsExport::write_batch):
#
#    file, assembly               CHROM       POS          ID
#    snps.vcf, de novo            locus ID    column       locus:column
#    snps.vcf, reference-aligned  chromosome  bp position  locus:column:strand
#    haps.vcf, de novo            locus ID    0            .
#    haps.vcf, reference-aligned  chromosome  locus bp     locus:1:strand
#
#  A haplotype VCF also has INFO snp_columns=... (the SNP columns joined into
#  each haplotype), writes a genotype as ./. when any of its SNPs is uncalled
#  (an N in the haplotype), and includes only loci with at least two
#  haplotypes. Because its alleles are the bases at every SNP of the tag, a
#  tag with one SNP has single-base alleles, and a tag with more has
#  multi-base ones. With locus_from = "auto", a de novo haplotype VCF (ID
#  ".") is grouped by the "window" rule; each record has its own CHROM (the
#  Stacks locus ID), so every record is its own RAD locus -- which is
#  correct -- and the progress message says "CHROM + POS".
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
#' **Assumptions to check on data from other pipelines.**
#' * Diploid calls only. Haploid (`0`) and polyploid (`0/0/1`) calls are read
#'   as missing, and the message lists them.
#' * The SNP-VCF analyses were designed for biallelic SNPs, which is what
#'   Stacks writes. Multi-allelic SNPs (e.g. `A` -> `C,T`, which GATK,
#'   freebayes and ipyrad can emit) are counted correctly by Ho, He, FIS and
#'   rarefaction, and do not change how the file is treated; but
#'   [kinship_check()] and [write_plink()] leave them out, and the report's
#'   note that allelic richness is "capped at 2" on SNPs assumes biallelic
#'   sites.
#' * An object is treated as haplotype data (one record spanning several
#'   sites) when any allele is longer than one base; otherwise every record is
#'   taken to be one site. A record with multi-base alleles (an indel or MNP)
#'   therefore makes a SNP VCF look like a haplotype VCF, which switches off
#'   the per-site conversion in [diversity_stats()]. Remove indels and MNPs
#'   from a SNP VCF first.
#' * Symbolic ALT alleles (`*` for a spanning deletion, or `<*>`, `<NON_REF>`,
#'   `<DEL>`, as GATK and bcftools can write) are not DNA sequences. They do
#'   not make a file count as haplotype data, but a genotype that uses one is
#'   still read as an allele. The summary message counts such records; remove
#'   them first (for example with bcftools view).
#' * With `locus_from = "auto"`, a VCF whose ID column is unique per SNP
#'   (e.g. `rs` numbers) makes every SNP its own locus, so linked SNPs are
#'   resampled as if independent. The message says how many records share each
#'   locus; if it is 1.00 for RAD data, use `locus_from = "window"` or
#'   `"CHROM"`.
#'
#' @section What Stacks writes:
#' The columns that identify a record differ between Stacks' two VCFs and
#' between de novo and reference-aligned runs (Stacks 2.68,
#' `export_formats.cc`):
#'
#' | file | CHROM | POS | ID |
#' |---|---|---|---|
#' | `populations.snps.vcf`, de novo | locus ID | column in the locus | `locus:column` |
#' | `populations.snps.vcf`, reference-aligned | chromosome | position | `locus:column:strand` |
#' | `populations.haps.vcf`, de novo | locus ID | `0` | `.` |
#' | `populations.haps.vcf`, reference-aligned | chromosome | locus position | `locus:1:strand` |
#'
#' * In a SNP VCF every record has an ID, so `locus_from = "auto"` takes the
#'   locus from the ID (the part before the first `:`), and SNPs of one RAD
#'   locus are grouped together.
#' * In a de novo haplotype VCF the ID is `.`, so `"auto"` falls back to
#'   `"window"`. Each record has its own CHROM (its locus ID), so each record
#'   becomes its own RAD locus, which is correct: the progress message then
#'   reads "RAD loci from CHROM + POS ... (1.00 per locus)". In a
#'   reference-aligned haplotype VCF the ID gives the locus.
#' * Each haplotype allele is the bases at the locus's SNP columns, joined
#'   (listed in the INFO field as `snp_columns=`). A locus with one SNP
#'   therefore has single-base alleles and a locus with several SNPs has
#'   multi-base alleles; the file is treated as haplotype data when any
#'   allele has more than one base.
#' * A haplotype genotype is written as missing (`./.`) when any of its SNPs
#'   was not called in that individual, and only loci with at least two
#'   haplotypes are written.
#' * The haplotype VCF has genotypes only (no `DP` or `AD`), so
#'   [filter_genotype_depth()] and [filter_low_conf_alt()] need the SNP VCF.
#'
#' @references
#' Danecek, P., Auton, A., Abecasis, G., et al. (2011) The variant call format
#' and VCFtools. *Bioinformatics* 27:2156-2158.
#'
#' Catchen, J., Hohenlohe, P.A., Bassham, S., Amores, A. & Cresko, W.A. (2013)
#' Stacks: an analysis tool set for population genomics. *Molecular Ecology*
#' 22:3124-3140.
#'
#' Rochette, N.C., Rivera-Colon, A.G. & Catchen, J.M. (2019) Stacks 2:
#' analytical methods for paired-end sequencing improve RADseq-based
#' population genomics. *Molecular Ecology* 28:4737-4754.
#'
#' @param path Path to the VCF file (`.vcf` or `.vcf.gz`).
#' @param locus_from How to tell which RAD locus each record belongs to. RAD
#'   loci are the unit that standard errors and bootstrap intervals resample,
#'   so that linked SNPs on one RAD tag move together.
#'   * `"auto"` (default): the ID column when every record has one, and
#'     `"window"` otherwise. For Stacks output both give one locus per Stacks
#'     locus; see "What Stacks writes" below.
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
#' @param chunk_lines Records read and parsed at a time. Default `20000`.
#'   Lower it if reading a very large VCF runs out of memory; the result does
#'   not depend on it.
#' @param verbose Print a short summary of what was read. Default `TRUE`.
#' @return An object of class `raddiv_vcf`: a list with elements `A1`, `A2`
#'   (allele-number matrices, one row per record, one column per sample;
#'   1 = REF, 2 = first ALT, ...; `NA` = missing), `locus` (unique record
#'   names), `locus_raw` (the RAD locus of each record), `alleles` (allele
#'   sequences per record), `n_alleles`, `samples`, `fields` (the VCF's nine
#'   fixed columns, CHROM ... FORMAT), `depth` (read depth of each call, from
#'   DP or the sum of AD; `NULL` if the VCF has neither), `ad_ref` and `ad_alt`
#'   (reads supporting REF and the call's ALT allele(s), from AD; `NULL` if
#'   the VCF has no AD), `n_records_read` and `n_loci_read` (records and RAD
#'   loci in the file, unchanged by the `filter_*()` functions), and
#'   `filter_log` (a data frame with one row per `filter_*()` call made on the
#'   object: `filter`, `setting`, `records_removed`, `loci_removed`,
#'   `calls_masked`, `samples_removed`; no rows straight after reading).
#'   Printing it shows a short summary, including the filters applied.
#'
#'   **Memory.** The genotype text is not kept, only integer matrices: about
#'   20 bytes per genotype call with depths, or 8 without. For example,
#'   100,000 SNPs x 200 samples is about 400 MB.
#' @examples
#' H <- read_stacks_vcf(system.file("extdata", "small.haps.vcf",
#'                                  package = "RADdiversity"))
#' H                  # a short summary
#' dim(H$A1)          # records x samples
#' H$alleles[[1]]     # the haplotype alleles of the first RAD locus
#' H$A1[1:3, 1:4]     # first allele of each genotype (1 = REF, 2 = first ALT, ...)
#' @export
read_stacks_vcf <- function(path, locus_from = "auto", window_bp = 1000, chunk_lines = 20000L,
                            verbose = TRUE) {
  .check_string(path, "path")
  .check_choice(locus_from, "locus_from", c("auto", "ID", "CHROM", "window"))
  .check_number(window_bp, "window_bp", min = 0)
  .check_flag(verbose, "verbose")
  if (!file.exists(path))
    stop("VCF file not found: ", path, "\n  Check the path and try again.", call. = FALSE)

  chunk_lines <- .check_count(chunk_lines, "chunk_lines", min = 1)

  ## The file is read `chunk_lines` records at a time. Each chunk's genotype
  ## text is turned into integer matrices (alleles, read depth, allele
  ## depths) and then discarded, so the whole text of a large VCF is never
  ## held in memory at once.
  vcf <- .open_vcf(path)
  on.exit(close(vcf$con), add = TRUE)
  header <- vcf$header
  if (length(header) < 10L)
    stop("The #CHROM header line has ", length(header), " columns, but a VCF with ",
         "genotypes needs at least 10 (CHROM ... FORMAT, then one per sample). ",
         "This file has no sample columns.", call. = FALSE)
  samples <- header[-(1:9)]

  parts <- list()
  n_read <- 0L
  lines <- vcf$pending
  repeat {
    lines <- lines[nzchar(lines)]
    if (length(lines)) {
      parts[[length(parts) + 1L]] <- .parse_vcf_chunk(lines, header, first_record = n_read + 1L)
      n_read <- n_read + length(lines)
    }
    lines <- readLines(vcf$con, n = chunk_lines)
    if (!length(lines)) break
  }
  if (!n_read) stop("No variant records found.", call. = FALSE)

  fields <- do.call(rbind, lapply(parts, `[[`, "fixed"))
  gt <- list(A1 = .bind_chunks(parts, "A1", length(samples)),
             A2 = .bind_chunks(parts, "A2", length(samples)),
             unrecognised = unique(unlist(lapply(parts, `[[`, "unrecognised"))))
  depth <- .bind_chunks(parts, "depth", length(samples))
  ad_ref <- .bind_chunks(parts, "ad_ref", length(samples))
  ad_alt <- .bind_chunks(parts, "ad_alt", length(samples))
  rm(parts)

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

  dimnames(gt$A1) <- dimnames(gt$A2) <- list(NULL, samples)
  if (!is.null(depth)) dimnames(depth) <- list(NULL, samples)
  if (!is.null(ad_ref)) dimnames(ad_ref) <- dimnames(ad_alt) <- list(NULL, samples)

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
    n_symbolic <- sum(vapply(alleles, function(a) any(.is_symbolic_allele(a)), logical(1)))
    if (n_symbolic)
      message("  ", .big(n_symbolic), " records have a symbolic ALT allele (*, <...>); ",
              "remove them first (see ?read_stacks_vcf, \"Assumptions\")")
  }

  ## `n_records_read` and `n_loci_read` are never changed by the filter_*()
  ## functions, so an analysis can tell that records or whole RAD loci were
  ## removed after reading (see the per-site checks in diversity_stats()).
  structure(list(A1 = gt$A1, A2 = gt$A2, locus = locus, locus_raw = locus_raw,
                 alleles = alleles, n_alleles = n_alleles, samples = samples,
                 fields = fields, depth = depth, ad_ref = ad_ref, ad_alt = ad_alt,
                 n_records_read = nrow(fields), n_loci_read = n_loci,
                 filter_log = .empty_filter_log()),
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
  steps <- .format_filter_log(x$filter_log)
  if (length(steps))
    cat("  filters applied since reading:\n", paste0("    ", steps, "\n"), sep = "")
  cat("  elements: $A1 $A2 (allele matrices), $locus, $locus_raw, $alleles, ",
      "$n_alleles, $samples, $fields",
      if (!is.null(x$depth)) ", $depth" else "",
      if (!is.null(x$ad_ref)) ", $ad_ref, $ad_alt" else "",
      if (!is.null(x$filter_log)) ", $filter_log" else "", "\n", sep = "")
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
  ## perl = TRUE: the same matches, several times faster on millions of cells.
  gt <- sub(":.*$", "", genotype_text, perl = TRUE)      # GT is the first subfield
  ok <- grepl("^[0-9]+[/|][0-9]+$", gt, perl = TRUE)
  A1 <- A2 <- matrix(NA_integer_, nrow(genotype_text), ncol(genotype_text))
  A1[ok] <- as.integer(sub("[/|].*$", "", gt[ok], perl = TRUE)) + 1L
  A2[ok] <- as.integer(sub("^.*[/|]", "", gt[ok], perl = TRUE)) + 1L
  unrecognised <- setdiff(unique(gt[!ok]), c("./.", ".|.", ".", "", "./", "/."))
  list(A1 = A1, A2 = A2, unrecognised = unrecognised)
}

## Not exported. Opens a VCF (plain or gzip-compressed; file() handles both)
## and reads up to its #CHROM line. Returns list(con = the open connection,
## header = the #CHROM line's column names, pending = any record lines read
## along with the header). The caller closes `con`.
.open_vcf <- function(path) {
  con <- file(path, "r")
  header <- NULL
  pending <- character(0)
  repeat {
    lines <- readLines(con, n = 1000L)
    if (!length(lines)) break
    header_line <- grep("^#CHROM", lines)
    if (length(header_line)) {
      header <- strsplit(sub("^#", "", lines[header_line[1]]), "\t", fixed = TRUE)[[1]]
      pending <- lines[-seq_len(header_line[1])]
      break
    }
  }
  if (is.null(header)) {
    close(con)
    stop("No '#CHROM' header line found. Is this a VCF?", call. = FALSE)
  }
  list(con = con, header = header, pending = pending)
}

## Not exported. The `position`-th colon-separated subfield of each cell of
## `cells` (a character vector). NA where the cell has fewer subfields (the VCF
## specification lets trailing subfields be dropped) or the subfield is "."
## or empty.
.subfield_at <- function(cells, position) {
  if (position < 2L) stop(".subfield_at() is for subfields after GT.", call. = FALSE)
  ## One pass: a cell with enough subfields loses at least "GT:" and so
  ## changes; a cell with too few does not match and comes back unchanged.
  out <- sub(paste0("^(?:[^:]*:){", position - 1L, "}([^:]*).*$"), "\\1", cells, perl = TRUE)
  out[out == cells | out %in% c(".", "")] <- NA_character_
  out
}

## Not exported. Parses one chunk of VCF record lines (`lines`, blank lines
## already removed). `first_record` is the number of the chunk's first record
## in the file, for error messages. Returns list(fixed = the 9 fixed columns,
## A1, A2, unrecognised (see .parse_gt()), depth, ad_ref, ad_alt), where:
##   depth   read depth of each call: its DP, or the sum of its AD where DP is
##           absent; NA when neither can be read
##   ad_ref  reads supporting REF (first AD value; "." inside AD counts as 0);
##           NA when the call has no AD
##   ad_alt  reads supporting the call's own ALT allele(s), each distinct ALT
##           allele counted once (a 1/1 call counts allele 1's reads once);
##           NA when the call has no AD
## depth, ad_ref and ad_alt are NULL when no FORMAT in the chunk has DP or AD.
## These are what filter_genotype_depth() and filter_low_conf_alt() use, so
## the genotype text itself need not be kept.
.parse_vcf_chunk <- function(lines, header, first_record) {
  second_header <- which(startsWith(lines, "#CHROM"))
  if (length(second_header))
    stop("Expected exactly one '#CHROM' header line; found a second one after record ",
         first_record + second_header[1] - 2L, ".", call. = FALSE)
  split_lines <- strsplit(lines, "\t", fixed = TRUE)
  n_fields <- lengths(split_lines)
  bad_len <- which(n_fields != length(header))
  if (length(bad_len)) {
    ## Filling a matrix from unequal-length records would shift every later
    ## column without any warning, so this must stop.
    stop(sprintf("Malformed VCF: record %d has %d tab-separated field(s), expected %d ",
                 first_record + bad_len[1] - 1L, n_fields[bad_len[1]], length(header)),
         "to match the #CHROM header line. The file is truncated, corrupted, or ",
         "not a VCF. Fix or regenerate the file and re-run.", call. = FALSE)
  }
  all_fields <- matrix(unlist(split_lines, use.names = FALSE), nrow = length(lines),
                       byrow = TRUE)
  rm(split_lines)
  fixed <- all_fields[, 1:9, drop = FALSE]
  colnames(fixed) <- header[1:9]
  .check_gt_first(fixed[, "FORMAT"], first_record)
  cells <- all_fields[, -(1:9), drop = FALSE]
  rm(all_fields)
  gt <- .parse_gt(cells)

  n_rec <- nrow(cells)
  n_samp <- ncol(cells)
  formats <- fixed[, "FORMAT"]
  unique_formats <- unique(formats)
  format_parts <- strsplit(unique_formats, ":", fixed = TRUE)
  has_tag <- function(tag) any(vapply(format_parts, function(f) tag %in% f, logical(1)))
  if (!has_tag("DP") && !has_tag("AD"))
    return(list(fixed = fixed, A1 = gt$A1, A2 = gt$A2, unrecognised = gt$unrecognised,
                depth = NULL, ad_ref = NULL, ad_alt = NULL))

  depth <- ad_ref <- ad_alt <- matrix(NA_integer_, n_rec, n_samp)
  for (f in seq_along(format_parts)) {
    rows <- which(formats == unique_formats[f])
    block <- cells[rows, , drop = FALSE]
    dp_pos <- match("DP", format_parts[[f]])
    ad_pos <- match("AD", format_parts[[f]])
    dp <- if (is.na(dp_pos)) rep(NA_integer_, length(block))
          else suppressWarnings(as.integer(.subfield_at(block, dp_pos)))
    if (!is.na(ad_pos)) {
      ad <- .subfield_at(block, ad_pos)
      usable <- !is.na(ad)
      ## One column per allele: the reads for allele a at each call; "." or
      ## a missing value inside AD counts as 0.
      n_values <- ifelse(usable, nchar(gsub("[^,]", "", ad)) + 1L, 0L)
      a1 <- gt$A1[rows, , drop = FALSE]
      a2 <- gt$A2[rows, , drop = FALSE]
      k <- max(1L, n_values)
      reads <- matrix(0L, length(ad), k)
      for (a in seq_len(k)) {
        at <- usable & n_values >= a
        if (!any(at)) next
        value <- sub(paste0("^(?:[^,]*,){", a - 1L, "}([^,]*).*$"), "\\1", ad[at], perl = TRUE)
        value <- suppressWarnings(as.integer(value))
        value[is.na(value)] <- 0L
        reads[at, a] <- value
      }
      reads_of <- function(allele) {
        out <- integer(length(allele))
        ok <- !is.na(allele) & allele >= 1L & allele <= k
        out[ok] <- reads[cbind(which(ok), allele[ok])]
        out
      }
      first_alt <- ifelse(!is.na(a1) & a1 > 1L, reads_of(a1), 0L)
      second_alt <- ifelse(!is.na(a2) & a2 > 1L & a2 != a1, reads_of(a2), 0L)
      ref_reads <- reads[, 1L]
      ref_reads[!usable] <- NA_integer_
      alt_reads <- first_alt + second_alt
      alt_reads[!usable] <- NA_integer_
      ad_ref[rows, ] <- ref_reads
      ad_alt[rows, ] <- alt_reads
      ad_total <- as.integer(rowSums(reads))
      dp <- ifelse(is.na(dp) & usable, ad_total, dp)
    }
    depth[rows, ] <- dp
  }
  if (!has_tag("AD")) ad_ref <- ad_alt <- NULL
  list(fixed = fixed, A1 = gt$A1, A2 = gt$A2, unrecognised = gt$unrecognised,
       depth = depth, ad_ref = ad_ref, ad_alt = ad_alt)
}

## Not exported. Stacks the records x samples matrix `name` of every chunk.
## A chunk without it (e.g. no DP or AD in its FORMAT) contributes NA rows;
## NULL when no chunk has it.
.bind_chunks <- function(parts, name, n_samp) {
  present <- !vapply(parts, function(p) is.null(p[[name]]), logical(1))
  if (!any(present)) return(NULL)
  do.call(rbind, lapply(parts, function(p)
    if (is.null(p[[name]])) matrix(NA_integer_, nrow(p$fixed), n_samp) else p[[name]]))
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
## previous record or more than `window_bp` beyond it. Each locus is labeled
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
#' The file is read with the same rules as Stacks: blank lines and lines
#' starting with `#` are skipped, and a third column (Stacks' optional
#' population group) is allowed and ignored. A line with only one column,
#' more than three, or columns separated by spaces instead of a TAB stops with
#' the line number and its text.
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
  pm <- .parse_popmap_lines(path, verbose)
  if (!nrow(pm)) stop("The popmap file has no samples.", call. = FALSE)
  ## Keep only samples in the VCF BEFORE grouping, so populations come out in
  ## the order they first appear among the retained samples.
  if (!is.null(samples)) {
    .report_popmap_not_in_vcf(pm$sample, samples, verbose)
    pm <- pm[pm$sample %in% samples, , drop = FALSE]
  }
  pops <- split(pm$sample, factor(pm$pop, levels = unique(pm$pop)))
  .clean_pops(pops, samples, verbose)
}

## Not exported. Reads a popmap file with the same rules as Stacks
## (MetaPopInfo::init_popmap(), Stacks 2.68): blank lines and lines starting
## with "#" are skipped, a Windows line ending is removed, fields are
## TAB-separated and trimmed, and a line has SAMPLE, POPULATION and an
## optional third GROUP column (which this package does not use). Anything
## else stops with the line number and its text, so a popmap separated by
## spaces, or with a missing population, is caught here with a message that
## says what is wrong. Every value stays text, so a sample ID such as "001"
## still matches the VCF. Returns data.frame(sample, pop).
.parse_popmap_lines <- function(path, verbose = TRUE) {
  lines <- sub("\r$", "", readLines(path, warn = FALSE))
  line_number <- seq_along(lines)
  content <- nzchar(trimws(lines)) & !startsWith(lines, "#")
  fields <- lapply(strsplit(lines[content], "\t", fixed = TRUE), trimws)
  fields <- lapply(fields, function(f) {
    while (length(f) && !nzchar(f[length(f)])) f <- f[-length(f)]   # trailing tabs
    f
  })
  n_fields <- lengths(fields)
  empty_part <- vapply(fields, function(f) length(f) >= 2 && !all(nzchar(f[1:2])), logical(1))
  bad <- which(n_fields < 2 | n_fields > 3 | empty_part)
  if (length(bad)) {
    shown <- utils::head(bad, 5)
    spaces <- grepl(" ", lines[content][shown]) & !grepl("\t", lines[content][shown])
    stop("Malformed popmap ", path, ": each line must be SAMPLE<TAB>POPULATION ",
         "(optionally <TAB>GROUP, as Stacks allows). Problem line(s):\n",
         paste0("  line ", line_number[content][shown], ": \"", lines[content][shown], "\"",
                ifelse(n_fields[shown] < 2, " (no population column)",
                       ifelse(n_fields[shown] > 3, sprintf(" (%d columns)", n_fields[shown]),
                              " (an empty sample or population)")),
                collapse = "\n"),
         if (any(spaces)) "\n  The columns look separated by spaces; they must be separated by a TAB."
         else "",
         if (length(bad) > 5) sprintf("\n  ... and %d more.", length(bad) - 5) else "",
         call. = FALSE)
  }
  if (any(n_fields == 3L))
    .inform(verbose, "  popmap has a third (group) column; it is not used.")
  data.frame(sample = vapply(fields, `[`, character(1), 1L),
             pop = vapply(fields, `[`, character(1), 2L), stringsAsFactors = FALSE)
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
  if (is.data.frame(pops))
    stop("`popmap` is a data frame; its columns would be read as populations. Pass the ",
         "popmap file's path, or a list of sample IDs per population, e.g. ",
         "split(as.character(df[[1]]), df[[2]]) for sample IDs in column 1 and ",
         "populations in column 2.", call. = FALSE)
  if (!is.list(pops) || is.null(names(pops)) || any(!nzchar(names(pops))))
    stop("`popmap` must be the path to a popmap file, or a named list of sample ",
         "IDs per population as returned by read_popmap() (got an object of class ",
         paste(class(pops), collapse = "/"), ").", call. = FALSE)
  not_text <- names(pops)[!vapply(pops, is.character, logical(1))]
  if (length(not_text))
    stop("`popmap`: each population must be a character vector of sample IDs, but ",
         paste(not_text, collapse = ", "), " is not (class ",
         paste(unique(vapply(pops[not_text], function(x) class(x)[1], character(1))),
               collapse = "/"),
         "). Convert with lapply(popmap, as.character).", call. = FALSE)
  if (anyDuplicated(names(pops)))
    stop("Population names must be unique; repeated: ",
         paste(unique(names(pops)[duplicated(names(pops))]), collapse = ", "), call. = FALSE)
  if (!is.null(samples)) {
    .report_popmap_not_in_vcf(unlist(pops, use.names = FALSE), samples, verbose)
    pops <- lapply(pops, function(ids) ids[ids %in% samples])
    pops <- pops[lengths(pops) > 0]
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

## Not exported. Compares the popmap's sample IDs (`ids`) with the VCF's
## (`samples`). A popmap sample missing from the VCF is left out of every
## analysis, which quietly makes its population smaller, and is most often a
## typo or a different popmap, so a message names up to 10 of them. If no
## popmap sample is in the VCF, stops and shows a few names from each, so a
## systematic difference (a suffix such as ".sorted", upper vs lower case)
## can be seen.
.report_popmap_not_in_vcf <- function(ids, samples, verbose = TRUE) {
  missing_ids <- unique(ids[!ids %in% samples])
  if (length(missing_ids) == length(unique(ids)))
    stop("No popmap sample names match the VCF.\n",
         "  VCF samples:    ", paste(utils::head(samples, 3), collapse = ", "),
         if (length(samples) > 3) ", ..." else "", "\n",
         "  popmap samples: ", paste(utils::head(unique(ids), 3), collapse = ", "),
         if (length(unique(ids)) > 3) ", ..." else "", "\n",
         "  The names must match exactly (spelling, upper/lower case, suffixes).", call. = FALSE)
  if (length(missing_ids))
    .inform(verbose, "  ", length(missing_ids), " popmap sample(s) not in the VCF, so left out: ",
            paste(utils::head(missing_ids, 10), collapse = ", "),
            if (length(missing_ids) > 10) ", ..." else "")
  invisible(missing_ids)
}

## Not exported. Stops, naming them, if any population has fewer than 2
## individuals. `why` says what the calling function needs 2 individuals for.
.stop_tiny_pops <- function(pops, why) {
  tiny <- names(pops)[lengths(pops) < 2]
  if (length(tiny))
    stop("Population(s) with fewer than 2 individuals: ", paste(tiny, collapse = ", "),
         "\n  ", why, " Drop these populations from the popmap or merge them.", call. = FALSE)
  invisible(NULL)
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

## Not exported. Is H a haplotype VCF (records that span several sites of a
## RAD tag) rather than a SNP VCF (one site per record)? Decided by allele
## LENGTH only: a record whose alleles are all single bases is one site,
## whatever the file is called and however many alleles it has, while Stacks
## writes a haplotype allele as the bases at every SNP of the tag ("AC",
## "CA"), so any tag with two or more SNPs gives multi-base alleles. The
## number of alleles is deliberately NOT used: a SNP VCF from GATK, freebayes
## or ipyrad can contain a few multi-allelic SNPs (e.g. A -> C,T), and
## treating those as haplotypes would switch off the per-site values. Indels
## and MNPs also have multi-base alleles, which is why ?read_stacks_vcf asks
## for them to be removed from a SNP VCF. Symbolic alleles ("*", "<NON_REF>")
## are codes, not sequences, so their length is ignored.
.is_haplotype_H <- function(H) {
  alleles <- unlist(H$alleles, use.names = FALSE)
  any(nchar(alleles[!.is_symbolic_allele(alleles)]) > 1L)
}

## Not exported. TRUE for a symbolic VCF allele, which is a code rather than
## a DNA sequence: "*" (an allele missing because of an overlapping deletion)
## or anything in angle brackets ("<*>", "<NON_REF>", "<DEL>").
.is_symbolic_allele <- function(alleles) alleles == "*" | grepl("^<.*>$", alleles)

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
## population number, then one column per locus with each genotype as one
## integer, 2 digits per allele (alleles 1 and 12 -> 112). Used only for
## Weir & Goudet's beta and the optional hierfstat cross-checks.
##
## WHY 2 DIGITS, AND WHY ALLELES ARE RENUMBERED. hierfstat does not know how
## many digits each allele has; its getal() (hierfstat 0.5.11) guesses from
## the numbers. With 3 digits per allele (1 and 12 -> 1012) and every number
## below 10000, it guesses 2 digits and reads 1012 as alleles 10 and 12. On
## haplotype data with 10 or more alleles that silently scrambled every
## genotype, and so beta. With 2 digits per allele its guess is always right:
## every number is 101-9999, and it picks 3 digits only when every allele is
## both >= 10 and < 10. So each record's alleles are first renumbered 1, 2,
## ... over the alleles these individuals actually carry (allele labels do not
## change any hierfstat statistic), which keeps them within 99. A record with
## more than 99 alleles needs 3 digits; if hierfstat would then guess wrong,
## this returns NULL and the caller skips hierfstat.
.to_hierfstat_df <- function(H, pops, rows = seq_len(nrow(H$A1))) {
  ids <- unlist(pops, use.names = FALSE)
  a <- H$A1[rows, ids, drop = FALSE]
  b <- H$A2[rows, ids, drop = FALSE]
  ## new_number[r, j]: allele j's new number at record r = how many alleles
  ## up to and including j are carried there.
  carried <- .allele_counts(a, b) > 0L                     # records x alleles
  new_number <- carried
  storage.mode(new_number) <- "integer"
  running <- integer(nrow(carried))
  for (j in seq_len(ncol(carried))) {
    running <- running + carried[, j]
    new_number[, j] <- running
  }
  renumber <- function(A) {
    typed <- !is.na(A)
    A[typed] <- new_number[cbind(row(A)[typed], A[typed])]
    A
  }
  a <- renumber(a)
  b <- renumber(b)
  digits <- if (max(c(0L, running)) <= 99L) 100L else 1000L
  genotype_codes <- t(pmin(a, b) * digits + pmax(a, b))    # individuals x loci
  if (digits == 1000L && !.hierfstat_reads_3_digits(genotype_codes)) return(NULL)
  dat <- data.frame(pop = rep(seq_along(pops), lengths(pops)), genotype_codes,
                    row.names = NULL)
  names(dat)[-1] <- paste0("L", seq_along(rows))
  dat
}

## Not exported. TRUE when hierfstat's getal() would read these
## 3-digits-per-allele genotype numbers with 3 digits (its rule, restated).
.hierfstat_reads_3_digits <- function(codes) {
  top <- max(codes, na.rm = TRUE)
  top >= 10000 ||
    (min(codes %/% 100, na.rm = TRUE) >= 10 && max(codes %% 100, na.rm = TRUE) < 10)
}
