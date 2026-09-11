###############################################################################
#
#  R/vcf_io.R  --  shared reader for Stacks `populations.haps.vcf`
#
#  Used by diversity_stats() and het_between_pops() so all callers parse the
#  file identically.
#
#  WHAT A HAPS VCF IS
#  ------------------
#  `populations.haps.vcf` holds one record per RAD locus. REF is one haplotype
#  (the concatenated variable sites of the tag) and ALT is a comma-separated
#  list of the others, so records are MULTI-ALLELIC. Genotypes are indices into
#  c(REF, ALT): 0/0, 0/1, 1/2, ... and ./. for missing.
#
#  This is why haplotype loci suit some analyses better than SNPs:
#    - one locus per RAD tag, so loci are independent by construction and no
#      block bootstrap is needed to handle within-tag linkage
#    - multi-allelic, so the biallelic pathology where a single minor allele
#      copy forces FIS to exactly zero does not arise
#    - allelic richness becomes a real quantity rather than 1-or-2
#  And why they suit others worse:
#    - one sequencing error anywhere in the tag creates a spurious haplotype,
#      inflating richness and gene diversity, worse at low depth
#    - a restriction-site null allele removes the whole haplotype, so allele
#      dropout is NOT mitigated
#    - paralogs generate extra haplotypes rather than just excess heterozygosity
#
#  Requires only base R.
#
###############################################################################

#' Read a Stacks haplotype or SNP VCF
#'
#' Parses `populations.haps.vcf` or `populations.snps.vcf` (gzip-compressed or
#' not) into genotype matrices plus locus/allele metadata. Used by
#' [diversity_stats()] and [het_between_pops()].
#'
#' @param path Path to the VCF file (`.vcf` or `.vcf.gz`).
#' @param verbose Print progress/summary messages. Default `TRUE`.
#' @return A list with elements `A1`, `A2` (allele-index matrices, one row per
#'   record, one column per sample), `locus` (unique per-record locus names),
#'   `locus_raw` (locus grouping, shared across SNPs on one RAD tag),
#'   `alleles` (list of allele strings per record), `n_alleles`, `samples`,
#'   and `fields` (the raw parsed VCF field matrix).
#' @export
read_haps_vcf <- function(path, verbose = TRUE) {

  if (grepl("\\.gz$", path)) {
    lines <- readLines(gzfile(path))
  } else {
    lines <- readLines(path)
  }
  hdr_i <- grep("^#CHROM", lines)
  if (length(hdr_i) != 1L)
    stop("Expected exactly one '#CHROM' header line; found ", length(hdr_i), ".")

  hdr  <- strsplit(sub("^#", "", lines[hdr_i]), "\t", fixed = TRUE)[[1]]
  body <- lines[(hdr_i + 1L):length(lines)]
  body <- body[nzchar(body)]
  if (!length(body)) stop("No variant records found.")

  ## Every record must have exactly as many tab-separated fields as the
  ## header. rbind() on an uneven list of vectors RECYCLES silently -- no
  ## warning at all when the short length divides evenly into the target --
  ## which would misalign every column after the truncated record with no
  ## indication anything went wrong. Check explicitly and fail loudly instead.
  split_body <- strsplit(body, "\t", fixed = TRUE)
  nf <- lengths(split_body)
  bad_len <- which(nf != length(hdr))
  if (length(bad_len))
    stop(sprintf(
      "Malformed VCF: record %d of %d has %d tab-separated field(s), expected %d ",
      bad_len[1], length(body), nf[bad_len[1]], length(hdr)),
      "to match the #CHROM header line. The file is truncated, corrupted, or ",
      "not a Stacks populations VCF. Fix or regenerate the file and re-run.")

  f <- do.call(rbind, split_body)
  colnames(f) <- hdr
  samples <- hdr[10:length(hdr)]

  ## The VCF spec requires GT, when present, to be the FIRST FORMAT subfield.
  ## Every genotype below is decoded by assuming that. Check it explicitly:
  ## otherwise a FORMAT field like "DP:GT" would be misread as valid-looking
  ## but wrong allele indices, most of which would slip past the overflow
  ## check further down.
  fmt <- f[, "FORMAT"]
  bad_fmt <- which(!grepl("^GT(:|$)", fmt))
  if (length(bad_fmt))
    stop(sprintf(
      "Malformed VCF: FORMAT field at record %d is \"%s\", not GT-first. ",
      bad_fmt[1], fmt[bad_fmt[1]]),
      "This parser requires GT to be the first FORMAT subfield, per the VCF ",
      "spec and standard Stacks output. Check how the file was produced.")

  ## Locus identity. De novo Stacks writes CHROM = "un" and puts the locus in
  ## ID; reference-aligned runs put a scaffold in CHROM. Either way one record
  ## should be one locus, but this is asserted rather than assumed.
  locus <- f[, "ID"]
  if (all(locus == ".") || any(!nzchar(locus))) locus <- paste0(f[, "CHROM"], "_", f[, "POS"])
  locus <- sub(":.*$", "", locus)
  ## locus_raw preserves the grouping: several SNPs on one RAD tag share an ID.
  ## `locus` is uniquified for use as marker names. ALWAYS group by locus_raw --
  ## grouping by `locus` makes one-SNP-per-locus thinning a silent no-op.
  locus_raw <- locus
  n_dup <- sum(duplicated(locus_raw))
  if (n_dup) {
    if (verbose) message(sprintf("  %s records share a locus ID with another (%s distinct loci) -- expected for a SNP VCF, NOT for a haps VCF",
                                 format(n_dup, big.mark = ","),
                                 format(length(unique(locus_raw)), big.mark = ",")))
    locus <- make.unique(locus)
  }

  ## Alleles per record: REF plus the comma-separated ALT list.
  alt      <- f[, "ALT"]
  alt_list <- ifelse(alt == "." | !nzchar(alt), "", alt)
  alleles  <- Map(function(r, a) if (nzchar(a)) c(r, strsplit(a, ",", fixed = TRUE)[[1]]) else r,
                  f[, "REF"], alt_list)
  n_alleles <- lengths(alleles)

  ## Genotypes -> pair of 1-based allele indices per individual.
  gt <- sub(":.*$", "", f[, 10:ncol(f), drop = FALSE])
  gt <- gsub("|", "/", gt, fixed = TRUE)

  parts <- strsplit(as.vector(gt), "/", fixed = TRUE)
  a1 <- suppressWarnings(as.integer(vapply(parts, function(x) x[1], character(1)))) + 1L
  a2 <- suppressWarnings(as.integer(vapply(parts, function(x) x[2], character(1)))) + 1L
  ## anything not a clean diploid call (./., ., haploid, phased-missing) -> NA
  bad <- is.na(a1) | is.na(a2) | lengths(parts) != 2L
  a1[bad] <- NA_integer_; a2[bad] <- NA_integer_
  A1 <- matrix(a1, nrow = nrow(gt), dimnames = list(NULL, samples))
  A2 <- matrix(a2, nrow = nrow(gt), dimnames = list(NULL, samples))

  ## An allele index above the number of declared alleles means the record and
  ## the genotypes disagree -- fail loudly rather than compute nonsense.
  over <- which(!is.na(A1) & (A1 > n_alleles) | !is.na(A2) & (A2 > n_alleles))
  if (length(over))
    stop(length(over), " genotype calls reference an allele index beyond REF+ALT. ",
         "The VCF is malformed or the FORMAT field is not GT-first.")

  if (verbose) {
    message(sprintf("  %s loci, %d samples", format(nrow(f), big.mark = ","), length(samples)))
    message(sprintf("  haplotypes per locus: median %d, max %d; %.1f%% of loci are multi-allelic (>2)",
                    as.integer(stats::median(n_alleles)), max(n_alleles),
                    100 * mean(n_alleles > 2)))
    message(sprintf("  missing genotype rate: %.2f%%", 100 * mean(is.na(A1))))
    unrec <- setdiff(unique(as.vector(gt)[bad]), c("./.", ".", "", "./", "/."))
    if (length(unrec))
      message("  treated as missing: ", paste(utils::head(unrec, 8), collapse = ", "))
  }

  list(A1 = A1, A2 = A2, locus = locus, locus_raw = locus_raw, alleles = alleles,
       n_alleles = n_alleles, samples = samples, fields = f)
}

## Not exported. Shared by diversity_stats() and het_between_pops() so both
## accept EITHER a path to a VCF file (the original behavior -- this reads
## it with read_haps_vcf()) OR an H list someone already built themselves,
## e.g. by running read_haps_vcf() and then one or more filter_*() functions
## from R/filter_loci.R. This is what lets a filtered/cleaned dataset go
## straight into either pipeline function without writing an intermediate
## VCF file to disk first.
##
## `vcf_file` is a "character" (plain text) value when it's a file path, and
## a "list" (the shape read_haps_vcf() returns) when it's already-parsed
## data -- is.character()/is.list() below just tell those two cases apart.
.resolve_H <- function(vcf_file, verbose = TRUE) {
  if (is.character(vcf_file)) return(read_haps_vcf(vcf_file, verbose = verbose))
  if (!is.list(vcf_file))
    stop("vcf_file must be either a path to a VCF file, or the list returned ",
         "by read_haps_vcf() (optionally passed through one or more filter_*() ",
         "functions first). Got an object of class: ", paste(class(vcf_file), collapse = "/"))
  ## A hand-built or corrupted list could be missing pieces read_haps_vcf()
  ## always includes -- check for those up front so a confusing error deep
  ## inside diversity_stats()/het_between_pops() doesn't happen instead.
  required <- c("A1", "A2", "locus", "locus_raw", "alleles", "n_alleles", "samples")
  missing_el <- setdiff(required, names(vcf_file))
  if (length(missing_el))
    stop("vcf_file looks like a list, but is missing element(s) that read_haps_vcf() ",
         "always includes: ", paste(missing_el, collapse = ", "),
         ". Pass the object returned by read_haps_vcf() (optionally filtered), not ",
         "something else.")
  if (!is.matrix(vcf_file$A1) || !is.matrix(vcf_file$A2) ||
      !identical(dim(vcf_file$A1), dim(vcf_file$A2)))
    stop("vcf_file$A1 and vcf_file$A2 must both be matrices of the same size ",
         "(as read_haps_vcf() always produces).")
  vcf_file
}

#' Read a Stacks-style popmap
#'
#' Parses a two-column, no-header popmap TSV (`sample_id <TAB> population`)
#' and splits sample IDs by population, restricted to samples present in
#' `samples`.
#'
#' @param path Path to the popmap TSV.
#' @param samples Character vector of sample IDs to keep (typically the VCF's
#'   sample columns).
#' @param verbose Print progress/summary messages. Default `TRUE`.
#' @return A named list of character vectors, one per population, each
#'   holding that population's sample IDs.
#' @export
read_popmap <- function(path, samples, verbose = TRUE) {
  pm <- utils::read.delim(path, header = FALSE, stringsAsFactors = FALSE,
                   col.names = c("sample", "pop"))
  pm <- pm[pm$sample %in% samples, , drop = FALSE]
  if (!nrow(pm)) stop("No popmap sample names match the VCF.")
  ## A sample listed twice ends up in BOTH populations and is counted twice,
  ## silently corrupting every between-population statistic. Fail loudly.
  dupd <- unique(pm$sample[duplicated(pm$sample)])
  if (length(dupd))
    stop(length(dupd), " sample(s) appear more than once in the popmap: ",
         paste(utils::head(dupd, 10), collapse = ", "),
         "\n  Each individual must be assigned to exactly one population.")
  drop <- setdiff(samples, pm$sample)
  if (length(drop) && verbose)
    message("  ", length(drop), " VCF samples absent from the popmap; excluded.")
  pops <- split(pm$sample, factor(pm$pop, levels = unique(pm$pop)))
  if (verbose)
    message("  populations: ",
            paste(sprintf("%s (n=%d)", names(pops), lengths(pops)), collapse = ", "))
  pops
}
