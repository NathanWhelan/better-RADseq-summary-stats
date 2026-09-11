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
#' `read_haps_vcf()` is this function's former name (it reads SNP VCFs too),
#' kept as an alias for one release.
#'
#' Other VCFs (ipyrad, dDocent/freebayes, GATK, ...) are read the same way,
#' provided GT is the first FORMAT field; `locus_from` below says how to
#' group their records into RAD loci.
#'
#' @param path Path to the VCF file (`.vcf` or `.vcf.gz`).
#' @param verbose Print progress/summary messages. Default `TRUE`.
#' @param locus_from How to tell which RAD locus each record belongs to --
#'   the unit the block jackknife and bootstrap resample, so that SNPs on one
#'   RAD tag move together. `"auto"` (the default) uses the ID column when
#'   every record has one (Stacks writes the locus number, or
#'   `locus:column:strand` in a SNP VCF), and `"window"` otherwise.
#'   `"ID"`: the ID column, up to its first `:`. `"CHROM"`: one locus per
#'   CHROM value, for de novo formats that put the locus there (some Stacks 2
#'   versions, ipyrad, dDocent). `"window"`: records on the same CHROM within
#'   `window_bp` of the previous record form one locus -- which gives one
#'   locus per CHROM for the de novo formats (a RAD locus is far shorter than
#'   1 kb) and one per RAD locus, or run of adjacent loci, for
#'   reference-aligned data. For a reference-aligned VCF whose IDs are
#'   unique per SNP (e.g. dbSNP rs numbers), choose `"window"` explicitly.
#' @param window_bp Largest gap, in base pairs, between consecutive records of
#'   one locus under `locus_from = "window"`. Default `1000`.
#' @return A list with elements `A1`, `A2` (allele-index matrices, one row per
#'   record, one column per sample), `locus` (unique per-record locus names),
#'   `locus_raw` (locus grouping, shared across SNPs on one RAD tag),
#'   `alleles` (list of allele strings per record), `n_alleles`, `samples`,
#'   and `fields` (the raw parsed VCF field matrix).
#' @examples
#' H <- read_stacks_vcf(system.file("extdata", "small.haps.vcf",
#'                                  package = "RADdiversity"))
#' dim(H$A1)          # records x samples
#' H$alleles[[1]]     # the haplotype alleles of the first RAD locus
#' H$A1[1:3, 1:4]     # first allele of each genotype (1 = REF, 2 = first ALT, ...)
#' @export
read_stacks_vcf <- function(path, verbose = TRUE, locus_from = "auto", window_bp = 1000) {
  if (!(length(locus_from) == 1L && locus_from %in% c("auto", "ID", "CHROM", "window")))
    stop("locus_from must be one of \"auto\", \"ID\", \"CHROM\", \"window\" (got: ",
         paste(locus_from, collapse = ", "), ").")
  if (!(is.numeric(window_bp) && length(window_bp) == 1L && window_bp >= 0))
    stop("window_bp must be a single non-negative number of base pairs.")

  if (grepl("\\.gz$", path)) {
    lines <- readLines(gzfile(path))
  } else {
    lines <- readLines(path)
  }
  hdr_i <- grep("^#CHROM", lines)
  if (length(hdr_i) != 1L)
    stop("Expected exactly one '#CHROM' header line; found ", length(hdr_i), ".")

  hdr  <- strsplit(sub("^#", "", lines[hdr_i]), "\t", fixed = TRUE)[[1]]
  ## Everything after the header. Not (hdr_i + 1):length(lines): when the
  ## header is the last line that sequence counts DOWN, picking up an NA and
  ## the header itself, and a header-only file then failed as "malformed"
  ## instead of "no variant records".
  body <- lines[-seq_len(hdr_i)]
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

  ## Which RAD locus each record belongs to (`locus_raw`) -- the unit the
  ## block jackknife and bootstrap resample, so SNPs on one RAD tag move
  ## together. Formats put it in different places:
  ##   Stacks (most versions)  ID = "locus:column:strand" (SNP VCF) or "locus"
  ##   Stacks 2, some builds   ID = "."; CHROM = catalog locus, POS = column
  ##   ipyrad                  ID = "."; CHROM = "RAD_<n>", one per locus
  ##   reference-aligned       CHROM = chromosome/scaffold, POS = genomic
  ## See `locus_from` in the documentation. Getting this wrong is silent: with
  ## every SNP its own "locus", the resampling treats linked SNPs as
  ## independent. `locus` (unique per record) is used as the marker name;
  ## ALWAYS group by locus_raw -- grouping by `locus` would make
  ## one-SNP-per-locus thinning a silent no-op.
  id <- f[, "ID"]
  id_ok <- all(nzchar(id) & id != ".")
  rule <- if (locus_from == "auto") (if (id_ok) "ID" else "window") else locus_from
  if (rule == "ID" && !id_ok)
    stop("locus_from = \"ID\", but some records have no ID (\".\"). Use ",
         "locus_from = \"CHROM\" or \"window\" (see ?read_stacks_vcf).")
  locus_raw <- switch(rule,
    ID     = sub(":.*$", "", id),
    CHROM  = f[, "CHROM"],
    window = .window_blocks(f[, "CHROM"], f[, "POS"], window_bp))
  locus <- if (id_ok) sub(":.*$", "", id) else paste0(f[, "CHROM"], "_", f[, "POS"])
  locus <- make.unique(locus)
  if (verbose) {
    n_loc <- length(unique(locus_raw))
    message(sprintf("  RAD loci from %s: %s records on %s loci (%.2f per locus)",
                    switch(rule, ID = "the ID column", CHROM = "CHROM",
                           window = sprintf("CHROM + POS (records within %s bp)",
                                            format(window_bp, big.mark = ","))),
                    format(nrow(f), big.mark = ","), format(n_loc, big.mark = ","),
                    nrow(f) / n_loc))
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

#' @rdname read_stacks_vcf
#' @export
read_haps_vcf <- function(path, verbose = TRUE, locus_from = "auto", window_bp = 1000)
  read_stacks_vcf(path, verbose = verbose, locus_from = locus_from, window_bp = window_bp)

## Not exported. Groups records into loci by position: sorted within each
## CHROM, a record starts a new locus when it is on a different CHROM from the
## previous record or more than `window_bp` beyond it. Each locus is labelled
## "CHROM:POS" of its first record; labels come back in the input order.
.window_blocks <- function(chrom, pos, window_bp) {
  pos <- suppressWarnings(as.numeric(pos))
  if (anyNA(pos))
    stop("Some POS values are not numbers, so records cannot be grouped by ",
         "position. Use locus_from = \"ID\" or \"CHROM\".")
  o <- order(chrom, pos)
  cs <- chrom[o]; ps <- pos[o]
  new_locus <- c(TRUE, cs[-1] != cs[-length(cs)] | diff(ps) > window_bp)
  first <- paste0(cs, ":", format(ps, scientific = FALSE, trim = TRUE))[new_locus]
  out <- character(length(chrom))
  out[o] <- first[cumsum(new_locus)]
  out
}

## Not exported. Shared by diversity_stats() and het_between_pops() so both
## accept EITHER a path to a VCF file (the original behavior -- this reads
## it with read_stacks_vcf()) OR an H list someone already built themselves,
## e.g. by running read_stacks_vcf() and then one or more filter_*() functions
## from R/filter_loci.R. This is what lets a filtered/cleaned dataset go
## straight into either pipeline function without writing an intermediate
## VCF file to disk first.
##
## `vcf_file` is a "character" (plain text) value when it's a file path, and
## a "list" (the shape read_stacks_vcf() returns) when it's already-parsed
## data -- is.character()/is.list() below just tell those two cases apart.
.resolve_H <- function(vcf_file, verbose = TRUE) {
  if (is.character(vcf_file)) return(read_stacks_vcf(vcf_file, verbose = verbose))
  if (!is.list(vcf_file))
    stop("vcf_file must be either a path to a VCF file, or the list returned ",
         "by read_stacks_vcf() (optionally passed through one or more filter_*() ",
         "functions first). Got an object of class: ", paste(class(vcf_file), collapse = "/"))
  ## A hand-built or corrupted list could be missing pieces read_stacks_vcf()
  ## always includes -- check for those up front so a confusing error deep
  ## inside diversity_stats()/het_between_pops() doesn't happen instead.
  required <- c("A1", "A2", "locus", "locus_raw", "alleles", "n_alleles", "samples")
  missing_el <- setdiff(required, names(vcf_file))
  if (length(missing_el))
    stop("vcf_file looks like a list, but is missing element(s) that read_stacks_vcf() ",
         "always includes: ", paste(missing_el, collapse = ", "),
         ". Pass the object returned by read_stacks_vcf() (optionally filtered), not ",
         "something else.")
  if (!is.matrix(vcf_file$A1) || !is.matrix(vcf_file$A2) ||
      !identical(dim(vcf_file$A1), dim(vcf_file$A2)))
    stop("vcf_file$A1 and vcf_file$A2 must both be matrices of the same size ",
         "(as read_stacks_vcf() always produces).")
  vcf_file
}

## Not exported. Input checks shared by diversity_stats(), het_between_pops()
## and differentiation_stats(), run before any (potentially slow) reading.
.check_run_inputs <- function(vcf_file, popmap_f, stem) {
  if (is.character(vcf_file) && !file.exists(vcf_file))
    stop("VCF file not found: ", vcf_file, "\n  Check the path and try again.")
  if (!is.character(vcf_file) && is.null(stem))
    stop("vcf_file is an already-parsed list rather than a file path, so its ",
         "filename can't be used to name the output files. Pass stem ",
         "explicitly, e.g. stem = \"haps\" or stem = \"snps\", matching which ",
         "VCF this data came from.")
  if (!file.exists(popmap_f))
    stop("Popmap file not found: ", popmap_f, "\n  Check the path and try again.")
  invisible(NULL)
}

## Not exported. Is H a haplotype VCF (one multi-allelic record per RAD tag)
## rather than a SNP VCF (one site per record)? Either signal is decisive:
## Stacks writes multi-nucleotide alleles ("AC", "CA") for haplotype records
## and single bases for SNP records, and haplotype records are often
## multi-allelic.
.is_haplotype_H <- function(H)
  max(H$n_alleles) > 2L || any(nchar(unlist(H$alleles, use.names = FALSE)) > 1L)

## Not exported. Typed (non-missing) individuals per record and population:
## an integer matrix, one row per record of H, one column per population.
.typed_by_pop <- function(H, pops) {
  out <- vapply(pops, function(ids) rowSums(!is.na(H$A1[, ids, drop = FALSE])),
                numeric(nrow(H$A1)))
  out <- matrix(out, nrow = nrow(H$A1), dimnames = list(NULL, names(pops)))
  storage.mode(out) <- "integer"
  out
}

## Not exported. Allele counts per record: a records x alleles matrix from
## records x individuals allele matrices `a1`/`a2` (1-based allele numbers,
## NA = missing), all in one tabulate() call. `k` is the number of allele
## columns (default: the largest allele number seen).
.allele_counts <- function(a1, a2, k = max(c(1L, a1, a2), na.rm = TRUE)) {
  nr <- nrow(a1)
  codes <- c(a1, a2)
  rows  <- rep(rep(seq_len(nr), ncol(a1)), 2L)
  ok    <- !is.na(codes)
  matrix(tabulate((codes[ok] - 1L) * nr + rows[ok], nbins = nr * k), nr, k)
}

## Not exported. H (records `rows`) as the data frame hierfstat expects: a
## population number, then one column per locus with each genotype as a
## 3-digits-per-allele integer (alleles 1 and 12 -> 1012).
.to_hierfstat_df <- function(H, pops, rows = seq_len(nrow(H$A1))) {
  ids <- unlist(pops, use.names = FALSE)
  a <- H$A1[rows, ids, drop = FALSE]; b <- H$A2[rows, ids, drop = FALSE]
  Gm <- t(pmin(a, b) * 1000L + pmax(a, b))           # individuals x loci
  dat <- data.frame(pop = rep(seq_along(pops), lengths(pops)), Gm, row.names = NULL)
  names(dat)[-1] <- paste0("L", seq_along(rows))
  dat
}

## Not exported. The text used in output filenames, e.g.
## diversity_per_population.<stem>.tsv. A caller-supplied `stem` always wins;
## otherwise it comes from the VCF's own name (populations.haps.vcf.gz ->
## "haps"), which is what lets the SNP and haplotype runs share an output
## directory without overwriting each other. Callers must require `stem`
## themselves when `vcf_file` is a list (there is no filename to use).
.derive_stem <- function(vcf_file, stem = NULL) {
  if (!is.null(stem)) return(stem)
  s <- sub("\\.gz$", "", basename(vcf_file))
  s <- sub("\\.vcf$", "", s)
  s <- sub("^.*\\.", "", s)                            # populations.haps -> haps
  if (!nzchar(s) || grepl("[^A-Za-z0-9_-]", s)) "out" else s
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
#' @examples
#' H <- read_stacks_vcf(system.file("extdata", "small.haps.vcf",
#'                                  package = "RADdiversity"), verbose = FALSE)
#' read_popmap(system.file("extdata", "small_popmap.tsv", package = "RADdiversity"),
#'             H$samples)
#' @export
read_popmap <- function(path, samples, verbose = TRUE) {
  ## colClasses = "character": without it read.delim() guesses each column's
  ## type, so all-numeric sample IDs such as "001" become the integer 1 (and
  ## no longer match the VCF's "001"), and names like "T"/"F" become
  ## TRUE/FALSE. IDs are labels, never numbers.
  pm <- utils::read.delim(path, header = FALSE, colClasses = "character",
                          col.names = c("sample", "pop"), strip.white = TRUE,
                          comment.char = "#", quote = "")
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
