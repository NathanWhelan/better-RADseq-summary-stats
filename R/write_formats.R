###############################################################################
#
#  R/write_formats.R -- write a (filtered) H out to the file formats other
#  RADseq/population-genetics tools expect.
#
#  Every function here takes `H` (the list read_haps_vcf() returns, possibly
#  already run through one or more filters from R/filter_loci.R) and writes
#  it to disk in one specific format, so a typical session looks like:
#
#    H <- read_haps_vcf("populations.snps.vcf")
#    H <- filter_call_rate(H, min_call = 0.8)
#    H <- filter_maf(H, min_maf = 0.05)
#    write_plink(H, "cleaned", pops = read_popmap("popmap.tsv", H$samples))
#
#  `pops` (a named list of sample-ID vectors, from read_popmap()) is
#  REQUIRED for the two formats that are always organized by population
#  (Genepop, FSTAT) and OPTIONAL for the rest (STRUCTURE and PLINK can both
#  carry a population label, but don't have to; plain VCF and RADpainter
#  have no notion of population at the file level at all).
#
#  Every format's exact layout below (locus-name placement, column order,
#  the missing-data code) was checked against that format's own
#  documentation or reference source code while writing this file, rather
#  than assumed from memory -- see each function's comments for what was
#  confirmed and against what.
#
###############################################################################

## Not exported. Turns `pops` (read_popmap()'s named-list-of-sample-IDs
## shape) into one population NUMBER per sample in H, in H$samples order --
## e.g. if pops = list(popA = c("a1","a2"), popB = c("b1")), then sample
## "b1" gets population number 2. Used by every writer below that needs a
## population column.
.pop_id_vector <- function(H, pops) {
  ids <- unlist(pops, use.names = FALSE)
  dupd <- unique(ids[duplicated(ids)])
  if (length(dupd))
    stop("These sample(s) appear in more than one element of `pops`, so it's ",
         "ambiguous which population they belong to: ", paste(dupd, collapse = ", "),
         ". Each sample must be assigned to exactly one population.")
  missing_samples <- setdiff(H$samples, ids)
  if (length(missing_samples))
    stop("These samples are in H but aren't assigned to any population in `pops`: ",
         paste(missing_samples, collapse = ", "),
         ". Every sample in H must appear in exactly one element of `pops` ",
         "(as read_popmap() already ensures) before this function can label them.")
  ## A named vector like c(a1 = 1, a2 = 1, b1 = 2) -- indexing it by sample
  ## NAME (not position) then gives each sample's population number, in
  ## whatever order H$samples happens to list them.
  pop_idx <- stats::setNames(rep(seq_along(pops), lengths(pops)), ids)
  unname(pop_idx[H$samples])
}

## Not exported. Genepop and FSTAT both encode each allele as a zero-padded
## number of a FIXED width for the whole file (2 digits if there are at most
## 99 possible alleles anywhere in the dataset, otherwise 3) -- picking the
## width once, for the whole file, is required by both formats' own specs
## (confirmed against the official Genepop manual and the FSTAT format
## description in hierfstat's documentation).
.digit_width <- function(H) if (max(H$n_alleles) <= 99L) 2L else 3L

## Not exported. A one-line wrapper around requireNamespace() so tests can
## force write_fstat() down its no-hierfstat fallback path (via testthat's
## mocking helpers) without actually needing hierfstat to be uninstalled --
## testthat can only override a binding that lives IN this package's own
## namespace, not one (like requireNamespace itself) inherited from base.
.hierfstat_available <- function() requireNamespace("hierfstat", quietly = TRUE)

## Not exported. Re-expresses H$A1/H$A2 (this package's small-integer allele
## numbers) as zero-padded text of a given width, e.g. allele 3 at width 2
## becomes "03". A missing allele (NA) becomes a run of zeros of that same
## width (e.g. "00") -- the standard missing-allele code both formats use.
.allele_code_matrix <- function(H, width = .digit_width(H)) {
  fmt  <- paste0("%0", width, "d")
  miss <- strrep("0", width)
  code_one <- function(A) {
    out <- sprintf(fmt, A)          # e.g. sprintf("%02d", 3) == "03"
    out[is.na(A)] <- miss
    matrix(out, nrow = nrow(A), dimnames = dimnames(A))
  }
  list(A1 = code_one(H$A1), A2 = code_one(H$A2), width = width)
}

#' Write a (filtered) H back out as a VCF file
#'
#' Serializes `H` -- including any changes made by [filter_maf()],
#' [filter_low_conf_alt()], or any other filter in this package -- back into
#' a valid VCF file, so it can be re-read with [read_haps_vcf()] or handed to
#' another VCF-based tool.
#'
#' Only the `GT` (genotype) values are ever rewritten from `H$A1`/`H$A2`
#' themselves. Everything else in each record (`QUAL`, `FILTER`, `INFO`, and
#' the original `AD`/`DP`/etc. per-sample values) is copied through
#' unchanged from the ORIGINAL file, because a filter like
#' [filter_low_conf_alt()] can change which genotype is recorded at a cell
#' without knowing (or needing to know) how to update every other per-cell
#' value that used to describe the old genotype -- writing them out anyway
#' would silently mislead anyone reading them back. For that reason, the
#' `FORMAT` column in the file this writes is always just `GT`.
#'
#' @param H A list as returned by [read_haps_vcf()]. Must still have its
#'   `fields` element (present by default; only missing if `H` was built by
#'   hand rather than from a real VCF file).
#' @param path Output file path. Ending it in `.gz` writes a gzip-compressed
#'   file, same as [read_haps_vcf()] can read back in.
#' @param verbose Print a short summary once written. Default `TRUE`.
#' @return `path`, invisibly.
#' @examples
#' vcf_lines <- c(
#'   "##fileformat=VCFv4.2",
#'   "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tind1\tind2",
#'   "un\t1\tlocus_1\tA\tC\t.\tPASS\t.\tGT\t0/0\t0/1"
#' )
#' in_file <- tempfile(fileext = ".vcf")
#' writeLines(vcf_lines, in_file)
#' H <- read_haps_vcf(in_file, verbose = FALSE)
#' out_file <- tempfile(fileext = ".vcf")
#' write_vcf(H, out_file, verbose = FALSE)
#' cat(readLines(out_file), sep = "\n")
#' @export
write_vcf <- function(H, path, verbose = TRUE) {
  if (is.null(H$fields))
    stop("write_vcf() needs H$fields (the raw VCF columns read_haps_vcf() keeps for ",
         "CHROM/POS/ID/REF/ALT/QUAL/FILTER/INFO) -- this H doesn't have it, e.g. ",
         "because it was built by hand rather than from read_haps_vcf().")
  n_rec  <- nrow(H$A1)
  n_samp <- ncol(H$A1)

  ## VCF numbers alleles starting at 0 (REF = 0, first ALT = 1, ...) while
  ## this package numbers them starting at 1 (REF = 1, first ALT = 2, ...)
  ## purely so 1-based indexing into H$alleles works naturally in R -- so
  ## converting back for the file just means subtracting 1. A missing
  ## allele (NA) becomes VCF's own missing genotype text, ".".
  a1_txt <- ifelse(is.na(H$A1), ".", as.character(H$A1 - 1L))
  a2_txt <- ifelse(is.na(H$A2), ".", as.character(H$A2 - 1L))
  gt <- matrix(paste0(a1_txt, "/", a2_txt), nrow = n_rec, dimnames = dimnames(H$A1))

  header <- c(
    "##fileformat=VCFv4.2",
    "##source=RADdiversity::write_vcf",
    paste(c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", H$samples),
          collapse = "\t")
  )

  ## Build the whole file body in one shot: a list of every column (as a
  ## plain vector, one entry per locus), then paste() glues matching
  ## positions from every column together with tabs in one vectorized pass
  ## -- much faster than looping row by row for a large VCF.
  cols <- c(
    list(H$fields[, "CHROM"], H$fields[, "POS"], H$fields[, "ID"], H$fields[, "REF"],
         H$fields[, "ALT"], H$fields[, "QUAL"], H$fields[, "FILTER"], H$fields[, "INFO"],
         rep("GT", n_rec)),
    lapply(H$samples, function(s) gt[, s])
  )
  body <- do.call(paste, c(cols, sep = "\t"))

  con <- if (grepl("\\.gz$", path)) gzfile(path, "w") else file(path, "w")
  on.exit(close(con), add = TRUE)
  writeLines(c(header, body), con)

  if (verbose)
    message(sprintf("Wrote %s (%s loci, %d samples, FORMAT=GT only -- see ?write_vcf for why)",
                    path, format(n_rec, big.mark = ","), n_samp))
  invisible(path)
}

#' Write a (filtered) H as PLINK text files
#'
#' Writes the classic PLINK 1.x text format: `<path_prefix>.map` (one line
#' per locus) and `<path_prefix>.ped` (one line per individual). Confirmed
#' against PLINK's own documented `.map`/`.ped` column layout and missing-
#' data code.
#'
#' PLINK's `.ped` format has exactly two allele columns per locus -- it has
#' no way to represent a locus with more than two alleles. Haplotype VCFs
#' from Stacks routinely do have loci with 3 or more alleles, so this
#' function checks for that using [locus_allele_stats()] before writing
#' anything.
#'
#' @param H A list as returned by [read_haps_vcf()].
#' @param path_prefix Output files are `<path_prefix>.map` and
#'   `<path_prefix>.ped`.
#' @param pops Optionally, a named list of sample-ID vectors (from
#'   [read_popmap()]) -- when given, each individual's PLINK family ID (FID)
#'   is their population name; otherwise FID is just their sample name.
#' @param drop_multiallelic If `FALSE` (the default), a locus with more than
#'   2 observed alleles causes an error (rather than writing a file PLINK
#'   can't correctly read). Set to `TRUE` to instead silently exclude just
#'   those loci and write everything else.
#' @param verbose Print a short summary once written. Default `TRUE`.
#' @return `path_prefix`, invisibly.
#' @examples
#' H <- list(
#'   A1 = matrix(c(1L, 1L, 1L, 2L), nrow = 1), A2 = matrix(c(1L, 2L, 1L, 2L), nrow = 1),
#'   locus = "locus_1", locus_raw = "locus_1", alleles = list(c("A", "C")),
#'   n_alleles = 2L, samples = c("a1", "a2", "a3", "a4")
#' )
#' out <- tempfile()
#' write_plink(H, out, verbose = FALSE)
#' cat(readLines(paste0(out, ".ped")), sep = "\n")
#' @export
write_plink <- function(H, path_prefix, pops = NULL, drop_multiallelic = FALSE, verbose = TRUE) {
  stats <- locus_allele_stats(H)
  bad <- stats$n_observed_alleles > 2L
  if (any(bad)) {
    if (!drop_multiallelic)
      stop(sum(bad), " locus/loci have more than 2 observed alleles, which PLINK's ",
           ".ped/.map format can't represent (it is strictly biallelic). Either set ",
           "drop_multiallelic = TRUE to exclude just those loci, or reduce to a ",
           "biallelic dataset first (e.g. filter_maf(), filter_thin_one_snp(), or use ",
           "a SNP rather than a haplotype VCF).")
    if (verbose)
      message(sprintf("write_plink(): excluding %s of %s multiallelic loci (drop_multiallelic = TRUE)",
                      format(sum(bad), big.mark = ","), format(nrow(stats), big.mark = ",")))
    H <- .subset_H(H, !bad)
  }
  n_rec <- nrow(H$A1)

  ## .map: chromosome, locus name, genetic distance (unused here, so always
  ## 0), base-pair position -- one line per locus. Guarded explicitly for
  ## the (unusual but possible, e.g. every locus excluded above) case of
  ## zero remaining loci: paste() mixing a zero-length vector with the
  ## constant "0" would otherwise produce ONE misleading blank-ish line
  ## instead of a genuinely empty file.
  if (n_rec == 0L) {
    map_lines <- character(0)
  } else {
    chrom <- if (!is.null(H$fields)) H$fields[, "CHROM"] else rep("0", n_rec)
    pos   <- if (!is.null(H$fields)) H$fields[, "POS"]   else as.character(seq_len(n_rec))
    map_lines <- paste(chrom, H$locus, "0", pos, sep = "\t")
  }
  writeLines(map_lines, paste0(path_prefix, ".map"))

  ## .ped: 6 standard columns (family ID, individual ID, father, mother,
  ## sex, phenotype -- PLINK's own codes for "unknown": 0, 0, 0, -9), then
  ## TWO allele columns per locus, one line per individual.
  fid <- if (!is.null(pops)) names(pops)[.pop_id_vector(H, pops)] else H$samples
  ped_prefix <- paste(fid, H$samples, "0", "0", "0", "-9")

  ## Each locus contributes its OWN two allele columns, holding the actual
  ## REF/ALT sequence text (from H$alleles) rather than this package's
  ## internal allele numbers -- PLINK's plain-text format accepts any
  ## token, including multi-letter ones. A missing allele becomes PLINK's
  ## missing code, "0".
  allele_cols <- vector("list", 2L * n_rec)
  k <- 1L
  for (j in seq_len(n_rec)) {
    al <- H$alleles[[j]]
    a1 <- al[H$A1[j, ]]; a1[is.na(H$A1[j, ])] <- "0"
    a2 <- al[H$A2[j, ]]; a2[is.na(H$A2[j, ])] <- "0"
    allele_cols[[k]] <- a1; allele_cols[[k + 1L]] <- a2
    k <- k + 2L
  }
  writeLines(do.call(paste, c(list(ped_prefix), allele_cols)), paste0(path_prefix, ".ped"))

  if (verbose)
    message(sprintf("Wrote %s.map and %s.ped (%s loci, %d individuals)",
                    path_prefix, path_prefix, format(n_rec, big.mark = ","), ncol(H$A1)))
  invisible(path_prefix)
}

#' Write a (filtered) H as a STRUCTURE input file
#'
#' Writes the standard two-rows-per-individual STRUCTURE format (one row per
#' allele copy). Confirmed against the official `structure` software
#' documentation, including its missing-data code (`-9`).
#'
#' @param H A list as returned by [read_haps_vcf()].
#' @param path Output file path.
#' @param pops Optionally, a named list of sample-ID vectors (from
#'   [read_popmap()]) giving each individual's population. Left `NULL`
#'   (the default), every individual is written with a placeholder
#'   population code of `1`.
#' @param verbose Print a short summary once written. Default `TRUE`.
#' @return `path`, invisibly.
#' @details This file always has BOTH a sample-label column and a
#'   population column, so the matching STRUCTURE run's `mainparams` should
#'   set `LABEL=1` and `POPDATA=1` (as well as `MARKERNAMES=1`, for the
#'   locus-name header row this function writes). If `pops` was left
#'   `NULL`, the population column is a constant placeholder -- `POPDATA=1`
#'   can still be set (it just won't distinguish any real groups), or the
#'   column can be deleted by hand if `POPDATA=0` is preferred instead.
#' @examples
#' H <- list(
#'   A1 = matrix(c(1L, 2L), nrow = 1), A2 = matrix(c(1L, 2L), nrow = 1),
#'   locus = "locus_1", n_alleles = 2L, samples = c("a1", "a2")
#' )
#' out <- tempfile()
#' write_structure(H, out, verbose = FALSE)
#' cat(readLines(out), sep = "\n")
#' @export
write_structure <- function(H, path, pops = NULL, verbose = TRUE) {
  n_samp <- ncol(H$A1)
  pop_id <- if (!is.null(pops)) .pop_id_vector(H, pops) else rep(1L, n_samp)
  if (is.null(pops) && verbose)
    message("write_structure(): no `pops` given -- every individual is written with a ",
            "placeholder population code of 1 (see ?write_structure for the mainparams ",
            "settings this file needs).")

  ## STRUCTURE's own documented missing-data code is -9 (this package's own
  ## NA has no meaning to STRUCTURE itself).
  recode <- function(A) { A[is.na(A)] <- -9L; A }
  A1c <- recode(H$A1); A2c <- recode(H$A2)

  header <- paste(H$locus, collapse = "\t")
  lines <- character(2L * n_samp)
  for (i in seq_len(n_samp)) {
    ## Two rows per individual -- one for each allele copy -- both starting
    ## with that individual's name and population code.
    lines[2L * i - 1L] <- paste(c(H$samples[i], pop_id[i], A1c[, i]), collapse = "\t")
    lines[2L * i]      <- paste(c(H$samples[i], pop_id[i], A2c[, i]), collapse = "\t")
  }
  writeLines(c(header, lines), path)
  if (verbose)
    message(sprintf("Wrote %s (%s loci, %d individuals, 2 rows each)",
                    path, format(nrow(H$A1), big.mark = ","), n_samp))
  invisible(path)
}

#' Write a (filtered) H as a Genepop input file
#'
#' Writes the standard Genepop text layout: a title line, one locus-name
#' line per locus, then each population introduced by a `POP` line followed
#' by one line per individual. Confirmed against the official Genepop
#' manual, including its missing-data code (all-zero allele codes, e.g.
#' `0000`) and comma-after-sample-ID convention.
#'
#' @param H A list as returned by [read_haps_vcf()].
#' @param path Output file path.
#' @param pops A named list of sample-ID vectors (from [read_popmap()]).
#'   Required -- Genepop files are always organized into population blocks,
#'   so unlike [write_structure()]/[write_plink()] there's no sensible
#'   single-population default.
#' @param title Text for the file's first line. Default
#'   `"RADdiversity export"`.
#' @param verbose Print a short summary once written. Default `TRUE`.
#' @return `path`, invisibly.
#' @examples
#' # Column (sample) names on A1/A2 matter here -- write_genepop() looks
#' # up each population's samples by name, exactly as read_haps_vcf()'s
#' # real output always allows.
#' H <- list(
#'   A1 = matrix(c(1L, 2L), nrow = 1, dimnames = list(NULL, c("a1", "a2"))),
#'   A2 = matrix(c(1L, 2L), nrow = 1, dimnames = list(NULL, c("a1", "a2"))),
#'   locus = "locus_1", n_alleles = 2L, samples = c("a1", "a2")
#' )
#' pops <- list(popA = "a1", popB = "a2")
#' out <- tempfile()
#' write_genepop(H, out, pops = pops, verbose = FALSE)
#' cat(readLines(out), sep = "\n")
#' @export
write_genepop <- function(H, path, pops, title = "RADdiversity export", verbose = TRUE) {
  if (is.null(pops))
    stop("write_genepop() needs `pops` (a named list of sample IDs per population, e.g. ",
         "from read_popmap()) -- Genepop files are always organized into population ",
         "blocks; there is no sensible single-population default the way there is for ",
         "write_structure()/write_plink().")
  width <- .digit_width(H)
  codes <- .allele_code_matrix(H, width)

  ## Combine each individual's two allele codes into ONE genotype code per
  ## locus, e.g. width 2, alleles 1 and 3 -> "0103". The smaller code is
  ## always written first so the same underlying genotype is always written
  ## the same way, regardless of which allele happened to be read first.
  ## (matrix() re-wraps pmin()/pmax()'s result with the original row/column
  ## layout explicitly, rather than relying on them to preserve it.)
  lo <- matrix(pmin(codes$A1, codes$A2), nrow = nrow(codes$A1), dimnames = dimnames(codes$A1))
  hi <- matrix(pmax(codes$A1, codes$A2), nrow = nrow(codes$A1), dimnames = dimnames(codes$A1))
  geno <- matrix(paste0(lo, hi), nrow = nrow(codes$A1), dimnames = dimnames(codes$A1))
  geno[is.na(H$A1) | is.na(H$A2)] <- strrep("0", 2L * width)  # Genepop's missing code

  lines <- c(title, H$locus)
  for (p in names(pops)) {
    lines <- c(lines, "POP")
    for (s in pops[[p]])
      lines <- c(lines, paste0(s, " ,\t", paste(geno[, s], collapse = "\t")))
  }
  writeLines(lines, path)
  if (verbose)
    message(sprintf("Wrote %s (%s loci, %d populations, %d individuals, %d-digit allele codes)",
                    path, format(nrow(H$A1), big.mark = ","), length(pops), ncol(H$A1), width))
  invisible(path)
}

#' Write a (filtered) H as an FSTAT input file
#'
#' Writes the FSTAT text format used by FSTAT itself and by the `hierfstat`
#' R package. When `hierfstat` is installed, this delegates to its own
#' `write.fstat()` (reusing the exact genotype-coding data frame shape
#' [diversity_stats()] already builds internally); otherwise an internal
#' fallback writes the format directly. Both confirmed against `hierfstat`'s
#' documented FSTAT format description.
#'
#' @param H A list as returned by [read_haps_vcf()].
#' @param path Output file path.
#' @param pops A named list of sample-ID vectors (from [read_popmap()]).
#'   Required -- FSTAT's first data column is always the population number.
#' @param verbose Print a short summary once written. Default `TRUE`.
#' @return `path`, invisibly.
#' @details When `hierfstat` is installed, this relies on its own
#'   `write.fstat()`, which has a known limitation of its own with a
#'   dataset of exactly ONE locus (it errors rather than writes a
#'   one-column file) -- not something a real RADseq dataset (always many
#'   loci) will ever run into, but worth knowing if you're experimenting
#'   with a tiny hand-built `H` for testing.
#' @examples
#' H <- list(
#'   A1 = matrix(c(1L, 2L, 1L, 1L), nrow = 2),
#'   A2 = matrix(c(1L, 2L, 1L, 2L), nrow = 2),
#'   locus = c("locus_1", "locus_2"), locus_raw = c("locus_1", "locus_2"),
#'   n_alleles = c(2L, 2L), samples = c("a1", "a2")
#' )
#' pops <- list(popA = "a1", popB = "a2")
#' out <- tempfile()
#' write_fstat(H, out, pops = pops, verbose = FALSE)
#' cat(readLines(out), sep = "\n")
#' @export
write_fstat <- function(H, path, pops, verbose = TRUE) {
  if (is.null(pops))
    stop("write_fstat() needs `pops` (see ?write_genepop for why) -- FSTAT's first ",
         "data column is always the population number.")
  pop_id <- .pop_id_vector(H, pops)
  n_rec  <- nrow(H$A1)
  n_samp <- ncol(H$A1)

  if (.hierfstat_available()) {
    ## Reuse hierfstat's own writer, fed the exact same genotype-coding
    ## scheme diversity_stats() already builds internally for its hierfstat
    ## calculations (pmin/pmax of the two allele numbers, glued together as
    ## a 3-digit-per-allele integer, e.g. alleles 1 and 12 -> 1012) --
    ## reusing an already-checked pattern rather than inventing a second one.
    Gm <- matrix(NA_integer_, n_samp, n_rec)
    for (j in seq_len(n_rec)) {
      a <- H$A1[j, ]; b <- H$A2[j, ]
      Gm[, j] <- pmin(a, b) * 1000L + pmax(a, b)
    }
    dat <- data.frame(pop = pop_id, Gm)
    names(dat)[-1] <- H$locus
    hierfstat::write.fstat(dat, fname = path)
  } else {
    ## No hierfstat installed -- write the FSTAT text format directly:
    ##   line 1: population count, locus count, largest allele count seen,
    ##           digit width
    ##   next nl lines: one locus name each
    ##   remaining lines: one per individual, population number then one
    ##           genotype code per locus (both alleles glued together)
    width <- .digit_width(H)
    codes <- .allele_code_matrix(H, width)
    lo <- matrix(pmin(codes$A1, codes$A2), nrow = n_rec, dimnames = dimnames(codes$A1))
    hi <- matrix(pmax(codes$A1, codes$A2), nrow = n_rec, dimnames = dimnames(codes$A1))
    geno <- matrix(paste0(lo, hi), nrow = n_rec, dimnames = dimnames(codes$A1))
    geno[is.na(H$A1) | is.na(H$A2)] <- strrep("0", 2L * width)
    nal <- max(locus_allele_stats(H)$n_observed_alleles, na.rm = TRUE)

    data_lines <- vapply(seq_len(n_samp), function(i)
      paste(pop_id[i], paste(geno[, i], collapse = " ")), character(1))
    writeLines(c(sprintf("%d %d %d %d", length(pops), n_rec, nal, width), H$locus, data_lines), path)
  }

  if (verbose)
    message(sprintf("Wrote %s (%s loci, %d populations, %d individuals)%s",
                    path, format(n_rec, big.mark = ","), length(pops), n_samp,
                    if (.hierfstat_available()) " via hierfstat::write.fstat()" else ""))
  invisible(path)
}

#' Write a (filtered) H as a RADpainter/fineRADstructure input file
#'
#' Writes fineRADstructure's haplotype input format for its `RADpainter`
#' program: a header line of sample names, then one line per RAD locus
#' giving each individual's two haplotype alleles. Confirmed directly
#' against two independent reference implementations of this format --
#' `hapsFromVCF.cpp` (the fineRADstructure project's own converter) and the
#' community `finerad_input.py` script -- which agree exactly on every
#' detail, including the easy-to-miss ones: no leading blank column before
#' the sample names, and a missing genotype is written as a completely
#' EMPTY field (nothing between its two surrounding tabs), not `-9`, `?`,
#' or any other placeholder text.
#'
#' RADpainter is a haplotype-based method, so this is only meaningful for a
#' haplotype-type `H` (one multi-allelic record per RAD tag, e.g. from
#' `populations.haps.vcf`) -- not a plain per-site SNP VCF.
#'
#' @param H A list as returned by [read_haps_vcf()], ideally from a
#'   haplotype VCF.
#' @param path Output file path.
#' @param verbose Print a short summary once written. Default `TRUE`.
#' @return `path`, invisibly.
#' @examples
#' H <- list(
#'   A1 = matrix(c(1L, 2L), nrow = 1), A2 = matrix(c(1L, 2L), nrow = 1),
#'   locus = "locus_1", alleles = list(c("AACGT", "AACGG")),
#'   n_alleles = 2L, samples = c("a1", "a2")
#' )
#' out <- tempfile()
#' write_radpainter(H, out, verbose = FALSE)
#' cat(readLines(out), sep = "\n")
#' @export
write_radpainter <- function(H, path, verbose = TRUE) {
  ## Same heuristic diversity_stats() already uses to tell a haplotype VCF
  ## from a SNP VCF: haplotype alleles are usually multi-letter sequences,
  ## and/or a locus commonly has more than 2 of them. If NEITHER is true
  ## anywhere in this dataset, it's very likely per-site SNP data instead --
  ## still writeable, just probably not what RADpainter is meant to analyze.
  looks_snp <- max(H$n_alleles) <= 2L &&
    all(nchar(unlist(H$alleles, use.names = FALSE)) == 1L)
  if (looks_snp && verbose)
    message("write_radpainter(): every locus here looks biallelic with single-",
            "nucleotide alleles, i.e. this looks like SNP data rather than RAD-tag ",
            "haplotypes. RADpainter is meant to run on haplotypes (e.g. from ",
            "populations.haps.vcf) -- double check this is the file you meant to use.")

  n_rec  <- nrow(H$A1)
  header <- paste(H$samples, collapse = "\t")

  lines <- character(n_rec)
  for (j in seq_len(n_rec)) {
    al <- H$alleles[[j]]
    a1 <- al[H$A1[j, ]]; a2 <- al[H$A2[j, ]]
    ## Each individual's cell is its two haplotype alleles joined by "/";
    ## a missing genotype becomes a totally empty cell, per the confirmed
    ## reference format (see the function's help page for the sources).
    cell <- paste0(a1, "/", a2)
    cell[is.na(H$A1[j, ]) | is.na(H$A2[j, ])] <- ""
    lines[j] <- paste(cell, collapse = "\t")
  }
  writeLines(c(header, lines), path)
  if (verbose)
    message(sprintf("Wrote %s (%s RAD loci, %d samples)", path, format(n_rec, big.mark = ","), ncol(H$A1)))
  invisible(path)
}
