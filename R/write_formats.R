###############################################################################
#
#  R/write_formats.R -- write a (filtered) H to the file formats other
#  population-genetics programs read.
#
#  Every function takes `H` (the object read_stacks_vcf() returns, possibly
#  run through filter_*() functions) and writes one format:
#
#    H <- read_stacks_vcf("populations.snps.vcf")
#    H <- filter_call_rate(H, min_call = 0.8)
#    write_plink(H, "cleaned", popmap = "popmap.tsv")
#
#  `popmap` (a popmap file path or the list from read_popmap()) is REQUIRED
#  for the formats organized by population (Genepop, FSTAT) and OPTIONAL for
#  STRUCTURE and PLINK, which can carry a population label. VCF and
#  RADpainter files have no population field.
#
#  Each format's layout (locus-name placement, column order, missing-data
#  code) follows that format's documentation or reference implementation, as
#  noted at each function. Checked by tests/testthat/test-write-formats.R.
#
###############################################################################

## Not exported. For each sample of H (in H$samples order), the NUMBER of its
## population in `pops`: with pops = list(popA = c("a1","a2"), popB = "b1"),
## sample "b1" gets 2. Stops if a sample of H has no population.
.pop_id_vector <- function(H, pops) {
  ids <- unlist(pops, use.names = FALSE)
  unassigned <- setdiff(H$samples, ids)
  if (length(unassigned))
    stop("These samples are in H but not in `popmap`: ", paste(unassigned, collapse = ", "),
         ". Every sample must be assigned to a population; remove the others from H first ",
         "with filter_samples(H, popmap).",
         call. = FALSE)
  pop_number <- stats::setNames(rep(seq_along(pops), lengths(pops)), ids)
  unname(pop_number[H$samples])
}

## Not exported. `popmap` for the writers: required for Genepop and FSTAT,
## optional for the others. Returns the population list, or NULL when
## `popmap` is NULL and not required.
.writer_pops <- function(H, popmap, required, format_name) {
  if (is.null(popmap)) {
    if (required)
      stop(format_name, " needs `popmap` (a popmap file path, or the list returned by ",
           "read_popmap()): ", format_name, " files are organized by population.",
           call. = FALSE)
    return(NULL)
  }
  .resolve_pops(popmap, H$samples, verbose = FALSE)
}

## Not exported. Stops, naming them, if any of `names` (sample or population
## names) contains `pattern`, a character that would break `format_name`'s file
## layout; `why` says what that character does there.
.stop_on_unwritable_names <- function(names, pattern, format_name, why) {
  bad <- unique(names[grepl(pattern, names)])
  if (length(bad))
    stop(format_name, " ", why, ", so these names can't be written: ",
         paste0("\"", bad, "\"", collapse = ", "),
         ". Rename them (in the VCF header or popmap) and re-run.", call. = FALSE)
  invisible(NULL)
}

## Not exported. Genepop and FSTAT write each allele as a zero-padded number
## of one FIXED width for the whole file: 2 digits when no record has more
## than 99 alleles, otherwise 3 (Genepop manual; FSTAT format description in
## hierfstat).
.digit_width <- function(H) if (max(H$n_alleles) <= 99L) 2L else 3L

## Not exported. Each genotype as the text Genepop and FSTAT use: both allele
## numbers, zero-padded to `width` and smaller first (alleles 3 and 1 at width
## 2 -> "0103"), with all zeros ("0000") for a missing genotype. Returns a
## records x samples character matrix.
.genotype_codes <- function(H, width = .digit_width(H)) {
  smaller <- pmin(H$A1, H$A2)
  larger <- pmax(H$A1, H$A2)
  codes <- matrix(paste0(sprintf(paste0("%0", width, "d"), smaller),
                         sprintf(paste0("%0", width, "d"), larger)),
                  nrow = nrow(H$A1), dimnames = dimnames(H$A1))
  codes[is.na(H$A1) | is.na(H$A2)] <- strrep("0", 2L * width)
  codes
}

#' Write a (filtered) H back out as a VCF file
#'
#' Writes `H`, including any changes made by the `filter_*()` functions, as a
#' VCF file that [read_stacks_vcf()] or other VCF tools can read.
#'
#' Only the genotypes (`GT`) are written from `H$A1`/`H$A2`. `QUAL`, `FILTER`
#' and `INFO` are copied from the original file, but per-sample fields such as
#' `AD` and `DP` are NOT: a filter such as [filter_low_conf_alt()] can change a
#' genotype without updating those values, and writing them would mislead
#' anyone reading the file. The `FORMAT` column is therefore always `GT`.
#'
#' The original file's `##` meta-information lines (`##contig`, `##INFO`,
#' `##FORMAT`, ...) are not kept: [read_stacks_vcf()] does not store them. Some
#' tools (e.g. `bcftools`) warn about INFO keys with no `##INFO` definition;
#' copy the header lines from the original file if a tool needs them.
#'
#' @references
#' Danecek, P., Auton, A., Abecasis, G., et al. (2011) The variant call format
#' and VCFtools. *Bioinformatics* 27:2156-2158.
#'
#' @param H The object returned by [read_stacks_vcf()] (optionally filtered).
#'   It must still contain `$fields`, the original VCF columns.
#' @param path Output file path. A name ending in `.gz` writes a
#'   gzip-compressed file.
#' @param verbose Print a short summary once written. Default `TRUE`.
#' @return `path`, invisibly.
#' @examples
#' H <- read_stacks_vcf(system.file("extdata", "small.haps.vcf",
#'                                  package = "RADdiversity"), verbose = FALSE)
#' out_file <- tempfile(fileext = ".vcf")
#' write_vcf(H, out_file)
#' head(readLines(out_file), 4)
#' @export
write_vcf <- function(H, path, verbose = TRUE) {
  H <- .resolve_H(H, verbose = FALSE)
  .check_string(path, "path")
  .check_flag(verbose, "verbose")
  if (is.null(H$fields))
    stop("write_vcf() needs H$fields (the original VCF columns CHROM ... INFO), but this ",
         "H has none -- for example because it was built by hand rather than read with ",
         "read_stacks_vcf().", call. = FALSE)
  n_rec <- nrow(H$A1)

  ## VCF numbers alleles from 0 (REF = 0); this package numbers them from 1,
  ## so writing subtracts 1. A missing allele is written as ".".
  allele_text <- function(A) ifelse(is.na(A), ".", as.character(A - 1L))
  gt <- matrix(paste0(allele_text(H$A1), "/", allele_text(H$A2)), nrow = n_rec,
               dimnames = dimnames(H$A1))

  header <- c("##fileformat=VCFv4.2",
              "##source=RADdiversity::write_vcf",
              paste(c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT",
                      H$samples), collapse = "\t"))
  ## One vector per output column, then paste() joins them row by row with
  ## tabs in a single vectorised call.
  columns <- c(lapply(c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"),
                      function(name) H$fields[, name]),
               list(rep("GT", n_rec)),
               lapply(H$samples, function(s) gt[, s]))
  body <- do.call(paste, c(columns, sep = "\t"))

  con <- if (grepl("\\.gz$", path)) gzfile(path, "w") else file(path, "w")
  on.exit(close(con), add = TRUE)
  writeLines(c(header, body), con)
  .inform(verbose, sprintf("Wrote %s (%s records, %d samples, FORMAT=GT only -- see ?write_vcf)",
                           path, .big(n_rec), ncol(H$A1)))
  invisible(path)
}

#' Write a (filtered) H as PLINK text files
#'
#' Writes the PLINK 1 text format: `<path_prefix>.map` (one line per record)
#' and `<path_prefix>.ped` (one line per individual), with PLINK's documented
#' column layout and missing-data code (`0`).
#'
#' PLINK's `.ped` format has exactly two allele columns per locus, so it
#' cannot represent a record with more than two alleles, which haplotype VCFs
#' often have. Such records stop the function unless `drop_multiallelic =
#' TRUE`.
#'
#' `.ped` columns are separated by white space, so a sample or population
#' name containing a space stops the function. The `.map` chromosome column
#' is the VCF's `CHROM`. Stacks de novo output uses names such as `un` or
#' locus numbers, which PLINK 1.9 and 2 reject unless run with
#' `--allow-extra-chr`.
#'
#' @references
#' Purcell, S., Neale, B., Todd-Brown, K., et al. (2007) PLINK: a tool set for
#' whole-genome association and population-based linkage analyses. *American
#' Journal of Human Genetics* 81:559-575.
#'
#' Chang, C.C., Chow, C.C., Tellier, L.C.A.M., Vattikuti, S., Purcell, S.M. &
#' Lee, J.J. (2015) Second-generation PLINK: rising to the challenge of larger
#' and richer datasets. *GigaScience* 4:7.
#'
#' @inheritParams write_vcf
#' @param path_prefix Output files are `<path_prefix>.map` and
#'   `<path_prefix>.ped`.
#' @param popmap Optional: a popmap file path or the list returned by
#'   [read_popmap()]. When given, each individual's PLINK family ID (FID) is
#'   its population name; otherwise FID is the sample name.
#' @param drop_multiallelic If `FALSE` (default), a record with more than 2
#'   observed alleles stops the function. If `TRUE`, such records are left
#'   out (with a message).
#' @return `path_prefix`, invisibly.
#' @examples
#' H <- read_stacks_vcf(system.file("extdata", "small.snps.vcf",
#'                                  package = "RADdiversity"), verbose = FALSE)
#' popmap <- system.file("extdata", "small_popmap.tsv", package = "RADdiversity")
#' out <- tempfile()
#' write_plink(H, out, popmap = popmap)
#' readLines(paste0(out, ".map"))[1:3]
#' @export
write_plink <- function(H, path_prefix, popmap = NULL, drop_multiallelic = FALSE, verbose = TRUE) {
  H <- .resolve_H(H, verbose = FALSE)
  .check_string(path_prefix, "path_prefix")
  .check_flag(drop_multiallelic, "drop_multiallelic")
  .check_flag(verbose, "verbose")
  pops <- .writer_pops(H, popmap, required = FALSE, "PLINK")
  ## .ped fields are whitespace-separated: a space inside a name would shift
  ## every later column.
  .stop_on_unwritable_names(c(H$samples, names(pops)), "[[:space:]]", "PLINK .ped files",
                            "separate columns by white space")

  multiallelic <- locus_allele_stats(H)$n_observed_alleles > 2L
  if (any(multiallelic)) {
    if (!drop_multiallelic)
      stop(sum(multiallelic), " record(s) have more than 2 observed alleles, which PLINK's ",
           ".ped/.map format can't represent (it is strictly biallelic). Set ",
           "drop_multiallelic = TRUE to leave them out, or use a biallelic dataset (e.g. ",
           "the SNP VCF rather than the haplotype VCF).", call. = FALSE)
    .inform(verbose, sprintf("write_plink(): excluding %s of %s multiallelic records (drop_multiallelic = TRUE)",
                             .big(sum(multiallelic)), .big(length(multiallelic))))
    H <- .subset_H(H, !multiallelic)
  }
  n_rec <- nrow(H$A1)

  ## .map: chromosome, marker name, genetic distance (unknown: 0), base-pair
  ## position. With no records left the file is empty (paste() would
  ## otherwise make one line out of the constant "0").
  map_lines <- character(0)
  if (n_rec > 0L) {
    chrom <- if (!is.null(H$fields)) H$fields[, "CHROM"] else rep("0", n_rec)
    pos <- if (!is.null(H$fields)) H$fields[, "POS"] else as.character(seq_len(n_rec))
    map_lines <- paste(chrom, H$locus, "0", pos, sep = "\t")
  }
  writeLines(map_lines, paste0(path_prefix, ".map"))

  ## .ped: family ID, individual ID, father, mother, sex, phenotype (PLINK's
  ## "unknown" codes 0, 0, 0, -9), then two allele columns per record holding
  ## the allele sequences; "0" = missing.
  family_id <- if (!is.null(pops)) names(pops)[.pop_id_vector(H, pops)] else H$samples
  ped_start <- paste(family_id, H$samples, "0", "0", "0", "-9")
  allele_columns <- vector("list", 2L * n_rec)
  for (j in seq_len(n_rec)) {
    sequences <- H$alleles[[j]]
    first <- sequences[H$A1[j, ]]
    second <- sequences[H$A2[j, ]]
    first[is.na(H$A1[j, ])] <- "0"
    second[is.na(H$A2[j, ])] <- "0"
    allele_columns[[2L * j - 1L]] <- first
    allele_columns[[2L * j]] <- second
  }
  writeLines(do.call(paste, c(list(ped_start), allele_columns)), paste0(path_prefix, ".ped"))

  .inform(verbose, sprintf("Wrote %s.map and %s.ped (%s records, %d individuals)",
                           path_prefix, path_prefix, .big(n_rec), ncol(H$A1)))
  invisible(path_prefix)
}

#' Write a (filtered) H as a STRUCTURE input file
#'
#' Writes the two-rows-per-individual STRUCTURE format (one row per allele
#' copy), with STRUCTURE's missing-data code `-9`.
#'
#' @details The file has a locus-name row, a sample-label column and a
#'   population column, so the STRUCTURE run's `mainparams` should set
#'   `MARKERNAMES=1`, `LABEL=1` and `POPDATA=1`. Without `popmap` the
#'   population column is a placeholder `1` for everyone. STRUCTURE separates
#'   columns by white space, so a sample name containing a space stops the
#'   function.
#' @references
#' Pritchard, J.K., Stephens, M. & Donnelly, P. (2000) Inference of population
#' structure using multilocus genotype data. *Genetics* 155:945-959.
#' @inheritParams write_vcf
#' @param popmap Optional: a popmap file path or the list returned by
#'   [read_popmap()]. Default `NULL`: every individual gets population `1`.
#' @return `path`, invisibly.
#' @examples
#' H <- read_stacks_vcf(system.file("extdata", "small.snps.vcf",
#'                                  package = "RADdiversity"), verbose = FALSE)
#' popmap <- system.file("extdata", "small_popmap.tsv", package = "RADdiversity")
#' out <- tempfile()
#' write_structure(H, out, popmap = popmap)
#' @export
write_structure <- function(H, path, popmap = NULL, verbose = TRUE) {
  H <- .resolve_H(H, verbose = FALSE)
  .check_string(path, "path")
  .check_flag(verbose, "verbose")
  pops <- .writer_pops(H, popmap, required = FALSE, "STRUCTURE")
  ## STRUCTURE reads columns separated by white space: a space inside a sample
  ## name would shift every genotype after it.
  .stop_on_unwritable_names(H$samples, "[[:space:]]", "STRUCTURE files",
                            "separate columns by white space")
  n_samp <- ncol(H$A1)
  pop_id <- if (!is.null(pops)) .pop_id_vector(H, pops) else rep(1L, n_samp)
  if (is.null(pops))
    .inform(verbose, "write_structure(): no `popmap` given -- every individual is written with a ",
            "placeholder population code of 1 (see ?write_structure for the mainparams settings).")

  missing_as_minus_9 <- function(A) {
    A[is.na(A)] <- -9L
    A
  }
  A1 <- missing_as_minus_9(H$A1)
  A2 <- missing_as_minus_9(H$A2)
  lines <- character(2L * n_samp)
  for (i in seq_len(n_samp)) {
    ## Two rows per individual, one per allele copy, each starting with the
    ## individual's name and population code.
    lines[2L * i - 1L] <- paste(c(H$samples[i], pop_id[i], A1[, i]), collapse = "\t")
    lines[2L * i] <- paste(c(H$samples[i], pop_id[i], A2[, i]), collapse = "\t")
  }
  writeLines(c(paste(H$locus, collapse = "\t"), lines), path)
  .inform(verbose, sprintf("Wrote %s (%s records, %d individuals, 2 rows each)",
                           path, .big(nrow(H$A1)), n_samp))
  invisible(path)
}

#' Write a (filtered) H as a Genepop input file
#'
#' Writes the Genepop text layout: a title line, one locus-name line per
#' record, then for each population a `POP` line followed by one line per
#' individual (`name ,` then its genotypes). Missing genotypes are all zeros
#' (e.g. `0000`), as in the Genepop manual. Genepop ends a name at the first
#' comma, so a sample name containing a comma stops the function.
#'
#' @references
#' Rousset, F. (2008) GENEPOP'007: a complete re-implementation of the GENEPOP
#' software for Windows and Linux. *Molecular Ecology Resources* 8:103-106.
#'
#' @inheritParams write_vcf
#' @param popmap A popmap file path or the list returned by [read_popmap()].
#'   Required: Genepop files are organized by population.
#' @param title Text for the file's first line. Default
#'   `"RADdiversity export"`.
#' @return `path`, invisibly.
#' @examples
#' H <- read_stacks_vcf(system.file("extdata", "small.haps.vcf",
#'                                  package = "RADdiversity"), verbose = FALSE)
#' popmap <- system.file("extdata", "small_popmap.tsv", package = "RADdiversity")
#' out <- tempfile()
#' write_genepop(H, out, popmap = popmap)
#' readLines(out)[1:3]
#' @export
write_genepop <- function(H, path, popmap = NULL, title = "RADdiversity export", verbose = TRUE) {
  H <- .resolve_H(H, verbose = FALSE)
  .check_string(path, "path")
  .check_string(title, "title")
  .check_flag(verbose, "verbose")
  pops <- .writer_pops(H, popmap, required = TRUE, "Genepop")
  .pop_id_vector(H, pops)                 # stops if a sample has no population
  ## Genepop ends an individual's name at the first comma (spaces are allowed).
  .stop_on_unwritable_names(H$samples, ",", "Genepop files",
                            "end each individual's name at a comma")
  width <- .digit_width(H)
  codes <- .genotype_codes(H, width)

  ## One POP block per population; each individual's line is its name, " ,",
  ## and its genotype codes.
  individual_lines <- lapply(names(pops), function(p) {
    c("POP", paste0(pops[[p]], " ,\t",
                    apply(codes[, pops[[p]], drop = FALSE], 2L, paste, collapse = "\t")))
  })
  writeLines(c(title, H$locus, unlist(individual_lines)), path)
  .inform(verbose, sprintf("Wrote %s (%s records, %d populations, %d individuals, %d-digit allele codes)",
                           path, .big(nrow(H$A1)), length(pops), ncol(H$A1), width))
  invisible(path)
}

#' Write a (filtered) H as an FSTAT input file
#'
#' Writes the FSTAT text format read by FSTAT and by the `hierfstat` R
#' package: a header line (number of populations, number of loci, largest
#' number of alleles, digits per allele), one locus name per line, then one
#' line per individual with its population number and genotype codes.
#'
#' @references
#' Goudet, J. (1995) FSTAT (version 1.2): a computer program to calculate
#' F-statistics. *Journal of Heredity* 86:485-486.
#'
#' Goudet, J. (2005) HIERFSTAT, a package for R to compute and test
#' hierarchical F-statistics. *Molecular Ecology Notes* 5:184-186.
#'
#' @inheritParams write_genepop
#' @param popmap A popmap file path or the list returned by [read_popmap()].
#'   Required: FSTAT's first data column is the population number.
#' @return `path`, invisibly.
#' @examples
#' H <- read_stacks_vcf(system.file("extdata", "small.haps.vcf",
#'                                  package = "RADdiversity"), verbose = FALSE)
#' popmap <- system.file("extdata", "small_popmap.tsv", package = "RADdiversity")
#' out <- tempfile()
#' write_fstat(H, out, popmap = popmap)
#' readLines(out)[1:3]
#' @export
write_fstat <- function(H, path, popmap = NULL, verbose = TRUE) {
  H <- .resolve_H(H, verbose = FALSE)
  .check_string(path, "path")
  .check_flag(verbose, "verbose")
  pops <- .writer_pops(H, popmap, required = TRUE, "FSTAT")
  pop_id <- .pop_id_vector(H, pops)
  n_rec <- nrow(H$A1)
  width <- .digit_width(H)
  codes <- .genotype_codes(H, width)
  ## FSTAT's third header number is the HIGHEST allele number used in the
  ## file, not how many alleles were observed: with alleles 2 and 3 present
  ## but not 1, it must be 3.
  max_alleles <- max(c(1L, H$A1, H$A2), na.rm = TRUE)
  data_lines <- vapply(seq_len(ncol(H$A1)), function(i)
    paste(pop_id[i], paste(codes[, i], collapse = " ")), character(1))
  writeLines(c(sprintf("%d %d %d %d", length(pops), n_rec, max_alleles, width),
               H$locus, data_lines), path)
  .inform(verbose, sprintf("Wrote %s (%s records, %d populations, %d individuals)",
                           path, .big(n_rec), length(pops), ncol(H$A1)))
  invisible(path)
}

#' Write a (filtered) H as a RADpainter / fineRADstructure input file
#'
#' Writes the haplotype input format of fineRADstructure's `RADpainter`: a
#' header line of sample names, then one line per RAD locus with each
#' individual's two haplotypes. The layout follows the fineRADstructure
#' project's `hapsFromVCF.cpp` and the community `finerad_input.py` script:
#' no blank column before the sample names, and a missing genotype written
#' as an EMPTY field (nothing between its two tabs).
#'
#' RADpainter uses haplotypes, so this is meant for a haplotype `H` (one
#' multi-allelic record per RAD tag, from `populations.haps.vcf`).
#'
#' @references
#' Malinsky, M., Trucchi, E., Lawson, D.J. & Falush, D. (2018) RADpainter and
#' fineRADstructure: population inference from RADseq data. *Molecular Biology
#' and Evolution* 35:1284-1290.
#'
#' @inheritParams write_vcf
#' @param H The object returned by [read_stacks_vcf()], ideally from a
#'   haplotype VCF.
#' @return `path`, invisibly.
#' @examples
#' H <- read_stacks_vcf(system.file("extdata", "small.haps.vcf",
#'                                  package = "RADdiversity"), verbose = FALSE)
#' out <- tempfile()
#' write_radpainter(H, out)
#' readLines(out)[1:2]
#' @export
write_radpainter <- function(H, path, verbose = TRUE) {
  H <- .resolve_H(H, verbose = FALSE)
  .check_string(path, "path")
  .check_flag(verbose, "verbose")
  if (!.is_haplotype_H(H))
    .inform(verbose, "write_radpainter(): every allele here is a single base, i.e. this ",
            "looks like SNP data rather than RAD-tag ",
            "haplotypes. RADpainter is meant for haplotypes (e.g. from ",
            "populations.haps.vcf) -- check this is the file you meant to use.")
  n_rec <- nrow(H$A1)
  lines <- character(n_rec)
  for (j in seq_len(n_rec)) {
    sequences <- H$alleles[[j]]
    ## An individual's cell is its two haplotypes joined by "/"; a missing
    ## genotype is an empty cell.
    cell <- paste0(sequences[H$A1[j, ]], "/", sequences[H$A2[j, ]])
    cell[is.na(H$A1[j, ]) | is.na(H$A2[j, ])] <- ""
    lines[j] <- paste(cell, collapse = "\t")
  }
  writeLines(c(paste(H$samples, collapse = "\t"), lines), path)
  .inform(verbose, sprintf("Wrote %s (%s RAD loci, %d samples)", path, .big(n_rec), ncol(H$A1)))
  invisible(path)
}
