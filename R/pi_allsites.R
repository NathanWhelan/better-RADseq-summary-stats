###############################################################################
#
#  R/pi_allsites.R -- nucleotide diversity (pi), divergence (dxy) and
#  per-individual heterozygosity per sequenced site, from an ALL-SITES VCF.
#
#  WHY. diversity_stats() gets per-site values by rescaling He from a
#  variant-only VCF by the number of sequenced sites (`sites`). That is an
#  approximation: it cannot know how many individuals were typed at each
#  INVARIANT site, so missing data at invariant sites are treated as if they
#  were not missing. An all-sites VCF -- Stacks' `populations --vcf-all`
#  (Stacks >= 2.62) -- records every sequenced site, variable or not, with
#  its own genotypes, and pi can then be computed the way pixy does it
#  (Korunes & Samuk 2021): each site contributes its own number of pairwise
#  differences and of pairwise comparisons among the gene copies actually
#  typed there, and
#
#       pi  = sum over sites of differences / sum over sites of comparisons
#
#  For one population at one site with allele counts c_a among m typed gene
#  copies: comparisons = C(m, 2), differences = C(m, 2) - sum_a C(c_a, 2).
#  Between populations x and y: comparisons = m_x m_y, differences =
#  m_x m_y - sum_a c_xa c_ya (dxy). With no missing data, pi is the average
#  over all sites of (m/(m-1))(1 - sum p^2) -- gene_div_2n_counts() -- i.e.
#  exactly He_2n_corr per sequenced site. Da = dxy - (pi_x + pi_y)/2 is Nei's
#  net divergence.
#
#  pi_nc. That pixy estimator counts the pair of gene copies inside each
#  individual as one of its comparisons, which is why it is the estimator the
#  rest of this package calls biased when FIS != 0 (R/estimators.R, section
#  1): inbreeding makes those two copies alike. Leaving them out -- for each
#  typed individual, one comparison fewer, and one difference fewer if it is
#  heterozygous -- gives the diversity between copies from DIFFERENT
#  individuals, which does not depend on FIS. At one site that is exactly
#  Nei & Chesser's Hs (hs_from_counts()):
#
#    [C(2n,2)(2n/(2n-1))(1 - sum p^2) - n Ho] / [C(2n,2) - n]
#        = (n/(n-1)) (1 - sum p^2 - Ho/(2n))
#
#  so pi_nc per sequenced site is He (Nei & Chesser) per sequenced site, with
#  missing data handled site by site. `pi` is kept for comparison with pixy,
#  VCFtools and Stacks' `Pi`.
#
#  Per individual: heterozygous sites / sites typed -- the individual
#  "autosomal heterozygosity" of Schmidt et al. (2021).
#
#  The file is read in chunks (an all-sites VCF can be many millions of
#  lines); only per-RAD-locus sums are kept, which is also all the block
#  jackknife and bootstrap over RAD loci need (R/resampling.R). Genotypes are
#  parsed with the same helpers as read_stacks_vcf() (R/vcf_io.R).
#  Checked by tests/testthat/test-pi-allsites.R.
#
###############################################################################

## Not exported. Streaming form of read_stacks_vcf()'s "window" rule: records
## are taken in file order (VCFs are sorted by position) and a record starts
## a new locus when its CHROM differs from the previous record's or it lies
## more than `window_bp` beyond it. `state` carries the last record of the
## previous chunk, so the answer does not depend on where chunks break.
.window_stream <- function(chrom, pos, window_bp, state) {
  n <- length(chrom)
  previous_chrom <- c(state$chrom, chrom[-n])
  previous_pos <- c(state$pos, pos[-n])
  starts_locus <- is.na(previous_chrom) | chrom != previous_chrom |
    pos - previous_pos > window_bp
  start_labels <- c(state$label,
                    paste0(chrom, ":", format(pos, scientific = FALSE, trim = TRUE))[starts_locus])
  labels <- start_labels[cumsum(starts_locus) + 1L]
  list(labels = labels, state = list(chrom = chrom[n], pos = pos[n], label = labels[n]))
}

#' Nucleotide diversity and divergence from an all-sites VCF
#'
#' Per-population nucleotide diversity (pi), between-population divergence
#' (dxy and net divergence Da) and per-individual heterozygosity per
#' sequenced site, from a VCF that holds every sequenced site -- variable or
#' not -- such as Stacks' `populations --vcf-all` output (Stacks 2.62 and
#' later). Missing data are handled site by site, the way pixy does it
#' (Korunes & Samuk 2021): every site contributes its own number of pairwise
#' differences and comparisons among the gene copies typed there.
#'
#' @details
#' **Why use it.** [diversity_stats()] turns He per variant record into a
#' per-site value by dividing by the number of sequenced sites. That cannot
#' account for missing genotypes at invariant sites. With an all-sites VCF,
#' this function computes pi directly, with no such approximation. With no
#' missing data the two agree: `pi` here equals the per-site average of the
#' 2n-corrected gene diversity (`He_2n_corr`) over every sequenced site, and
#' `pi_nc` the per-site average of Nei & Chesser's `He`.
#'
#' **No allele-frequency filter.** Stacks' `--min-mac`/`--min-maf` turn a
#' failing SNP into a fixed site, so its diversity is lost while the site is
#' still counted. Under a neutral site-frequency spectrum, sites with minor
#' allele count 2 or less carry about 4/(N - 1) of pi for N gene copies (21%
#' at 10 diploids, 10% at 20), so leave those flags off for this run.
#'
#' **Estimators.** For a population at a site with allele counts `c_a` among
#' `m` typed gene copies: comparisons `C(m, 2)`, differences
#' `C(m, 2) - sum C(c_a, 2)`; pi is total differences over total comparisons.
#' Between populations: comparisons `m_x m_y`, differences
#' `m_x m_y - sum c_xa c_ya` (dxy); `Da = dxy - (pi_x + pi_y) / 2`. Per
#' individual, heterozygous sites over typed sites (Schmidt et al. 2021's
#' autosomal heterozygosity).
#'
#' **`pi` and `pi_nc`.** `pi` is the pixy estimator above. It counts the two
#' gene copies inside each individual as a comparison, so, like Stacks' `Pi`,
#' it runs low when individuals are inbred (FIS > 0): its expectation is
#' the true value times `1 - FIS / (2n - 1)` for n individuals. `pi_nc` leaves those within-individual
#' comparisons out and compares only copies from different individuals; at
#' each site it equals Nei & Chesser's (1983) gene diversity, the `He` of
#' [diversity_stats()], and it stays unbiased at any FIS. Report `pi_nc` as
#' the estimate and `pi` when comparing with pixy, VCFtools or Stacks. Standard errors are a delete-one-block
#' jackknife and intervals a block bootstrap, both over RAD loci (see
#' `locus_from` in [read_stacks_vcf()]).
#'
#' **Memory.** The file is read `chunk_lines` records at a time and only
#' per-RAD-locus sums are kept, so it is never held in memory whole.
#'
#' @inheritParams diversity_stats
#' @inheritParams read_stacks_vcf
#' @param vcf Path to an all-sites VCF (`.vcf` or `.vcf.gz`). Unlike the other
#'   analysis functions, this one needs the file itself, because it reads it
#'   in chunks.
#' @param nboot Bootstrap replicates over RAD loci. Default `1000`; `0` for
#'   none (the jackknife SE is always given).
#' @param complete_sites If `TRUE`, a population uses only the sites at which
#'   every one of its individuals is typed (Schmidt et al. 2021's
#'   recommendation). Default `FALSE`: every site counts with the individuals
#'   typed there, which is the unbiased pixy estimator.
#' @param chunk_lines Records read at a time. Default `50000`.
#' @return An object of class `raddiv_pi`, a list of:
#'   \describe{
#'     \item{pi}{One row per population: `pi`, `pi_se`, `pi_lo`, `pi_hi`, the
#'       same for `pi_nc` (see Details), and `sites` (sites with at least two
#'       typed gene copies).}
#'     \item{dxy}{One row per pair: `dxy` and `da`, each with SE and interval;
#'       `NULL` for a single population.}
#'     \item{individual}{Per individual: `het_sites`, `called_sites`,
#'       `het_per_site`.}
#'     \item{settings}{The settings of this run.}
#'   }
#'   Values are stored at full precision. With `outdir`, the three tables are
#'   written to `pi_allsites_pi.<stem>.tsv`, `pi_allsites_dxy.<stem>.tsv` and
#'   `pi_allsites_individual.<stem>.tsv`.
#' @references
#' Korunes, K.L. & Samuk, K. (2021) pixy: Unbiased estimation of nucleotide
#' diversity and divergence in the presence of missing data. *Molecular
#' Ecology Resources* 21:1359-1368.
#'
#' Nei, M. & Li, W.-H. (1979) Mathematical model for studying genetic
#' variation in terms of restriction endonucleases. *PNAS* 76:5269-5273.
#'
#' Schmidt, T.L., Jasper, M.-E., Weeks, A.R. & Hoffmann, A.A. (2021) Unbiased
#' population heterozygosity estimates from genome-wide sequence data.
#' *Methods in Ecology and Evolution* 12:1888-1898.
#' @examples
#' # A tiny all-sites VCF: two RAD loci (in CHROM, as in Stacks 2 de novo
#' # output), invariant sites included, one missing genotype.
#' vcf <- tempfile(fileext = ".vcf")
#' writeLines(c(
#'   "##fileformat=VCFv4.2",
#'   "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ta1\ta2\tb1\tb2",
#'   "1\t1\t.\tA\tG\t.\tPASS\t.\tGT\t0/0\t0/1\t1/1\t1/1",
#'   "1\t2\t.\tC\t.\t.\tPASS\t.\tGT\t0/0\t0/0\t0/0\t0/0",
#'   "2\t1\t.\tT\tC\t.\tPASS\t.\tGT\t0/1\t./.\t0/0\t0/1"), vcf)
#' popmap <- list(popA = c("a1", "a2"), popB = c("b1", "b2"))
#' res <- pi_allsites(vcf, popmap, nboot = 0)
#' res$pi
#' res$dxy
#' @export
pi_allsites <- function(vcf, popmap, locus_from = "auto", window_bp = 1000,
                        nboot = 1000L, complete_sites = FALSE, chunk_lines = 50000L,
                        outdir = NULL, stem = NULL, seed = NULL, verbose = TRUE) {
  ## ---- 1. Check the arguments ----------------------------------------------
  if (!(is.character(vcf) && length(vcf) == 1L))
    stop("`vcf` must be the path to an all-sites VCF (e.g. Stacks' populations --vcf-all ",
         "output); pi_allsites() reads the file in chunks.", call. = FALSE)
  .check_run_inputs(vcf, popmap, stem, outdir)
  .check_choice(locus_from, "locus_from", c("auto", "ID", "CHROM", "window"))
  .check_number(window_bp, "window_bp", min = 0)
  nboot <- .check_count(nboot, "nboot")
  chunk_lines <- .check_count(chunk_lines, "chunk_lines", min = 1)
  .check_flag(complete_sites, "complete_sites")
  .check_flag(verbose, "verbose")
  .check_seed(seed)
  if (!is.null(seed)) {
    restore_rng <- .save_rng_state()
    on.exit(restore_rng(), add = TRUE)
    set.seed(seed)
  }

  ## ---- 2. Header ------------------------------------------------------------
  ## file() reads .gz files transparently.
  con <- file(vcf, "r")
  on.exit(close(con), add = TRUE)
  .inform(verbose, "Reading ", vcf, " ...")
  ## Skip the ## lines; the #CHROM line names the samples. Records read along
  ## with the header are processed as the first chunk.
  header <- NULL
  chunk <- character(0)
  repeat {
    lines <- readLines(con, n = 1000L)
    if (!length(lines)) break
    header_line <- grep("^#CHROM", lines)
    if (length(header_line)) {
      header <- strsplit(sub("^#", "", lines[header_line[1]]), "\t", fixed = TRUE)[[1]]
      chunk <- lines[-seq_len(header_line[1])]
      break
    }
  }
  if (is.null(header)) stop("No '#CHROM' header line found. Is this a VCF?", call. = FALSE)
  if (length(header) < 10L)
    stop("The #CHROM header line has no sample columns.", call. = FALSE)
  samples <- header[-(1:9)]
  pops <- .resolve_pops(popmap, samples, verbose = verbose)
  pop_names <- names(pops)
  n_pops <- length(pops)
  sample_columns <- lapply(pops, function(ids) match(ids, samples))
  pairs <- if (n_pops >= 2) utils::combn(n_pops, 2, simplify = FALSE) else list()
  n_pairs <- length(pairs)

  ## Per-record pieces, in columns:
  ##   1 .. 2r                    differences, comparisons for each population
  ##   2r + 1 .. 2r + 2P          differences, comparisons for each pair
  ##   2r + 2P + 1 .. 3r + 2P     1 if the site has a comparison in the population
  ##   3r + 2P + 1 .. 5r + 2P     the same differences and comparisons among gene
  ##                              copies of DIFFERENT individuals only (pi_nc)
  ## (r = populations, P = pairs)
  n_pieces <- 5 * n_pops + 2 * n_pairs
  pop_diff_col <- 2 * seq_len(n_pops) - 1
  pop_comp_col <- 2 * seq_len(n_pops)
  pair_diff_col <- 2 * n_pops + 2 * seq_len(n_pairs) - 1
  pair_comp_col <- 2 * n_pops + 2 * seq_len(n_pairs)
  site_used_col <- 2 * n_pops + 2 * n_pairs + seq_len(n_pops)
  nc_diff_col <- 3 * n_pops + 2 * n_pairs + 2 * seq_len(n_pops) - 1
  nc_comp_col <- 3 * n_pops + 2 * n_pairs + 2 * seq_len(n_pops)

  ## ---- 3. Read chunks, keeping only per-locus sums --------------------------
  locus_sums <- NULL
  het_sites <- called_sites <- stats::setNames(numeric(length(samples)), samples)
  n_rec <- 0
  rule <- NULL
  window_state <- list(chrom = NA_character_, pos = -Inf, label = NA_character_)
  repeat {
    chunk <- chunk[nzchar(chunk)]
    if (length(chunk)) {
      fields <- strsplit(chunk, "\t", fixed = TRUE)
      n_fields <- lengths(fields)
      if (any(n_fields != length(header))) {
        bad <- which(n_fields != length(header))[1]
        stop(sprintf("Malformed VCF: record %d has %d tab-separated fields, expected %d.",
                     n_rec + bad, n_fields[bad], length(header)), call. = FALSE)
      }
      fields <- do.call(rbind, fields)
      .check_gt_first(fields[, 9], first_record = n_rec + 1)

      ## RAD locus of each record (see locus_from in ?read_stacks_vcf). The
      ## rule is decided on the first chunk and then applied to every chunk.
      if (is.null(rule)) rule <- .choose_locus_rule(locus_from, fields[, 3])
      labels <- switch(rule,
        ID = {
          .choose_locus_rule("ID", fields[, 3])       # stops if an ID is missing
          sub(":.*$", "", fields[, 3])
        },
        CHROM = fields[, 1],
        window = {
          pos <- suppressWarnings(as.numeric(fields[, 2]))
          if (anyNA(pos)) stop("Some POS values are not numbers.", call. = FALSE)
          streamed <- .window_stream(fields[, 1], pos, window_bp, window_state)
          window_state <- streamed$state
          streamed$labels
        })

      gt <- .parse_gt(fields[, -(1:9), drop = FALSE])
      k <- max(c(1L, gt$A1, gt$A2), na.rm = TRUE)
      n_chunk <- nrow(fields)

      pieces <- matrix(0, n_chunk, n_pieces)
      counts <- vector("list", n_pops)
      for (p in seq_len(n_pops)) {
        a1 <- gt$A1[, sample_columns[[p]], drop = FALSE]
        a2 <- gt$A2[, sample_columns[[p]], drop = FALSE]
        typed <- !is.na(a1)
        use <- if (complete_sites) rowSums(typed) == ncol(a1) else rep(TRUE, n_chunk)
        C <- .allele_counts(a1, a2, k)
        C[!use, ] <- 0L
        counts[[p]] <- C
        copies <- rowSums(C)
        comparisons <- copies * (copies - 1) / 2
        differences <- comparisons - rowSums(C * (C - 1) / 2)
        pieces[, pop_diff_col[p]] <- differences
        pieces[, pop_comp_col[p]] <- comparisons
        pieces[, site_used_col[p]] <- comparisons > 0
        typed_used <- typed & use
        ## pi_nc: leave out the one pair of gene copies inside each typed
        ## individual (a difference exactly when it is heterozygous).
        pieces[, nc_diff_col[p]] <- differences - rowSums(typed_used & (a1 != a2), na.rm = TRUE)
        pieces[, nc_comp_col[p]] <- comparisons - rowSums(typed_used)
        het_sites[sample_columns[[p]]] <- het_sites[sample_columns[[p]]] +
          colSums(typed_used & (a1 != a2))
        called_sites[sample_columns[[p]]] <- called_sites[sample_columns[[p]]] +
          colSums(typed_used)
      }
      for (q in seq_len(n_pairs)) {
        C_x <- counts[[pairs[[q]][1]]]
        C_y <- counts[[pairs[[q]][2]]]
        comparisons <- rowSums(C_x) * rowSums(C_y)
        pieces[, pair_diff_col[q]] <- comparisons - rowSums(C_x * C_y)
        pieces[, pair_comp_col[q]] <- comparisons
      }
      chunk_sums <- rowsum(pieces, labels, reorder = FALSE)
      locus_sums <- if (is.null(locus_sums)) chunk_sums
                    else rowsum(rbind(locus_sums, chunk_sums),
                                c(rownames(locus_sums), rownames(chunk_sums)), reorder = FALSE)
      n_rec <- n_rec + n_chunk
    }
    chunk <- readLines(con, n = chunk_lines)
    if (!length(chunk)) break
  }
  if (is.null(locus_sums)) stop("No records found after the VCF header.", call. = FALSE)

  S <- unname(locus_sums)
  n_loci <- nrow(S)
  .inform(verbose, sprintf("  %s sites on %s RAD loci (loci from %s)", .big(n_rec), .big(n_loci),
                           .describe_locus_rule(rule, window_bp)))

  ## ---- 4. Statistics from sums, jackknife and bootstrap ---------------------
  pair_labels <- vapply(pairs, function(q) paste(pop_names[q], collapse = "__"), character(1))
  stat_names <- c(paste0("pi_", pop_names), paste0("pinc_", pop_names),
                  if (n_pairs) c(paste0("dxy_", pair_labels), paste0("da_", pair_labels)))
  stats_from_sums <- function(totals) {
    totals <- matrix(totals, ncol = n_pieces)
    pi_by_pop <- totals[, pop_diff_col, drop = FALSE] / totals[, pop_comp_col, drop = FALSE]
    pi_nc <- totals[, nc_diff_col, drop = FALSE] / totals[, nc_comp_col, drop = FALSE]
    out <- cbind(pi_by_pop, pi_nc)
    if (n_pairs) {
      dxy <- totals[, pair_diff_col, drop = FALSE] / totals[, pair_comp_col, drop = FALSE]
      mean_pi <- matrix(vapply(pairs, function(q) (pi_by_pop[, q[1]] + pi_by_pop[, q[2]]) / 2,
                               numeric(nrow(totals))), nrow = nrow(totals))
      out <- cbind(out, dxy, dxy - mean_pi)
    }
    out[!is.finite(out)] <- NA_real_
    colnames(out) <- stat_names
    out
  }
  point <- stats_from_sums(colSums(S))[1L, ]
  se <- if (n_loci > 1) .jack_block_sums(S, stats_from_sums) else point * NA
  boot_matrix <- NULL
  if (nboot > 0 && n_loci > 1) {
    .inform(verbose, sprintf("Bootstrapping %s replicates over %s RAD loci ...",
                             .big(nboot), .big(n_loci)))
    boot_matrix <- .boot_block_sums(S, nboot, stats_from_sums)
  }
  ci <- .percentile_ci(boot_matrix, stat_names)

  ## ---- 5. Result tables ------------------------------------------------------
  pi_names <- paste0("pi_", pop_names)
  nc_names <- paste0("pinc_", pop_names)
  pi_table <- data.frame(population = pop_names, pi = unname(point[pi_names]),
                         pi_se = unname(se[pi_names]), pi_lo = unname(ci[pi_names, "lo"]),
                         pi_hi = unname(ci[pi_names, "hi"]),
                         pi_nc = unname(point[nc_names]), pi_nc_se = unname(se[nc_names]),
                         pi_nc_lo = unname(ci[nc_names, "lo"]), pi_nc_hi = unname(ci[nc_names, "hi"]),
                         sites = as.integer(colSums(S[, site_used_col, drop = FALSE])),
                         row.names = NULL)
  dxy_table <- NULL
  if (n_pairs) {
    dxy_names <- paste0("dxy_", pair_labels)
    da_names <- paste0("da_", pair_labels)
    dxy_table <- data.frame(
      pop1 = vapply(pairs, function(q) pop_names[q[1]], character(1)),
      pop2 = vapply(pairs, function(q) pop_names[q[2]], character(1)),
      dxy = unname(point[dxy_names]), dxy_se = unname(se[dxy_names]),
      dxy_lo = unname(ci[dxy_names, "lo"]), dxy_hi = unname(ci[dxy_names, "hi"]),
      da = unname(point[da_names]), da_se = unname(se[da_names]),
      da_lo = unname(ci[da_names, "lo"]), da_hi = unname(ci[da_names, "hi"]),
      row.names = NULL)
  }
  ids <- unlist(pops, use.names = FALSE)
  individual <- data.frame(
    sample = ids, population = rep(pop_names, lengths(pops)),
    het_sites = as.integer(het_sites[ids]), called_sites = as.integer(called_sites[ids]),
    het_per_site = unname(ifelse(called_sites[ids] > 0, het_sites[ids] / called_sites[ids], NA_real_)),
    row.names = NULL)

  written <- character(0)
  if (!is.null(outdir)) {
    stem <- .derive_stem(vcf, stem)
    tables <- list(pi = pi_table, dxy = dxy_table, individual = individual)
    written <- .write_tables(lapply(tables, .round_pi_table), outdir,
                             stats::setNames(sprintf("pi_allsites_%s.%s.tsv", names(tables), stem),
                                             names(tables)))
  }
  settings <- list(n_sites = n_rec, n_loci = n_loci, nboot = nboot,
                   complete_sites = complete_sites, seed = seed, files = written)
  structure(list(pi = pi_table, dxy = dxy_table, individual = individual, settings = settings),
            class = "raddiv_pi")
}

## Not exported. pi_allsites() tables are shown and written to 4 significant
## digits (counts of sites are left as they are).
.round_pi_table <- function(df) {
  if (is.null(df)) return(NULL)
  .round_table(df, signif_cols = stats::setNames(rep(4L, ncol(df)), names(df)))
}

#' @rdname pi_allsites
#' @param x,object A `raddiv_pi` object, as returned by `pi_allsites()`.
#' @param ... Ignored.
#' @export
print.raddiv_pi <- function(x, ...) {
  st <- x$settings
  cat(sprintf("Nucleotide diversity from an all-sites VCF: %s sites on %s RAD loci\n",
              .big(st$n_sites), .big(st$n_loci)))
  cat("  pi and dxy per sequenced site, missing data handled site by site (pixy); _se: jackknife over RAD loci\n")
  cat("\n$pi\n")
  .print_table(.round_pi_table(x$pi))
  if (!is.null(x$dxy)) {
    cat("\n$dxy\n")
    .print_table(.round_pi_table(utils::head(x$dxy, 10)))
    if (nrow(x$dxy) > 10) cat(sprintf("  (first 10 of %d pairs; the full table is x$dxy)\n", nrow(x$dxy)))
  }
  cat("\nPer-individual heterozygosity per site is in $individual.\n")
  cat("summary() prints the full report.\n")
  invisible(x)
}

#' @rdname pi_allsites
#' @export
summary.raddiv_pi <- function(object, ...) {
  structure(list(result = object), class = "summary.raddiv_pi")
}

#' @export
print.summary.raddiv_pi <- function(x, ...) {
  res <- x$result
  st <- res$settings
  note <- function(...) cat(paste0("  ", c(...), "\n"), sep = "")
  rule <- strrep("=", 69)
  cat("\n", rule, "\n", sep = "")
  cat("  NUCLEOTIDE DIVERSITY FROM AN ALL-SITES VCF --", .big(st$n_sites),
      "sites on", .big(st$n_loci), "RAD loci\n")
  note(.report_text$pi_definition)
  if (st$complete_sites)
    note("complete_sites: each population uses only sites typed in all its individuals")
  if (st$nboot > 0)
    cat("  95% CI from ", .big(st$nboot), " bootstrap replicates over RAD loci\n", sep = "")
  cat(rule, "\n\n", sep = "")
  .print_table(.round_pi_table(res$pi))
  note(.report_text$pi_columns)
  if (!is.null(res$dxy)) {
    cat("\nBetween populations: dxy (absolute divergence) and Da = dxy - mean(pi) (net)\n")
    .print_table(.round_pi_table(res$dxy))
  }
  cat("\nPer-individual heterozygosity per sequenced site (Schmidt et al. 2021),\n")
  cat("population means:\n")
  mean_het <- tapply(res$individual$het_per_site, res$individual$population, mean, na.rm = TRUE)
  .print_table(data.frame(population = names(mean_het),
                          mean_het_per_site = signif(as.numeric(mean_het), 4)))
  note("(the full table is x$individual)")
  if (length(st$files)) cat("\nWrote ", paste(st$files, collapse = ", "), "\n", sep = "")
  cat("\n")
  invisible(x)
}
