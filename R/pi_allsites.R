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
#  Per individual: heterozygous sites / sites typed -- the individual
#  "autosomal heterozygosity" of Schmidt et al. (2021).
#
#  The file is read in chunks (an all-sites VCF can be many millions of
#  lines); only per-RAD-locus sums are kept, which is also all the block
#  jackknife and bootstrap over RAD loci need (R/resampling.R).
#
###############################################################################

## Not exported. Streaming form of read_stacks_vcf()'s "window" rule: records
## are taken in file order (VCFs are sorted by position) and a record starts
## a new locus when its CHROM differs from the previous record's or it lies
## more than `window_bp` beyond it. `state` carries the last record of the
## previous chunk, so the answer does not depend on where chunks break.
.window_stream <- function(chrom, pos, window_bp, state) {
  n <- length(chrom)
  prev_c <- c(state$chrom, chrom[-n]); prev_p <- c(state$pos, pos[-n])
  new <- is.na(prev_c) | chrom != prev_c | pos - prev_p > window_bp
  starts <- c(state$label,
              paste0(chrom, ":", format(pos, scientific = FALSE, trim = TRUE))[new])
  lab <- starts[cumsum(new) + 1L]
  list(labels = lab, state = list(chrom = chrom[n], pos = pos[n], label = lab[n]))
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
#' missing data the two agree: pi here equals the per-site average of the
#' 2n-corrected gene diversity (`He_2n_corr`) over every sequenced site.
#'
#' **Estimators.** For a population at a site with allele counts `c_a` among
#' `m` typed gene copies: comparisons `C(m, 2)`, differences
#' `C(m, 2) - sum C(c_a, 2)`; pi is total differences over total comparisons.
#' Between populations: comparisons `m_x m_y`, differences
#' `m_x m_y - sum c_xa c_ya` (dxy); `Da = dxy - (pi_x + pi_y) / 2`. Per
#' individual, heterozygous sites over typed sites (Schmidt et al. 2021's
#' autosomal heterozygosity). Standard errors are a delete-one-block
#' jackknife and intervals a block bootstrap, both over RAD loci (see
#' `locus_from` in [read_stacks_vcf()]).
#'
#' **Memory.** The file is read `chunk_lines` records at a time and only
#' per-RAD-locus sums are kept, so it is never held in memory whole.
#'
#' @inheritParams diversity_stats
#' @inheritParams read_stacks_vcf
#' @param vcf_file Path to an all-sites VCF (`.vcf` or `.vcf.gz`).
#' @param nboot Bootstrap replicates over RAD loci. Default `1000L`; `0` for
#'   none (the jackknife SE is always given).
#' @param complete_sites If `TRUE`, a population uses only the sites at which
#'   every one of its individuals is typed (Schmidt et al. 2021's
#'   recommendation). Default `FALSE`: every site counts with the individuals
#'   typed there, which is the unbiased pixy estimator.
#' @param chunk_lines Records read at a time. Default `50000L`.
#' @return An object of class `raddiv_pi` (printing it shows a report): a
#'   list with `pi` (one row per population: `pi`, `pi_se`, `pi_lo`,
#'   `pi_hi`, `sites` used), `dxy` (one row per pair: `dxy` and `da`, each
#'   with SE and CI; `NULL` for a single population) and `individual`
#'   (`het_sites`, `called_sites`, `het_per_site` per individual). With
#'   `outdir`, the three tables are written to `pi_allsites_pi.<stem>.tsv`,
#'   `pi_allsites_dxy.<stem>.tsv` and `pi_allsites_individual.<stem>.tsv`.
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
#' popmap <- tempfile()
#' writeLines(c("a1\tpopA", "a2\tpopA", "b1\tpopB", "b2\tpopB"), popmap)
#' res <- pi_allsites(vcf, popmap, nboot = 0)
#' res$pi
#' res$dxy
#' @export
pi_allsites <- function(vcf_file, popmap_f, locus_from = "auto", window_bp = 1000,
                        nboot = 1000L, complete_sites = FALSE, chunk_lines = 50000L,
                        outdir = NULL, stem = NULL, seed = 2024, verbose = TRUE) {
  if (!(is.character(vcf_file) && length(vcf_file) == 1L))
    stop("vcf_file must be the path to an all-sites VCF (e.g. Stacks' populations --vcf-all output).")
  .check_run_inputs(vcf_file, popmap_f, stem)
  if (!(length(locus_from) == 1L && locus_from %in% c("auto", "ID", "CHROM", "window")))
    stop("locus_from must be one of \"auto\", \"ID\", \"CHROM\", \"window\".")
  nboot <- as.integer(nboot); chunk_lines <- as.integer(chunk_lines)
  if (is.na(nboot) || nboot < 0) stop("nboot must be a non-negative integer.")
  if (is.na(chunk_lines) || chunk_lines < 1) stop("chunk_lines must be a positive integer.")
  restore_rng <- .save_rng_state()
  on.exit(restore_rng(), add = TRUE)
  set.seed(seed)

  con <- if (grepl("\\.gz$", vcf_file)) gzfile(vcf_file, "r") else file(vcf_file, "r")
  on.exit(close(con), add = TRUE)
  if (verbose) message("Reading ", vcf_file, " ...")

  ## Header: skip the ## lines; the #CHROM line names the samples. Any records
  ## read along with it are processed as the first chunk.
  hdr <- NULL; chunk <- character(0)
  repeat {
    x <- readLines(con, n = 1000L)
    if (!length(x)) break
    h <- grep("^#CHROM", x)
    if (length(h)) {
      hdr <- strsplit(sub("^#", "", x[h[1]]), "\t", fixed = TRUE)[[1]]
      chunk <- x[-seq_len(h[1])]
      break
    }
  }
  if (is.null(hdr)) stop("No '#CHROM' header line found. Is this a VCF?")
  samples <- hdr[-(1:9)]
  pops <- read_popmap(popmap_f, samples, verbose = verbose)
  r <- length(pops)
  cols <- lapply(pops, function(ids) match(ids, samples))
  pairs <- if (r >= 2) utils::combn(r, 2, simplify = FALSE) else list()
  npair <- length(pairs)
  ## Per-record pieces, in columns: for each population (differences,
  ## comparisons), then for each pair (differences, comparisons), then for
  ## each population 1 if the site has any comparison there (sites used).
  n_piece <- 2 * r + 2 * npair + r

  acc <- NULL
  het <- called <- stats::setNames(numeric(length(samples)), samples)
  n_rec <- 0; rule <- NULL
  state <- list(chrom = NA_character_, pos = -Inf, label = NA_character_)
  repeat {
    chunk <- chunk[nzchar(chunk)]
    if (length(chunk)) {
      f <- strsplit(chunk, "\t", fixed = TRUE)
      nf <- lengths(f)
      if (any(nf != length(hdr))) {
        bad <- which(nf != length(hdr))[1]
        stop(sprintf("Malformed VCF: record %d has %d tab-separated fields, expected %d.",
                     n_rec + bad, nf[bad], length(hdr)))
      }
      f <- do.call(rbind, f)
      if (any(!grepl("^GT(:|$)", f[, 9])))
        stop("Malformed VCF: every record's FORMAT must start with GT.")

      ## RAD locus of each record (see locus_from in ?read_stacks_vcf).
      if (is.null(rule))
        rule <- if (locus_from != "auto") locus_from
                else if (all(nzchar(f[, 3]) & f[, 3] != ".")) "ID" else "window"
      lab <- switch(rule,
        ID = {
          if (any(!nzchar(f[, 3]) | f[, 3] == "."))
            stop("Using the ID column for RAD loci, but some records have no ID ",
                 "(\".\"). Use locus_from = \"CHROM\" or \"window\".")
          sub(":.*$", "", f[, 3])
        },
        CHROM = f[, 1],
        window = {
          pos <- suppressWarnings(as.numeric(f[, 2]))
          if (anyNA(pos)) stop("Some POS values are not numbers.")
          st <- .window_stream(f[, 1], pos, window_bp, state)
          state <- st$state
          st$labels
        })

      ## Genotypes: anything that is not a clean diploid call is missing.
      g <- sub(":.*$", "", f[, -(1:9), drop = FALSE])
      ok <- grepl("^[0-9]+[/|][0-9]+$", g)
      a1 <- a2 <- matrix(NA_integer_, nrow(g), ncol(g))
      a1[ok] <- as.integer(sub("[/|].*$", "", g[ok])) + 1L
      a2[ok] <- as.integer(sub("^.*[/|]", "", g[ok])) + 1L
      k <- max(c(1L, a1, a2), na.rm = TRUE)
      nr <- nrow(g)

      pieces <- matrix(0, nr, n_piece)
      cnts <- vector("list", r)
      for (p in seq_len(r)) {
        A1 <- a1[, cols[[p]], drop = FALSE]; A2 <- a2[, cols[[p]], drop = FALSE]
        typed <- !is.na(A1)
        use <- if (complete_sites) rowSums(typed) == ncol(A1) else rep(TRUE, nr)
        cnt <- .allele_counts(A1, A2, k)
        cnt[!use, ] <- 0
        cnts[[p]] <- cnt
        m <- rowSums(cnt)
        comps <- m * (m - 1) / 2
        pieces[, 2 * p - 1] <- comps - rowSums(cnt * (cnt - 1) / 2)
        pieces[, 2 * p]     <- comps
        pieces[, 2 * r + 2 * npair + p] <- comps > 0
        tu <- typed & use
        het[cols[[p]]]    <- het[cols[[p]]] + colSums(tu & (A1 != A2))
        called[cols[[p]]] <- called[cols[[p]]] + colSums(tu)
      }
      for (q in seq_len(npair)) {
        cx <- cnts[[pairs[[q]][1]]]; cy <- cnts[[pairs[[q]][2]]]
        comps <- rowSums(cx) * rowSums(cy)
        pieces[, 2 * r + 2 * q - 1] <- comps - rowSums(cx * cy)
        pieces[, 2 * r + 2 * q]     <- comps
      }
      part <- rowsum(pieces, lab, reorder = FALSE)
      acc <- if (is.null(acc)) part
             else rowsum(rbind(acc, part), c(rownames(acc), rownames(part)), reorder = FALSE)
      n_rec <- n_rec + nr
    }
    chunk <- readLines(con, n = chunk_lines)
    if (!length(chunk)) break
  }
  if (is.null(acc)) stop("No records found after the VCF header.")

  S <- unname(acc); nL <- nrow(S)
  if (verbose)
    message(sprintf("  %s sites on %s RAD loci (%s)", format(n_rec, big.mark = ","),
                    format(nL, big.mark = ","),
                    switch(rule, ID = "loci from the ID column", CHROM = "loci from CHROM",
                           window = sprintf("loci from CHROM + POS within %s bp",
                                            format(window_bp, big.mark = ",")))))
  pn <- names(pops)
  pair_lab <- vapply(pairs, function(q) paste(pn[q], collapse = "__"), character(1))
  stat_names <- c(paste0("pi_", pn),
                  if (npair) c(paste0("dxy_", pair_lab), paste0("da_", pair_lab)))
  stat_mat <- function(tot) {
    tot <- matrix(tot, ncol = n_piece)
    pi <- tot[, 2 * seq_len(r) - 1, drop = FALSE] / tot[, 2 * seq_len(r), drop = FALSE]
    out <- pi
    if (npair) {
      dxy <- tot[, 2 * r + 2 * seq_len(npair) - 1, drop = FALSE] /
             tot[, 2 * r + 2 * seq_len(npair), drop = FALSE]
      mean_pi <- matrix(vapply(pairs, function(q) (pi[, q[1]] + pi[, q[2]]) / 2,
                               numeric(nrow(tot))), nrow = nrow(tot))
      out <- cbind(out, dxy, dxy - mean_pi)
    }
    out[!is.finite(out)] <- NA_real_
    colnames(out) <- stat_names
    out
  }
  point <- stat_mat(colSums(S))[1, ]
  se <- if (nL > 1) .jack_block_sums(S, stat_mat) else point * NA
  ci <- matrix(NA_real_, length(point), 2, dimnames = list(names(point), NULL))
  if (nboot > 0 && nL > 1) {
    if (verbose) message(sprintf("Bootstrapping %s replicates over %s RAD loci ...",
                                 format(nboot, big.mark = ","), format(nL, big.mark = ",")))
    bt <- .boot_block_sums(S, nboot, stat_mat)
    ci <- t(apply(bt, 2, stats::quantile, c(0.025, 0.975), na.rm = TRUE))
  }

  sites_used <- colSums(S[, 2 * r + 2 * npair + seq_len(r), drop = FALSE])
  pi_tab <- data.frame(population = pn,
                       pi = signif(point[paste0("pi_", pn)], 4),
                       pi_se = signif(se[paste0("pi_", pn)], 4),
                       pi_lo = signif(ci[paste0("pi_", pn), 1], 4),
                       pi_hi = signif(ci[paste0("pi_", pn), 2], 4),
                       sites = as.integer(sites_used), row.names = NULL)
  dxy_tab <- NULL
  if (npair) {
    d <- paste0("dxy_", pair_lab); a <- paste0("da_", pair_lab)
    dxy_tab <- data.frame(
      pop1 = vapply(pairs, function(q) pn[q[1]], ""), pop2 = vapply(pairs, function(q) pn[q[2]], ""),
      dxy = signif(point[d], 4), dxy_se = signif(se[d], 4),
      dxy_lo = signif(ci[d, 1], 4), dxy_hi = signif(ci[d, 2], 4),
      da = signif(point[a], 4), da_se = signif(se[a], 4),
      da_lo = signif(ci[a, 1], 4), da_hi = signif(ci[a, 2], 4), row.names = NULL)
  }
  ids <- unlist(pops, use.names = FALSE)
  ind <- data.frame(sample = ids, population = rep(pn, lengths(pops)),
                    het_sites = het[ids], called_sites = called[ids],
                    het_per_site = signif(ifelse(called[ids] > 0, het[ids] / called[ids], NA_real_), 4),
                    row.names = NULL)

  written <- character(0)
  if (!is.null(outdir)) {
    dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
    stem <- .derive_stem(vcf_file, stem)
    out_tabs <- list(pi = pi_tab, dxy = dxy_tab, individual = ind)
    for (nm in names(out_tabs)) {
      if (is.null(out_tabs[[nm]])) next
      fl <- file.path(outdir, sprintf("pi_allsites_%s.%s.tsv", nm, stem))
      utils::write.table(out_tabs[[nm]], fl, sep = "\t", quote = FALSE, row.names = FALSE)
      written <- c(written, fl)
    }
  }
  structure(list(pi = pi_tab, dxy = dxy_tab, individual = ind),
            class = "raddiv_pi",
            report = list(n_sites = n_rec, n_loci = nL, nboot = nboot,
                          complete_sites = complete_sites, files = written))
}

#' @rdname pi_allsites
#' @param x A `raddiv_pi` object, as returned by `pi_allsites()`.
#' @param ... Ignored.
#' @export
print.raddiv_pi <- function(x, ...) {
  rp <- attr(x, "report")
  old <- options(width = max(200L, getOption("width")))
  on.exit(options(old), add = TRUE)
  cat("\n=====================================================================\n")
  cat("  NUCLEOTIDE DIVERSITY FROM AN ALL-SITES VCF --", format(rp$n_sites, big.mark = ","),
      "sites on", format(rp$n_loci, big.mark = ","), "RAD loci\n")
  cat("  pi = sum of pairwise differences / sum of pairwise comparisons, site by\n")
  cat("  site, among the gene copies typed there (pixy; Korunes & Samuk 2021)\n")
  if (rp$complete_sites) cat("  complete_sites: each population uses only sites typed in all its individuals\n")
  if (rp$nboot > 0) cat("  95% CI from ", format(rp$nboot, big.mark = ","),
                        " bootstrap replicates over RAD loci\n", sep = "")
  cat("=====================================================================\n\n")
  print(x$pi, row.names = FALSE)
  cat("  _se: delete-one-block jackknife over RAD loci; sites: sites with at least\n")
  cat("  two typed gene copies in that population.\n")
  if (!is.null(x$dxy)) {
    cat("\nBetween populations: dxy (absolute divergence) and Da = dxy - mean(pi) (net)\n")
    print(x$dxy, row.names = FALSE)
  }
  cat("\nPer-individual heterozygosity per sequenced site (Schmidt et al. 2021),\n")
  cat("population means:\n")
  mh <- tapply(x$individual$het_per_site, x$individual$population, mean, na.rm = TRUE)
  print(data.frame(population = names(mh), mean_het_per_site = signif(as.numeric(mh), 4)),
        row.names = FALSE)
  cat("  (the full table is x$individual)\n")
  if (length(rp$files)) cat("\nWrote ", paste(rp$files, collapse = ", "), "\n", sep = "")
  cat("\n")
  invisible(x)
}
