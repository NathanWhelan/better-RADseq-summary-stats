###############################################################################
#
#  R/diversity_table.R -- one table for a manuscript from the two
#  diversity_stats() runs (SNP VCF and haplotype VCF), with a caption.
#
#  The cells are exactly the ones the short summary() of each run prints
#  (summary(x)$tables), so the table holds the recommended standard error for
#  each statistic and nothing else. The two runs are joined by population
#  name, never by row order, and the rows follow the popmap.
#
###############################################################################

#' A diversity table for a manuscript
#'
#' Joins the results of the two [diversity_stats()] runs into one table, one
#' row per population, and writes a caption to go with it. The table holds
#' exactly what the short `summary()` of each run shows: each statistic as
#' `estimate (SE)`, with the one standard error the package recommends.
#'
#' | column | from | SE |
#' |---|---|---|
#' | `Ho (SE)`, `He (SE)` | SNP VCF | RAD loci and individuals |
#' | `Ho_autosomal (SE)`, `He_autosomal (SE)` (per sequenced site; only with `sites`) | SNP VCF | RAD loci and individuals |
#' | `Fis (SE)` | haplotype VCF | RAD loci and individuals |
#' | `Ar (SE)`, `rarefied private alleles (SE)` | haplotype VCF | RAD loci only |
#'
#' With `se_individuals = FALSE` in [diversity_stats()], the SEs count RAD
#' loci only, and the caption says so.
#'
#' Other columns (the other SEs and intervals, `pct_poly`, `privAr`, and the
#' counts behind each value) stay in the results and the `outdir` files, for a
#' supplement. Report them only if readers need them to judge the table.
#'
#' Give both runs, or only one: with only a SNP VCF (e.g. from ipyrad or
#' dDocent), pass `snps` alone.
#'
#' @param snps The result of [diversity_stats()] on the SNP VCF
#'   (`populations.snps.vcf`). Optional.
#' @param haps The result of [diversity_stats()] on the haplotype VCF
#'   (`populations.haps.vcf`). Optional.
#' @return A data frame (class `raddiv_table`) with one row per population, in
#'   popmap order. Printing it shows the table and then the caption; the
#'   caption is also in `attr(x, "caption")`. `write.csv()` and
#'   `knitr::kable()` work on it as on any data frame.
#' @seealso [diversity_stats()], "Tables for a paper".
#' @examples
#' # Toy data: 2 populations of 4 and 3 individuals (far too few for a study).
#' snps   <- system.file("extdata", "small.snps.vcf", package = "RADdiversity")
#' haps   <- system.file("extdata", "small.haps.vcf", package = "RADdiversity")
#' popmap <- system.file("extdata", "small_popmap.tsv", package = "RADdiversity")
#' res_snps <- diversity_stats(snps, popmap, g = 4, nboot = 0, se_individuals = FALSE,
#'                             verbose = FALSE)
#' res_haps <- diversity_stats(haps, popmap, g = 4, nboot = 0, se_individuals = FALSE,
#'                             verbose = FALSE)
#' table1 <- diversity_table(res_snps, res_haps)
#' table1                        # the table, then its caption
#' attr(table1, "caption")       # the caption on its own
#' out <- tempfile(fileext = ".csv")
#' write.csv(table1, out, row.names = FALSE)
#' @export
diversity_table <- function(snps = NULL, haps = NULL) {
  if (is.null(snps) && is.null(haps))
    stop("Give at least one diversity_stats() result: `snps` (from the SNP VCF), ",
         "`haps` (from the haplotype VCF), or both.", call. = FALSE)
  .check_diversity_run(snps, "snps", haplotype = FALSE)
  .check_diversity_run(haps, "haps", haplotype = TRUE)

  first <- if (!is.null(snps)) snps else haps
  pops <- first$per_population$population
  if (!is.null(snps) && !is.null(haps)) .check_same_populations(snps, haps)

  ## The cells of each run's short summary, joined by population name.
  pieces <- list()
  if (!is.null(snps)) {
    tables <- summary(snps)$tables
    pieces <- c(pieces, list(tables$results), if (!is.null(tables$per_site)) list(tables$per_site))
  }
  if (!is.null(haps)) {
    results <- summary(haps)$tables$results
    if (!is.null(snps)) results$n <- NULL          # `n` is in the SNP table already
    pieces <- c(pieces, list(results))
  }
  out <- data.frame(population = pops, stringsAsFactors = FALSE)
  for (piece in pieces)
    out <- cbind(out, piece[match(pops, piece$population), setdiff(names(piece), "population"),
                            drop = FALSE])
  rownames(out) <- NULL
  structure(out, caption = .diversity_caption(snps, haps),
            class = c("raddiv_table", "data.frame"))
}

#' @rdname diversity_table
#' @param x A `raddiv_table`, as returned by `diversity_table()`.
#' @param ... Ignored.
#' @export
print.raddiv_table <- function(x, ...) {
  df <- x
  class(df) <- "data.frame"
  attr(df, "caption") <- NULL
  .print_compact(df)
  caption <- attr(x, "caption")
  if (!is.null(caption)) {
    cat("\nCaption:\n")
    cat(strwrap(caption, width = 76, indent = 2, exdent = 2), sep = "\n")
  }
  invisible(x)
}

## Not exported. Stops unless `x` is NULL or a diversity_stats() result from
## the right kind of VCF.
.check_diversity_run <- function(x, name, haplotype) {
  if (is.null(x)) return(invisible(NULL))
  if (!inherits(x, "raddiv_diversity"))
    stop("`", name, "` must be the result of diversity_stats() (got an object of class ",
         paste(class(x), collapse = "/"), ").", call. = FALSE)
  if (!identical(isTRUE(x$settings$is_haplotype), haplotype)) {
    other <- if (haplotype) "snps" else "haps"
    stop("`", name, "` is a ", if (haplotype) "SNP" else "haplotype", "-VCF result, but `", name,
         "` should come from the ", if (haplotype) "haplotype" else "SNP", " VCF. Did you swap `",
         name, "` and `", other, "`?", call. = FALSE)
  }
  invisible(NULL)
}

## Not exported. The two runs must describe the same populations: stops and
## names any population found in only one run. Warns when a population has a
## different number of individuals in the two runs (a different popmap).
.check_same_populations <- function(snps, haps) {
  ps <- snps$per_population
  ph <- haps$per_population
  only_snps <- setdiff(ps$population, ph$population)
  only_haps <- setdiff(ph$population, ps$population)
  if (length(only_snps) || length(only_haps))
    stop("The two runs have different populations.",
         if (length(only_snps)) paste0("\n  Only in the SNP run: ", paste(only_snps, collapse = ", ")),
         if (length(only_haps)) paste0("\n  Only in the haplotype run: ", paste(only_haps, collapse = ", ")),
         "\n  Run diversity_stats() on both VCFs with the same popmap.", call. = FALSE)
  n_haps <- ph$n[match(ps$population, ph$population)]
  differ <- ps$n != n_haps
  if (any(differ))
    warning("The two runs used different individuals: ",
            paste(sprintf("%s has %d in the SNP run and %d in the haplotype run",
                          ps$population[differ], ps$n[differ], n_haps[differ]), collapse = "; "),
            ". The table shows the SNP run's n. Use the same popmap for both runs.", call. = FALSE)
  invisible(NULL)
}

## Not exported. "a", "a and b", "a, b and c".
.and_list <- function(x) {
  if (length(x) <= 1L) return(paste(x, collapse = ""))
  paste(paste(x[-length(x)], collapse = ", "), "and", x[length(x)])
}

## Not exported. The caption of a diversity_table(): which VCF each column came
## from, the rarefaction size g, and which standard error each column holds,
## all read from the settings of the runs.
.diversity_caption <- function(snps, haps) {
  per_site <- !is.null(snps) && !is.null(summary(snps)$tables$per_site)
  sources <- c(
    if (!is.null(snps))
      paste(if (per_site) "Ho, He and the per-site Ho and He are" else "Ho and He are",
            "from the SNP VCF"),
    if (!is.null(haps)) "Fis, Ar and rarefied private alleles are from the haplotype VCF")
  text <- paste0(paste(sources, collapse = "; "), ".")
  if (!is.null(haps)) {
    private_n <- unique(haps$richness$privAr_n)
    text <- paste0(text, sprintf(
      paste(" Ar and rarefied private alleles are counted in g = %d gene copies at each locus;",
            "the rarefied private alleles are the expected number found in only one",
            "population, added up over %s."),
      haps$settings$g,
      if (length(private_n) == 1L)
        sprintf("the %s loci typed at g copies in every population", .big(private_n))
      else "the loci typed at g copies in every population"))
  }
  text <- paste(text, "Values are estimate (SE).")

  ## Which SE each column holds (see .diversity_results()).
  combined <- function(x) .diversity_results(x)$combined
  both <- loci <- character(0)
  if (!is.null(snps)) {
    snp_stats <- c("Ho", "He", if (per_site) "the per-site values")
    if (combined(snps)) both <- c(both, snp_stats) else loci <- c(loci, snp_stats)
  }
  if (!is.null(haps)) {
    if (combined(haps)) both <- c(both, "Fis") else loci <- c(loci, "Fis")
    loci <- c(loci, "Ar", "rarefied private alleles")
  }
  says <- function(stats, what)
    if (length(stats) == 1L) sprintf("the SE for %s counts variation among %s", stats, what)
    else sprintf("SEs for %s count variation among %s", .and_list(stats), what)
  se <- c(if (length(both)) says(both, "RAD loci and among individuals"),
          if (length(loci)) says(loci, "RAD loci only"))
  se <- paste(se, collapse = "; ")
  paste0(text, " ", toupper(substr(se, 1, 1)), substring(se, 2), ".")
}
