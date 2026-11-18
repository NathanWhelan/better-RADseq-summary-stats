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
#' | `rarefied private SNPs (SE)` | SNP VCF | RAD loci only |
#' | `absolute private SNPs (no SE)` (only with `PA_count = "absolute"` or `"both"`) | SNP VCF | none: a count |
#' | `Fis (SE)` | haplotype VCF | RAD loci and individuals |
#' | `Ar (SE)`, `rarefied private alleles (SE)` | haplotype VCF | RAD loci only |
#' | `absolute private alleles (no SE)` (only with `PA_count = "absolute"` or `"both"`) | haplotype VCF | none: a count |
#'
#' With `se_individuals = FALSE` in [diversity_stats()], the SEs count RAD
#' loci only, and the caption says so.
#'
#' **Private SNPs and private alleles.** A private SNP allele (from the SNP
#' VCF) and a private allele (a haplotype, from the haplotype VCF) count
#' different things, so the table names each one. `private` picks which kinds
#' are shown and `PA_count` which count: the rarefied count (the expected
#' number in `g` gene copies, the fair comparison between populations), the
#' absolute count (Stacks' rule: what Stacks' `Private` column reports), or
#' both. The rarefied count is often higher than the absolute one; see
#' "Private alleles and Stacks' `Private` column" in [diversity_stats()].
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
#' @param private Which private alleles to show: `"snps"` (private SNPs, from
#'   `snps`), `"alleles"` (private alleles, from `haps`) or `"both"`. Default:
#'   every kind the runs given hold, so `"both"` when both runs are given.
#' @param PA_count Which count of each: `"rarefied"` (the default, with its
#'   SE), `"absolute"` (no SE) or `"both"`.
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
#' # Private alleles only, rarefied and absolute side by side:
#' diversity_table(res_snps, res_haps, private = "alleles", PA_count = "both")
#' out <- tempfile(fileext = ".csv")
#' write.csv(table1, out, row.names = FALSE)
#' @export
diversity_table <- function(snps = NULL, haps = NULL, private = NULL, PA_count = "rarefied") {
  if (is.null(snps) && is.null(haps))
    stop("Give at least one diversity_stats() result: `snps` (from the SNP VCF), ",
         "`haps` (from the haplotype VCF), or both.", call. = FALSE)
  .check_diversity_run(snps, "snps", haplotype = FALSE)
  .check_diversity_run(haps, "haps", haplotype = TRUE)
  .check_choice(PA_count, "PA_count", .pa_count_choices)
  ## By default, every kind of private allele the runs given hold.
  if (is.null(private))
    private <- if (!is.null(snps) && !is.null(haps)) "both" else if (!is.null(snps)) "snps" else "alleles"
  .check_choice(private, "private", c("snps", "alleles", "both"))
  show_snps <- private %in% c("snps", "both")
  show_alleles <- private %in% c("alleles", "both")
  if (show_snps && is.null(snps))
    stop("Private SNPs come from the SNP VCF: give `snps` (the diversity_stats() result of ",
         "populations.snps.vcf), or use private = \"alleles\".", call. = FALSE)
  if (show_alleles && is.null(haps))
    stop("Private alleles come from the haplotype VCF: give `haps` (the diversity_stats() ",
         "result of populations.haps.vcf), or use private = \"snps\".", call. = FALSE)

  first <- if (!is.null(snps)) snps else haps
  pops <- first$per_population$population
  if (!is.null(snps) && !is.null(haps)) .check_same_populations(snps, haps)

  ## The cells of each run's short summary, joined by population name; the
  ## private-allele columns of a kind not asked for are left out.
  without_private <- function(tab) tab[!grepl("private (SNPs|alleles) \\((no )?SE\\)$", names(tab))]
  pieces <- list()
  ## A run whose private columns are left out is asked only for the rarefied
  ## count, so a result saved by an older version (no priv_absolute) still
  ## works there.
  if (!is.null(snps)) {
    tables <- summary(snps, PA_count = if (show_snps) PA_count else "rarefied")$tables
    results <- if (show_snps) tables$results else without_private(tables$results)
    ## the per-site table goes right after Ho and He, before the private SNPs
    private_cols <- grepl("private SNPs", names(results))
    pieces <- c(pieces, list(results[!private_cols]),
                if (!is.null(tables$per_site)) list(tables$per_site),
                if (any(private_cols)) list(results[c("population", names(results)[private_cols])]))
  }
  if (!is.null(haps)) {
    results <- summary(haps, PA_count = if (show_alleles) PA_count else "rarefied")$tables$results
    if (!show_alleles) results <- without_private(results)
    if (!is.null(snps)) results$n <- NULL          # `n` is in the SNP table already
    pieces <- c(pieces, list(results))
  }
  out <- data.frame(population = pops, stringsAsFactors = FALSE)
  for (piece in pieces)
    out <- cbind(out, piece[match(pops, piece$population), setdiff(names(piece), "population"),
                            drop = FALSE])
  rownames(out) <- NULL
  structure(out, caption = .diversity_caption(snps, haps, show_snps, show_alleles, PA_count),
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
## from, the rarefaction size g, what the private-allele columns count, and
## which standard error each column holds, all read from the settings of the
## runs. `show_snps` and `show_alleles` say whether the private SNP and
## private allele columns are in the table; `PA_count` which counts.
.diversity_caption <- function(snps, haps, show_snps = !is.null(snps),
                               show_alleles = !is.null(haps), PA_count = "rarefied") {
  per_site <- !is.null(snps) && !is.null(summary(snps)$tables$per_site)
  rarefied <- PA_count != "absolute"
  absolute <- PA_count != "rarefied"
  snp_stats <- if (per_site && show_snps) "Ho, He, the per-site Ho and He, and private SNPs"
               else if (per_site) "Ho, He and the per-site Ho and He"
               else if (show_snps) "Ho, He and private SNPs"
               else "Ho and He"
  sources <- c(
    if (!is.null(snps)) paste(snp_stats, "are from the SNP VCF"),
    if (!is.null(haps))
      paste(.and_list(c("Fis", "Ar", if (show_alleles) "private alleles")),
            "are from the haplotype VCF"))
  text <- paste0(paste(sources, collapse = "; "), ".")

  ## What the rarefied counts are added up over: the same words as the short
  ## summary's legend (.private_over()).
  if (show_snps && rarefied)
    text <- paste0(text, sprintf(
      paste(" Rarefied private SNPs are the expected number of SNP alleles found in only",
            "one population in g = %d gene copies at each SNP, added up over %s."),
      snps$settings$g, .private_over(snps)))
  if (!is.null(haps)) {
    text <- paste0(text, if (show_alleles && rarefied) sprintf(
      paste(" Ar and rarefied private alleles are counted in g = %d gene copies at each locus;",
            "the rarefied private alleles are the expected number found in only one",
            "population, added up over %s."), haps$settings$g, .private_over(haps))
      else sprintf(" Ar is counted in g = %d gene copies at each locus.", haps$settings$g))
  }
  absolute_kinds <- c(if (show_snps) "private SNPs", if (show_alleles) "private alleles")
  if (absolute && length(absolute_kinds))
    text <- paste0(text, sprintf(
      paste(" Absolute %s count the alleles found in only one population, counting every",
            "individual (as Stacks does); they have no SE."), .and_list(absolute_kinds)))
  text <- paste(text, "Values are estimate (SE).")

  ## Which SE each column holds (see .diversity_results()).
  combined <- function(x) .diversity_results(x)$combined
  both <- loci <- character(0)
  if (!is.null(snps)) {
    snp_stats <- c("Ho", "He", if (per_site) "the per-site values")
    if (combined(snps)) both <- c(both, snp_stats) else loci <- c(loci, snp_stats)
    if (show_snps && rarefied) loci <- c(loci, "rarefied private SNPs")
  }
  if (!is.null(haps)) {
    if (combined(haps)) both <- c(both, "Fis") else loci <- c(loci, "Fis")
    loci <- c(loci, "Ar", if (show_alleles && rarefied) "rarefied private alleles")
  }
  says <- function(stats, what)
    if (length(stats) == 1L) sprintf("the SE for %s counts variation among %s", stats, what)
    else sprintf("SEs for %s count variation among %s", .and_list(stats), what)
  se <- c(if (length(both)) says(both, "RAD loci and among individuals"),
          if (length(loci)) says(loci, "RAD loci only"))
  se <- paste(se, collapse = "; ")
  paste0(text, " ", toupper(substr(se, 1, 1)), substring(se, 2), ".")
}
