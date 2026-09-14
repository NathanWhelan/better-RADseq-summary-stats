###############################################################################
#
#  R/diversity_stats.R -- Ho, He, FIS, rarefied allelic richness and rarefied
#  private allelic richness from a Stacks VCF, in one run.
#
#  WHY THE NUMBERS DIFFER FROM STACKS, in one line each. The full explanation,
#  with every Stacks column matched to its counterpart here, is in
#  vignette("rationale").
#
#    He     Nei & Chesser (1983), not Stacks' `Pi`, whose correction assumes
#           FIS = 0. Both are reported side by side.     [-> section 2]
#    FIS    a ratio of sums over loci, not Stacks' mean of per-locus ratios,
#           and not divided by Pi.                       [-> section 2]
#    CIs    a bootstrap over RAD LOCI (boot = "loci", the default), not over
#           SNP rows, because SNPs on one RAD tag are linked. boot =
#           "individuals"/"both" are for comparison only (see
#           vignette("rationale"), "Why not bootstrap individuals?"). None of
#           these intervals is a test between populations; that is
#           het_between_pops().
#    records  by default a population uses a record once it has >= min_n
#           genotyped individuals there (available data, per population).
#           complete_case = TRUE uses only records genotyped in every
#           individual of every population (Schmidt et al. 2021).
#
#  This file holds the exported function (an outline of the steps), its
#  print() and summary() methods, and the handling of the `sites` argument.
#  The steps themselves are in R/diversity_internals.R and the estimators in
#  R/estimators.R. Checked by tests/testthat/test-diversity-stats.R,
#  test-regression.R (golden values) and diversity_core_selftest().
#
###############################################################################

#' Diversity statistics from a Stacks VCF
#'
#' Computes Ho, He, FIS (Nei & Chesser 1983, ratio of sums), rarefied allelic
#' richness and rarefied private allelic richness from a Stacks 2
#' `populations` VCF and popmap, with a delete-one-locus jackknife standard
#' error and a bootstrap confidence interval over RAD loci. Printing the
#' result shows the main tables; `summary()` shows the full report with notes
#' on interpretation; `outdir` also writes the tables to TSV files.
#'
#' Run it twice: `populations.snps.vcf` gives Ho, He, `pct_poly` and the
#' per-sequenced-site values; `populations.haps.vcf` gives FIS, Ar and
#' privAr. The printed result says which numbers to take from each file. See
#' `vignette("rationale")`, section 3.
#'
#' The SNP-VCF side was designed for biallelic SNPs, as Stacks writes them.
#' A VCF whose alleles are all single bases is treated as one site per record
#' even if a few SNPs are multi-allelic (as GATK or freebayes can emit), and
#' those sites are counted correctly; a VCF with any multi-base allele is
#' treated as haplotype data (see "Assumptions" in [read_stacks_vcf()]).
#'
#' @param vcf Path to `populations.snps.vcf` or `populations.haps.vcf`
#'   (optionally gzip-compressed), or the object returned by
#'   [read_stacks_vcf()], optionally passed through `filter_*()` functions.
#' @param popmap Path to a two-column, no-header popmap file (`sample_id
#'   <TAB> population`), or the list returned by [read_popmap()].
#' @param g Rarefaction size in GENE COPIES (10 diploid individuals = 20).
#'   Required; at most twice the size of the smallest population.
#' @param nboot Bootstrap replicates. Default `10000`; `0` for none (the
#'   jackknife standard errors are always computed).
#' @param boot What the bootstrap resamples: `"loci"` (whole RAD loci, the
#'   default), `"individuals"` (within each population, as diveRsity does) or
#'   `"both"` (Owen & Eckles 2012). Matched exactly. `"individuals"` and
#'   `"both"` are for comparison only: in a coverage simulation their
#'   intervals badly missed the true He and FIS (`vignette("rationale")`,
#'   "Why not bootstrap individuals?"). For uncertainty over individuals use
#'   `se_individuals = TRUE` instead.
#' @param sites Total sequenced sites, to convert Ho and He to per-site values
#'   (`Ho_autosomal`, `He_autosomal`; SNP VCF only). Default `NULL`: no
#'   conversion. One of:
#'   * one number, used for every population;
#'   * a named numeric vector, one value per population (names must match
#'     the popmap's populations);
#'   * a data frame with a `population` (or `pop`) column and a `sites`
#'     column;
#'   * the path to Stacks' `populations.sumstats_summary.tsv`, from which each
#'     population's own `Sites` ("All positions (variant and fixed)" block) is
#'     read with [read_sumstats_summary()].
#'
#'   A population's own value matters because Stacks counts a site for a
#'   population only where that population has a genotyped individual. For
#'   the same reason, each population's He is multiplied by its own number of
#'   variant records: `Variant_Sites` from the same file when `sites` is that
#'   file, otherwise the records of this VCF where the population has a
#'   genotyped individual (column `variant_records` of the result).
#'
#'   Per-site values assume the variant records are all the variant sites
#'   inside `sites`. Records removed in R with a `filter_*()` function are
#'   still counted in `sites`: pass the `populations.sumstats_summary.tsv`
#'   path (whose `Variant_Sites` is unfiltered) and a warning is given.
#'   Allele-frequency filters ([filter_maf()], [filter_mac()], Stacks'
#'   `--min-mac`/`--min-maf`) remove rare variants, which lowers per-site
#'   values by an amount that grows as samples get smaller; report the
#'   threshold used (see `vignette("rationale")`, section 1).
#' @param min_n Minimum genotyped individuals a population needs at a record
#'   to use that record. Default `2`, the smallest number for which He and
#'   FIS are defined. Ignored when `complete_case = TRUE`.
#' @param complete_case Use a record only if every individual of every
#'   population is genotyped there (Schmidt et al. 2021, recommendation (b)).
#'   Default `FALSE`.
#' @param se_individuals Also report delete-one-INDIVIDUAL jackknife standard
#'   errors (`Ho_se_ind`, `He_se_ind`, `Fis_se_ind`, `Ar_se_ind`,
#'   `privAr_se_ind`): the uncertainty from which individuals were sampled,
#'   which the locus-based `_se`, `_lo` and `_hi` leave out. Report these when
#'   individuals differ in inbreeding or include relatives (check with
#'   [identity_disequilibrium()]: g2 > 0). The Ar/privAr ones need
#'   `g <= 2 * (n - 1)` for every population and are `NA` otherwise. Two
#'   caveats for those two, from simulation (`vignette("rationale")`, section
#'   4): `privAr_se_ind` deletes only that population's individuals and holds
#'   the other populations' fixed, although private richness depends on them
#'   too, so it runs too small (by up to about 20%); and when many records
#'   have fewer than `g + 2` gene copies (missing data with `g` close to
#'   `2 * (n - 1)`), deleting one individual drops them from some jackknife
#'   replicates, and both come out too large (up to about twice the true
#'   value); a warning then suggests a smaller `g` for them. Default `FALSE`.
#' @param hierfstat_check If `TRUE` and the `hierfstat` package is installed,
#'   also compute per-record Ho/Hs with `hierfstat::basic.stats()` and
#'   allelic richness with `hierfstat::allelic.richness()`, and report the
#'   largest difference from this package's values (with a warning if they
#'   disagree). Default `FALSE`: the package tests already make this
#'   comparison, and on large datasets the hierfstat calls are slow.
#' @param outdir Directory to write the result tables to as TSV files
#'   (created if needed). Default `NULL`: write nothing.
#' @param stem Text in the output file names
#'   (`diversity_per_population.<stem>.tsv`). By default it comes from the
#'   VCF's file name (`populations.haps.vcf` gives `"haps"`), so runs on the
#'   haplotype and SNP VCFs can share `outdir` without overwriting each other.
#'   Required when writing files from an already-read object, which has no
#'   file name.
#' @param seed Random seed for the bootstrap. Default `NULL`: use R's current
#'   random-number stream, so call `set.seed()` first for reproducible
#'   intervals. A number makes the result reproducible on its own and leaves
#'   your session's random-number stream as it was.
#' @param verbose Print progress messages. Default `TRUE`.
#' @return An object of class `raddiv_diversity`, a list of:
#'   \describe{
#'     \item{per_population}{Ho, He, FIS and % polymorphic records per
#'       population, each with its jackknife SE (`_se`) and 95% bootstrap
#'       interval (`_lo`, `_hi`).}
#'     \item{richness}{Ar, privAr and priv_total per population, with the
#'       number of records each rests on (`Ar_n`, `privAr_n`).}
#'     \item{autosomal}{Per-sequenced-site Ho/He, with the `sites_used` and
#'       `variant_records` each population was scaled by; `NULL` unless
#'       `sites` gave at least one population a plausible value, and always
#'       `NULL` for a haplotype VCF.}
#'     \item{estimator_comparison}{He and FIS from Nei & Chesser and from the
#'       estimator behind Stacks' `Pi`, side by side.}
#'     \item{he_difference}{For two populations with `nboot > 0`: the He
#'       difference with its locus-bootstrap interval (not a test; see
#'       [het_between_pops()]). Otherwise `NULL`.}
#'     \item{fis_by_call_rate}{A check for null alleles and allele dropout:
#'       each population's records grouped by that population's call rate
#'       (`100%`, `90-99%`, `75-89%`, `<75%`), with the number of records, He,
#'       FIS and FIS's locus-jackknife SE in each group. FIS that rises as
#'       call rate falls points to null alleles or dropout rather than
#'       inbreeding (Gautier et al. 2013).}
#'     \item{settings}{The settings and record counts of this run.}
#'   }
#'   Values are stored at full precision; `print()`, `summary()` and the TSV
#'   files round them. With `outdir`, `per_population`, `richness` and
#'   `autosomal` are written to `diversity_<table>.<stem>.tsv`.
#' @references
#' Nei, M. & Chesser, R.K. (1983) Estimation of fixation indices and gene
#' diversities. *Annals of Human Genetics* 47:253-259.
#' \doi{10.1111/j.1469-1809.1983.tb00993.x}
#'
#' Weir, B.S. & Cockerham, C.C. (1984) Estimating F-statistics for the
#' analysis of population structure. *Evolution* 38:1358-1370.
#' \doi{10.1111/j.1558-5646.1984.tb05657.x}
#'
#' Hurlbert, S.H. (1971) The nonconcept of species diversity: a critique and
#' alternative parameters. *Ecology* 52:577-586.
#'
#' El Mousadik, A. & Petit, R.J. (1996) High level of genetic differentiation
#' for allelic richness among populations of the argan tree. *Theoretical and
#' Applied Genetics* 92:832-839.
#'
#' Kalinowski, S.T. (2004) Counting alleles with rarefaction: private alleles
#' and hierarchical sampling designs. *Conservation Genetics* 5:539-543.
#' \doi{10.1023/B:COGE.0000041021.91777.1a}
#'
#' Szpiech, Z.A., Jakobsson, M. & Rosenberg, N.A. (2008) ADZE: a rarefaction
#' approach for counting alleles private to combinations of populations.
#' *Bioinformatics* 24:2498-2504.
#'
#' Schmidt, T.L., Jasper, M.-E., Weeks, A.R. & Hoffmann, A.A. (2021) Unbiased
#' population heterozygosity estimates from genome-wide sequence data.
#' *Methods in Ecology and Evolution* 12:1888-1898.
#' \doi{10.1111/2041-210X.13659}
#'
#' Weir, B.S. (1996) *Genetic Data Analysis II*. Sinauer, Sunderland, MA.
#' (Delete-one jackknife over loci.)
#'
#' Owen, A.B. & Eckles, D. (2012) Bootstrapping data arrays of arbitrary
#' order. *Annals of Applied Statistics* 6:895-927. (`boot = "both"`.)
#'
#' Keenan, K., McGinnity, P., Cross, T.F., Crozier, W.W. & Prodohl, P.A.
#' (2013) diveRsity: an R package for the estimation and exploration of
#' population genetics parameters and their associated errors. *Methods in
#' Ecology and Evolution* 4:782-788. (`boot = "individuals"`.)
#'
#' Goudet, J. (2005) HIERFSTAT, a package for R to compute and test
#' hierarchical F-statistics. *Molecular Ecology Notes* 5:184-186.
#' (`hierfstat_check`.)
#'
#' Gautier, M., Gharbi, K., Cezard, T., et al. (2013) The effect of RAD
#' allele dropout on the estimation of genetic variation within and between
#' populations. *Molecular Ecology* 22:3165-3178. (`fis_by_call_rate`.)
#' @examples
#' # A toy dataset shipped with the package: 80 RAD loci, 2 populations of 4
#' # and 3 individuals (real studies need far more individuals).
#' haps   <- system.file("extdata", "small.haps.vcf", package = "RADdiversity")
#' snps   <- system.file("extdata", "small.snps.vcf", package = "RADdiversity")
#' popmap <- system.file("extdata", "small_popmap.tsv", package = "RADdiversity")
#'
#' # Haplotype VCF: FIS, allelic richness, private allelic richness.
#' set.seed(1)
#' res <- diversity_stats(haps, popmap, g = 4, nboot = 100)
#' res                    # the main tables
#' res$per_population     # one table
#'
#' # SNP VCF with each population's own sequenced-site count, read from Stacks'
#' # populations.sumstats_summary.tsv: Ho and He per sequenced site.
#' sumstats <- system.file("extdata", "populations.sumstats_summary.tsv",
#'                         package = "RADdiversity")
#' res_snps <- diversity_stats(snps, popmap, g = 4, nboot = 0, sites = sumstats,
#'                             verbose = FALSE)
#' res_snps$autosomal
#'
#' # The full report, with notes on interpreting each number:
#' summary(res_snps)
#' @export
diversity_stats <- function(vcf, popmap, g, nboot = 10000L, boot = "loci", sites = NULL,
                            min_n = 2L, complete_case = FALSE, se_individuals = FALSE,
                            hierfstat_check = FALSE, outdir = NULL, stem = NULL,
                            seed = NULL, verbose = TRUE) {

  ## ---- 1. Check the arguments (before any slow reading) --------------------
  if (missing(g))
    stop("`g` is missing: give the rarefaction size in gene copies, e.g. g = 20 ",
         "for 10 diploid individuals.", call. = FALSE)
  g <- .check_count(g, "g", min = 1)
  if (g < 2) stop("g must be at least 2 gene copies.", call. = FALSE)
  nboot <- .check_count(nboot, "nboot")
  .check_choice(boot, "boot", c("loci", "individuals", "both"))
  min_n <- .check_count(min_n, "min_n")
  if (min_n < 2)
    stop("min_n must be at least 2 (He and FIS are undefined below n = 2).", call. = FALSE)
  .check_flag(complete_case, "complete_case")
  .check_flag(se_individuals, "se_individuals")
  .check_flag(hierfstat_check, "hierfstat_check")
  .check_flag(verbose, "verbose")
  .check_seed(seed)
  .resolve_sites(sites)      # shape only; matched to populations in step 3
  .check_run_inputs(vcf, popmap, stem, outdir)

  if (!is.null(seed)) {
    restore_rng <- .save_rng_state()
    on.exit(restore_rng(), add = TRUE)
    set.seed(seed)
  }
  use_hierfstat <- hierfstat_check && .hierfstat_available()
  if (hierfstat_check && !use_hierfstat)
    .inform(verbose, "hierfstat_check = TRUE but hierfstat is not installed; skipping the cross-check.")

  ## ---- 2. Read the data -----------------------------------------------------
  H <- .resolve_H(vcf, verbose = verbose)
  pops <- .resolve_pops(popmap, H$samples, verbose = verbose)
  pop_names <- names(pops)
  pop_sizes <- lengths(pops)
  n_rec <- nrow(H$A1)

  ## ---- 3. Checks that need the populations ----------------------------------
  if (length(pops) < 2) stop("Need at least 2 populations.", call. = FALSE)
  ## With one individual, He is undefined and every statistic would read as 0.
  tiny <- pop_names[pop_sizes < 2]
  if (length(tiny))
    stop("Population(s) with fewer than 2 individuals: ", paste(tiny, collapse = ", "),
         "\n  Expected heterozygosity is undefined at n = 1, and rarefaction to ",
         "2 gene\n  copies is meaningless. Drop these populations or merge them.", call. = FALSE)
  if (g > 2L * min(pop_sizes))
    stop("g = ", g, " gene copies exceeds the smallest population's ",
         2L * min(pop_sizes), ". g is in GENE COPIES: 10 diploids = 20.", call. = FALSE)
  sites_by_pop <- .resolve_sites(sites, pop_names)

  ## Data may have been filtered before reaching this package; say what it
  ## looks filtered at (see .prior_filter_signals() in R/filter_loci.R).
  prior_filters <- .prior_filter_signals(H, pops)

  is_haplotype <- .is_haplotype_H(H)
  .inform(verbose, sprintf("  record type: %s", if (is_haplotype)
    "HAPLOTYPE (one multi-allelic locus per RAD tag)" else "SNP (one site per record)"))
  .inform(verbose, sprintf("  %s records on %s RAD loci (%.2f per locus)", .big(n_rec),
                           .big(length(unique(H$locus_raw))),
                           n_rec / length(unique(H$locus_raw))))

  ## ---- 4. Which records are used --------------------------------------------
  n_typed_all <- .typed_by_pop(H, pops)
  inclusion <- .locus_inclusion(n_typed_all, pop_sizes, min_n, complete_case, verbose)
  rows <- which(inclusion$keep)
  ## Under complete_case every used record is fully genotyped, so only the
  ## mathematical floor of 2 individuals applies.
  min_n_used <- if (complete_case) 2L else min_n
  ## The RAD locus of each used record, as a number 1..n_loci. Standard errors
  ## and the bootstrap resample these loci.
  locus_names <- unique(H$locus_raw[rows])
  locus_index <- match(H$locus_raw[rows], locus_names)
  n_loci <- length(locus_names)

  ## ---- 5. Per-record statistics ---------------------------------------------
  .inform(verbose, "Engine: internal (Nei & Chesser Hs, R/estimators.R)",
          if (use_hierfstat) "; hierfstat cross-check requested" else "")
  .inform(verbose, sprintf("Rarefying to g = %d gene copies (smallest population has %d)",
                           g, 2L * min(pop_sizes)))
  counts <- .pop_counts(H, pops, rows)
  n_typed <- n_typed_all[rows, , drop = FALSE]
  n_het <- .het_by_pop(H, pops, rows)
  rs <- .record_stats(counts, n_typed, n_het, min_n_used, g)
  if (use_hierfstat) .hierfstat_crosscheck(H, pops, rows, rs, n_typed, g, min_n_used, verbose)

  n_undefined <- sum(!(is.finite(rs$Ho) & is.finite(rs$Hs)))
  if (n_undefined)
    .inform(verbose, sprintf("  %s record-by-population cells have no defined Ho/He and are ",
                             .big(n_undefined)),
            "excluded from that population's means (not counted as zero).")
  if (!any(rs$n_present > 1))
    stop("No locus is polymorphic in any population, so He, FIS and allelic\n",
         "  richness are all degenerate. Check that genotypes were parsed: the\n",
         "  missing-genotype rate reported above should not be near 100%.", call. = FALSE)

  ## ---- 6. Point estimates and jackknife over RAD loci -----------------------
  ## One row per RAD locus of summed pieces (see .diversity_pieces()).
  S <- rowsum(.diversity_pieces(rs), locus_index, reorder = TRUE)
  stats_from_sums <- function(totals) .diversity_from_sums(totals, pop_names)
  point <- stats_from_sums(colSums(S))[1L, ]
  ## Delete-one-locus jackknife: deterministic, needs no nboot (Weir 1996).
  se_loci <- .jack_block_sums(S, stats_from_sums)

  ## ---- 7. Optional: jackknife over individuals ------------------------------
  se_ind <- NULL
  jackknife_boundary <- NULL
  if (se_individuals) {
    .inform(verbose, "Delete-one-individual jackknife over ", .big(sum(pop_sizes)),
            " individuals ...")
    jack <- .jack_individuals(H, pops, rows, counts, n_typed, n_het, min_n_used, g)
    se_ind <- jack$se
    ## Built-in check that the two code paths agree.
    gap <- suppressWarnings(max(abs(jack$full - point[names(jack$full)]), na.rm = TRUE))
    if (is.finite(gap) && gap > 1e-8)
      warning("Individual jackknife: its full-sample estimates differ from the point ",
              "estimates by ", signif(gap, 3), ". Please report this.", call. = FALSE)
    if (g > 2L * (min(pop_sizes) - 1L))
      .inform(verbose, "  Ar/privAr individual SEs are NA for populations with fewer than ",
              "g/2 + 1 individuals: removing one leaves fewer than g = ", g,
              " gene copies. Use g <= ", 2L * (min(pop_sizes) - 1L), " to get them ",
              "(with missing data, a smaller g still; see the warning if one is given).")
    .inform(verbose, "  Note: privAr_se_ind holds the OTHER populations' individuals fixed, so it ",
            "leaves out part of the uncertainty (up to about 20% too small in simulation; ",
            "see ?diversity_stats).")
    ## Ar/privAr SEs are inflated when many records sit within 2 gene copies
    ## of g (see .jackknife_boundary()). Only populations whose SEs exist.
    boundary <- .jackknife_boundary(n_typed, rs$Ar, g, min_n_used, .jackknife_boundary_share)
    has_se <- is.finite(se_ind[paste0("Ar_", pop_names)])
    jackknife_boundary <- boundary$share
    inflated <- pop_names[has_se & is.finite(boundary$share) &
                            boundary$share > .jackknife_boundary_share]
    if (length(inflated))
      warning("Ar_se_ind and privAr_se_ind are probably too large for ",
              paste(sprintf("%s (%.0f%% of its records with a defined Ar)", inflated,
                            100 * boundary$share[inflated]), collapse = ", "),
              ": those records have fewer than g + 2 = ", g + 2L, " gene copies, so deleting one ",
              "individual drops them from some jackknife replicates, which made these SEs up ",
              "to about twice the true value in simulation (vignette(\"rationale\"), section 4). ",
              if (!is.na(boundary$suggest_g) && boundary$suggest_g >= 2L)
                paste0("For these two SEs, re-run with g = ", boundary$suggest_g, " or less. ")
              else "For these two SEs, re-run with a smaller g. ",
              "The Ho, He and FIS SEs are not affected.", call. = FALSE)
  }

  ## ---- 8. Bootstrap ---------------------------------------------------------
  boot_matrix <- NULL
  if (nboot > 0) {
    .inform(verbose, sprintf("Bootstrapping %s replicates over %s (boot = \"%s\") ...",
                             .big(nboot), switch(boot,
      loci        = paste(.big(n_loci), "RAD loci"),
      individuals = "individuals within each population",
      both        = paste(.big(n_loci), "RAD loci AND individuals within each population")),
      boot))
    boot_matrix <- switch(boot,
      loci        = .boot_block_sums(S, nboot, stats_from_sums),
      individuals = .boot_individuals(H, pops, rows, locus_index, n_loci, min_n_used, g,
                                      nboot, resample_loci = FALSE),
      both        = .boot_individuals(H, pops, rows, locus_index, n_loci, min_n_used, g,
                                      nboot, resample_loci = TRUE))
  }
  ci <- .percentile_ci(boot_matrix, names(point))

  ## ---- 9. Result tables ------------------------------------------------------
  ## Per-site values divide by each population's own sequenced sites and
  ## multiply by its own count of variant records (see .autosomal_counts()).
  ## A sites value is plausible if it is at least that count.
  counts_autosomal <- .autosomal_counts(n_typed_all, sites_by_pop, H$n_records_read, is_haplotype)
  variant_records <- counts_autosomal$variant_records
  for (note_text in counts_autosomal$notes) warning(note_text, call. = FALSE)
  ok_sites <- is.finite(sites_by_pop) & sites_by_pop > 0 & sites_by_pop >= variant_records
  tables <- .diversity_tables(pop_names, pop_sizes, point, se_loci, ci, se_ind,
                              ar_n = colSums(!is.na(rs$Ar)), pr_n = colSums(!is.na(rs$Pr)),
                              sites = sites_by_pop, ok_sites = ok_sites,
                              variant_records = variant_records,
                              is_haplotype = is_haplotype, boot_matrix = boot_matrix)
  tables$fis_by_call_rate <- .fis_by_call_rate(rs, n_typed, pop_sizes, locus_index)

  ## ---- 10. Files (only with outdir) -----------------------------------------
  ## A non-default boot mode is in the file name, so a comparison run cannot
  ## overwrite the default output.
  written <- character(0)
  if (!is.null(outdir)) {
    stem <- .derive_stem(vcf, stem)
    boot_suffix <- if (boot != "loci") paste0(".boot-", boot) else ""
    to_write <- c("per_population", "richness", "autosomal")
    written <- .write_tables(
      lapply(stats::setNames(to_write, to_write),
             function(nm) .round_diversity_table(tables[[nm]], nm)),
      outdir,
      stats::setNames(sprintf("diversity_%s.%s%s.tsv", to_write, stem, boot_suffix), to_write))
  }

  settings <- list(n_pops = length(pops), n_individuals = pop_sizes, n_records = n_rec,
                   n_records_used = length(rows), n_loci = n_loci, g = g,
                   complete_case = complete_case, min_n = min_n,
                   cells_used = inclusion$cells_used, cells_total = inclusion$cells_total,
                   nboot = nboot, boot = boot, se_individuals = se_individuals,
                   jackknife_boundary = jackknife_boundary,
                   is_haplotype = is_haplotype, sites = as.vector(sites_by_pop),
                   variant_records = variant_records, ok_sites = ok_sites,
                   prior_filters = prior_filters, seed = seed, files = written)
  names(settings$sites) <- pop_names
  structure(c(tables, list(settings = settings)), class = "raddiv_diversity")
}

#' @rdname diversity_stats
#' @param x,object A `raddiv_diversity` object, as returned by
#'   `diversity_stats()`.
#' @param ... Ignored.
#' @export
print.raddiv_diversity <- function(x, ...) {
  st <- x$settings
  cat(sprintf("Diversity statistics: %d populations, %s records on %s RAD loci (%s VCF)\n",
              st$n_pops, .big(st$n_records_used), .big(st$n_loci),
              if (st$is_haplotype) "haplotype" else "SNP"))
  cat(sprintf("  He: Nei & Chesser (1983); FIS: 1 - sum(Ho)/sum(He); Ar, privAr: rarefied to g = %d gene copies\n",
              st$g))
  cat("  _se: jackknife over RAD loci",
      if (st$nboot > 0) sprintf("; _lo/_hi: 95%% bootstrap interval, %s replicates (boot = \"%s\")",
                                .big(st$nboot), st$boot) else "; no bootstrap (nboot = 0)",
      "\n", sep = "")

  cat("\n$per_population\n")
  .print_table(.round_diversity_table(x$per_population, "per_population"))
  cat("\n$richness\n")
  .print_table(.round_diversity_table(x$richness, "richness"))
  if (!is.null(x$autosomal)) {
    cat("\n$autosomal\n")
    .print_table(.round_diversity_table(x$autosomal, "autosomal"))
  }

  cat("\n")
  if (st$is_haplotype) {
    cat("Take from this haplotype VCF: Fis, Ar, privAr. (Ho and He: run the SNP VCF.)\n")
  } else {
    cat("Take from this SNP VCF: Ho, He, pct_poly", if (!is.null(x$autosomal)) ", Ho_autosomal, He_autosomal",
        ". (Fis, Ar, privAr: run the haplotype VCF.)\n", sep = "")
  }
  for (w in .diversity_warnings(x)) cat("NOTE:", w, "\n")
  cat("summary() prints the full report, with notes on interpreting each number.\n")
  invisible(x)
}

## Not exported. The short warnings print.raddiv_diversity() shows; the full
## report in summary() explains each of them.
.diversity_warnings <- function(x) {
  st <- x$settings
  out <- character(0)
  if (st$boot != "loci")
    out <- c(out, sprintf("boot = \"%s\" is a comparison mode; its intervals undercover (see ?diversity_stats).", st$boot))
  if (min(x$richness$privAr_n) < 0.5 * st$n_records_used)
    out <- c(out, "privAr rests on under half of the records for some population; consider a smaller g.")
  if (st$complete_case && st$n_records_used / st$n_records < 0.5)
    out <- c(out, "complete_case dropped most records for missing data; check poor libraries.")
  if (any(st$ok_sites) && st$is_haplotype)
    out <- c(out, "`sites` was ignored: per-site values need the SNP VCF.")
  if (!any(st$ok_sites) && any(st$sites > 0))
    out <- c(out, "`sites` was ignored: every value is smaller than the number of variant records.")
  if (isTRUE(st$se_individuals) && anyNA(x$richness$Ar_se_ind))
    out <- c(out, "some Ar/privAr individual SEs are NA: g is too large to remove one individual.")
  if (isTRUE(st$se_individuals)) {
    share <- st$jackknife_boundary
    inflated <- names(share)[is.finite(share) & share > .jackknife_boundary_share &
                               is.finite(x$richness$Ar_se_ind)]
    if (length(inflated))
      out <- c(out, sprintf("Ar_se_ind/privAr_se_ind are probably too large for %s (records within 2 gene copies of g); use a smaller g for them.",
                            paste(inflated, collapse = ", ")))
    out <- c(out, "privAr_se_ind holds the other populations' individuals fixed (up to about 20% too small in simulation).")
  }
  out
}

#' @rdname diversity_stats
#' @export
summary.raddiv_diversity <- function(object, ...) {
  structure(list(result = object), class = "summary.raddiv_diversity")
}

#' @export
print.summary.raddiv_diversity <- function(x, ...) {
  .diversity_report(x$result)
  invisible(x)
}

## Not exported. The full diversity report (summary() of a result, and what
## the command-line script prints). Numbers come from the result object; the
## longer fixed paragraphs are in .report_text (R/report_text.R).
.diversity_report <- function(x) {
  st <- x$settings
  n_used <- st$n_records_used
  n_rec <- st$n_records
  note <- function(...) cat(paste0("  ", c(...), "\n"), sep = "")
  see  <- function(section) cat(paste0("  -> vignette(\"rationale\"): \"", section, "\"\n"))
  rule <- function(char = "=") cat(strrep(char, 69), "\n", sep = "")

  cat("\n")
  rule()
  cat("  DIVERSITY --", st$n_pops, "populations,", .big(n_used), "records on",
      .big(st$n_loci), "RAD loci\n")
  cat("  He = Nei & Chesser (1983) Hs, the hierfstat quantity\n")
  cat("  FIS = 1 - sum(Ho)/sum(Hs), a ratio of sums, monomorphic loci included\n")
  cat("  Ar / private Ar rarefied to", st$g, "GENE COPIES (Kalinowski 2004)\n")
  if (st$complete_case) {
    cat(sprintf("  %s of %s variant records have complete data and were used (%.1f%%)\n",
                .big(n_used), .big(n_rec), 100 * n_used / n_rec))
    if (n_used / n_rec < 0.5) note(.report_text$diversity_complete_case_dropped)
  } else {
    cat(sprintf("  %s of %s record-by-population cells meet min_n = %d (%.1f%%)\n",
                .big(st$cells_used), .big(st$cells_total), st$min_n,
                100 * st$cells_used / st$cells_total))
    note(.report_text$diversity_available_data)
  }
  if (st$nboot > 0)
    cat("  95% CI from ", .big(st$nboot), " replicates, boot=\"", st$boot,
        "\" (see vignette(\"rationale\"), \"Why not bootstrap individuals?\")\n", sep = "")
  rule()
  cat("\n")

  .print_table(.round_diversity_table(x$per_population, "per_population"))
  note(.report_text$diversity_se_columns)
  if (isTRUE(st$se_individuals)) note(.report_text$diversity_se_ind_columns)
  cat("\n")

  .print_table(.round_diversity_table(x$richness, "richness"))
  note("_se columns: delete-one-block jackknife over RAD loci, same as above.",
       paste("Ar_n / privAr_n: how many of the", .big(n_used),
             "loci had a defined value for that population"),
       .report_text$diversity_ar_n)
  if (min(x$richness$privAr_n) < 0.5 * n_used)
    note(sprintf("WARNING: privAr_n is below 50%% of loci for at least one population --"),
         "its rarefied private richness rests on a minority of loci. Consider a",
         sprintf("smaller g (rarefaction target), currently %d gene copies.", st$g))
  note("Ar and privAr are PER LOCUS (the HP-RARE / ADZE convention);",
       paste("priv_total is the dataset-wide sum over privAr_n loci. Private = absent from all",
             st$n_pops - 1, "others,"),
       "so it is not comparable to a dataset with a different number of populations.")
  ar <- x$richness$Ar
  if (all(is.na(ar))) {
    note(sprintf("NOTE: Ar is undefined for every population -- no locus has any population at >= g = %d", st$g),
         "gene copies (see Ar_n above). That is too little data at this rarefaction target,",
         "not a biallelic-SNP ceiling; try a smaller g.")
  } else if (max(ar, na.rm = TRUE) <= 2.001) {
    note(.report_text$diversity_ar_capped)
  }

  cat("\nTwo estimators of the same parameter (gene diversity per variant record)\n")
  note(.report_text$diversity_two_estimators)
  .print_table(.round_diversity_table(x$estimator_comparison, "estimator_comparison"))
  n_ind <- st$n_individuals
  note(.report_text$diversity_fis_sensitivity,
       sprintf("Compare against `Pi`, never `Exp_Het`: Pi = Exp_Het * 2n/(2n-1) = %s.",
               paste(sprintf("x%.3f at n=%d", (2 * n_ind) / (2 * n_ind - 1), n_ind),
                     collapse = ", ")))
  see("Appendix A. Stacks columns and RADdiversity")

  pf <- st$prior_filters
  if (!is.null(pf)) {
    cat("\nWhat the data look filtered at (before or during this run)",
        if (!is.null(pf$n_samples)) sprintf(", over the %d individuals in the popmap", pf$n_samples),
        "\n", sep = "")
    note(sprintf("rarest allele in any variable record: %s copies; %.1f%% of %s variable records have an allele seen once or twice",
                 pf$min_allele_count, 100 * pf$rare_share, .big(pf$n_variable)),
         sprintf("lowest record call rate: %.2f pooled; within populations %s",
                 pf$min_call_pooled,
                 paste(sprintf("%s %.2f", names(pf$min_call_by_pop), pf$min_call_by_pop),
                       collapse = ", ")),
         sprintf("highest record Ho (pooled): %.2f", pf$max_ho))
    if (isTRUE(pf$looks_mac_filtered)) note(.report_text$diversity_looks_mac_filtered)
    note(.report_text$diversity_prior_filters)
  }

  if (!is.null(x$fis_by_call_rate)) {
    cat("\nFIS by call rate -- a check for null alleles and allele dropout\n")
    .print_table(.round_diversity_table(x$fis_by_call_rate, "fis_by_call_rate"))
    note(.report_text$diversity_fis_by_call_rate)
  }

  cat("\nWhat the Ho and He numbers are\n")
  sites <- st$sites
  ok_sites <- st$ok_sites
  mode_text <- if (st$complete_case) "complete-case" else sprintf("available-data, min_n = %d", st$min_n)
  if (any(ok_sites) && st$is_haplotype) {
    note(.report_text$diversity_sites_on_haplotypes)
    see("3. Which Stacks file for which statistic")
  } else if (any(ok_sites)) {
    retained <- n_used / n_rec
    ## Older result objects (before per-population counts) have no
    ## variant_records; they were scaled by every record.
    records <- if (is.null(st$variant_records)) stats::setNames(rep(n_rec, length(sites)), names(sites))
               else st$variant_records
    site_range <- range(sites[ok_sites])
    note("Scaled per population: He x (variant records with data in that population)",
         sprintf("/ (its sequenced sites) -- see `variant_records` and `sites_used` (%s sites).",
                 if (site_range[1] == site_range[2]) .big(site_range[1])
                 else paste0(.big(site_range[1]), "-", .big(site_range[2]))),
         sprintf("He estimated using %s of %s records (%.1f%% of them, %s).",
                 .big(n_used), .big(n_rec), 100 * retained, mode_text))
    if (any(!ok_sites)) {
      reasons <- vapply(names(sites)[!ok_sites], function(p) {
        if (sites[[p]] > 0)
          sprintf("%-22s sites = %s, less than its %s variant records",
                  p, .big(sites[[p]]), .big(records[[p]]))
        else sprintf("%-22s no sites value supplied", p)
      }, character(1))
      note("No autosomal value for:", reasons)
    }
    .print_table(.round_diversity_table(x$autosomal, "autosomal"))
    note(.report_text$diversity_autosomal_meaning)
    see("Appendix A. Stacks columns and RADdiversity")
    if (st$complete_case && retained < 0.5)
      note(sprintf("NOTE: only %.0f%% of variant records have complete data, so this rests on",
                   100 * retained),
           .report_text$diversity_complete_case_representative)
  } else {
    if (any(sites > 0))
      note("`sites` was supplied but every value is LESS than that population's variant records",
           .report_text$diversity_sites_implausible)
    else
      note(.report_text$diversity_per_ascertained_record)
    see("Heterozygosity per sequenced site")
  }

  cat("\n")
  rule("-")
  note(if (st$is_haplotype) .report_text$diversity_take_haplotype else .report_text$diversity_take_snp)
  note("", .report_text$diversity_take_common)
  see("3. Which Stacks file for which statistic")
  rule("-")

  cat("\nFor Weir & Cockerham FST, Jost's D and Weir & Goudet's beta between these\n")
  cat("populations, run differentiation_stats() on the same file.\n")

  hd <- x$he_difference
  if (!is.null(hd)) {
    cat(sprintf("\nHe difference (%s - %s): %+.4f  95%% CI [%+.4f, %+.4f]\n",
                hd$pop1, hd$pop2, hd$He_diff, hd$lo, hd$hi))
    if (st$boot == "loci") note(.report_text$diversity_he_difference_loci)
    else note(sprintf("This interval also resamples INDIVIDUALS (boot=\"%s\"), unlike a plain", st$boot),
              .report_text$diversity_he_difference_individuals)
  }
  if (length(st$files)) cat(sprintf("\nWrote %s\n", paste(st$files, collapse = ", ")))
  cat("\n")
  invisible(NULL)
}

## ---------------------------------------------------------------------------
## The `sites` argument
## ---------------------------------------------------------------------------

## Not exported. Turns diversity_stats()' `sites` argument into a numeric
## vector with one value per population, named and ordered as `pop_names`.
## Accepted shapes:
##   - NULL or 0: no per-site conversion (every population gets 0)
##   - one non-negative number: used for every population
##   - a named numeric vector, names = population names
##   - a data frame or matrix with a population column ("population"/"pop",
##     any case) and a sites column ("sites"/"site")
##   - a path to populations.sumstats_summary.tsv: read with
##     read_sumstats_summary(); uses the "All positions" block's `sites`, and
##     also returns that block's `variant_sites` as the attribute
##     "variant_sites" (see .autosomal_counts())
## Called twice by diversity_stats(): first with pop_names = NULL, before the
## VCF is read, only to fail fast on a malformed value; then with the real
## population names, to match values to populations.
.resolve_sites <- function(sites, pop_names = NULL) {
  variant_sites <- NULL
  if (is.null(sites)) {
    values <- 0
  } else if (is.character(sites)) {
    if (length(sites) != 1L) stop("sites: a file path must be a single string.", call. = FALSE)
    if (!file.exists(sites)) stop("sites: file not found: ", sites, call. = FALSE)
    all_positions <- read_sumstats_summary(sites)$all_positions
    values <- stats::setNames(all_positions$sites, all_positions$population)
    variant_sites <- stats::setNames(all_positions$variant_sites, all_positions$population)
  } else if (is.data.frame(sites) || is.matrix(sites)) {
    df <- as.data.frame(sites, stringsAsFactors = FALSE)
    column_names <- tolower(names(df))
    pop_col <- which(column_names %in% c("population", "pop"))
    site_col <- which(column_names %in% c("sites", "site"))
    if (length(pop_col) != 1L || length(site_col) != 1L)
      stop("sites: a data frame/matrix must have exactly one population column ",
           "(\"population\" or \"pop\") and exactly one sites column (\"sites\"). ",
           "Found columns: ", paste(names(df), collapse = ", "), call. = FALSE)
    values <- stats::setNames(suppressWarnings(as.numeric(df[[site_col]])),
                              as.character(df[[pop_col]]))
  } else if (is.numeric(sites)) {
    values <- sites
  } else {
    stop("sites must be NULL, a non-negative number, a named numeric vector or data ",
         "frame keyed by population, or a path to populations.sumstats_summary.tsv ",
         "(got class: ", paste(class(sites), collapse = "/"), ").", call. = FALSE)
  }

  if (anyNA(values) || any(values < 0))
    stop("sites: every value must be a non-negative, non-missing number.", call. = FALSE)
  one_for_all <- length(values) == 1L && is.null(names(values))
  if (!one_for_all && (is.null(names(values)) || any(!nzchar(names(values)))))
    stop("sites: more than one value must be named (or keyed) by population.", call. = FALSE)

  if (is.null(pop_names)) return(invisible(NULL))
  if (one_for_all) values <- stats::setNames(rep(as.numeric(values), length(pop_names)), pop_names)
  out <- .match_sites_to_pops(values, pop_names)
  if (!is.null(variant_sites)) attr(out, "variant_sites") <- variant_sites[pop_names]
  out
}

## Not exported. The two counts behind each population's per-site values,
## He_autosomal = He * variant_records / sites (see autosomal_het()).
##
## `sites` (Stacks' `Sites`, All-positions block) counts only the sites where
## that population has at least one genotyped individual (Stacks 2.68,
## SumStatsSummary::accumulate()). The matching numerator is therefore the
## number of VARIANT sites where that population has at least one genotyped
## individual -- not every record in the VCF. The two differ whenever a
## population is absent from some records: Stacks' `-r` filter blanks a
## population site by site while keeping the site for the others (with 3 or
## more populations and `-p` below their number), and every such record
## would otherwise inflate that population's He_autosomal.
##
## Where the count comes from:
##   * `sites` read from populations.sumstats_summary.tsv: that file's own
##     `Variant_Sites`, the exact count Stacks paired with `Sites`. It stays
##     right if records were removed in R afterwards (the mean He over the
##     remaining records still estimates the mean over all of them).
##   * `sites` given as numbers: the records of this VCF where the population
##     has at least one genotyped individual.
## Returns list(variant_records = named per-population counts, notes =
## warning texts for diversity_stats() to raise).
.autosomal_counts <- function(n_typed_all, sites_by_pop, n_records_read, is_haplotype) {
  present_in_vcf <- colSums(n_typed_all > 0)
  n_rec <- nrow(n_typed_all)
  from_file <- attr(sites_by_pop, "variant_sites")
  notes <- character(0)
  supplied <- any(sites_by_pop > 0) && !is_haplotype

  if (!is.null(from_file) && all(is.finite(from_file))) {
    variant_records <- stats::setNames(as.numeric(from_file), names(sites_by_pop))
    gap <- abs(present_in_vcf - variant_records) / pmax(variant_records, 1)
    if (supplied && any(gap > 0.005))
      notes <- c(notes, paste0(
        "Per-site values: the VCF has a different number of variant records per population ",
        "than populations.sumstats_summary.tsv's Variant_Sites (",
        paste(sprintf("%s: %s vs %s", names(gap)[gap > 0.005], present_in_vcf[gap > 0.005],
                      variant_records[gap > 0.005]), collapse = "; "),
        "). Ho_autosomal/He_autosomal use Variant_Sites. Check that both files come from the ",
        "same populations run; if records were removed in R, see ?diversity_stats (`sites`)."))
  } else {
    variant_records <- stats::setNames(as.numeric(present_in_vcf), names(sites_by_pop))
    if (supplied && !is.null(n_records_read) && n_rec < n_records_read)
      notes <- c(notes, paste0(
        "Per-site values: ", n_records_read - n_rec, " of ", n_records_read, " records were ",
        "removed after reading (filter_*()), but the removed variant sites are still in ",
        "`sites`, so Ho_autosomal/He_autosomal are too low. Pass `sites` as the path to ",
        "populations.sumstats_summary.tsv (which supplies the unfiltered count), or compute ",
        "per-site values on the unfiltered VCF."))
  }
  list(variant_records = variant_records, notes = notes)
}

## Not exported. Matches named sites values to `pop_names` in both
## directions: every population needs a value (error naming the missing
## ones), and a name that matches no population is reported, since it is
## most likely a typo.
.match_sites_to_pops <- function(values, pop_names) {
  missing_pops <- setdiff(pop_names, names(values))
  if (length(missing_pops))
    stop("sites: no value given for population(s): ", paste(missing_pops, collapse = ", "),
         ". Every population in the popmap needs a sites value.", call. = FALSE)
  extra <- setdiff(names(values), pop_names)
  if (length(extra))
    message("sites: ", length(extra), " name(s) not among this run's populations, ignored: ",
            paste(extra, collapse = ", "))
  values[pop_names]
}
