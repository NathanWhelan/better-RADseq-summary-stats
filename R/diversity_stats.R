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
#' privAr. Haplotype Ar and privAr depend on how many SNPs each RAD locus
#' holds, which changes with read length, enzyme and SNP calling. Compare them
#' only among populations analyzed together (the same run, filters and `g`):
#' which population is richer is reliable, but the size of the difference
#' also depends on SNPs per locus. The printed result says which numbers to
#' take from each file. See `vignette("rationale")`, section 3.
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
#'   `se_individuals = TRUE` instead (see "How standard errors are
#'   calculated" below).
#' @param sites Number of sequenced sites, to turn Ho and He into per-site
#'   values (`Ho_autosomal`, `He_autosomal`; SNP VCF only). Default `NULL`: no
#'   per-site values. Give one of:
#'   * the path to Stacks' `populations.sumstats_summary.tsv`: each
#'     population's own `Sites` (from the "All positions (variant and fixed)"
#'     block) is read with [read_sumstats_summary()];
#'   * one number, used for every population;
#'   * a named numeric vector, one value per population (names as in the
#'     popmap);
#'   * a data frame with a `population` (or `pop`) column and a `sites` column.
#'
#'   **How per-site values are calculated.** For each population, He summed
#'   over the SNP records in the data, divided by `sites`. (In practice: the
#'   mean He over the records the population uses, times the number of records
#'   where it has at least one genotyped individual, shown as
#'   `variant_records`.) Fixed sites add nothing to the sum but count in
#'   `sites`. Ho is done the same way.
#'
#'   **What `sites` must count.** Every sequenced site, variant or fixed, of
#'   the RAD loci in the data, counting a site for a population only where at
#'   least one of its individuals is genotyped there. This is how Stacks
#'   counts `Sites`. Use the same popmap that was given to Stacks.
#'
#'   **When Stacks' `Sites` is right.** Only for the data exactly as Stacks'
#'   `populations` wrote them. Filtering afterwards (in R or with another
#'   program) changes things:
#'   * SNPs removed, for example by [filter_mac()]: their He drops out of the
#'     sum, so per-site values go down. For MAC and MAF filters this is real
#'     rare-variant diversity being lost, as with Stacks' own `--min-mac`;
#'     report the threshold (`vignette("rationale")`, section 1). A RAD locus
#'     whose every SNP was removed this way was still sequenced, so its sites
#'     stay in `sites`.
#'   * Whole RAD loci removed as unreliable (low call rate, paralogs): their
#'     sites should leave `sites` too, but Stacks' `Sites` still counts them,
#'     so per-site values come out too low. Give `sites` for the loci you kept,
#'     or use [pi_allsites()] on an all-sites VCF.
#'   * A VCF thinned to one SNP per RAD locus (Stacks' `--write-single-snp`,
#'     [filter_thin_one_snp()]) no longer holds all the SNPs. Compute per-site
#'     values from the unthinned SNP VCF.
#'
#'   A warning is given when records or whole loci were removed after reading
#'   (the `filter_*()` functions record which filter removed them, so loci
#'   emptied by [filter_maf()] or [filter_mac()] give no loci warning), or when
#'   the VCF's SNP count does not match `Variant_Sites` in the Stacks file.
#' @param min_n Minimum genotyped individuals a population needs at a record
#'   to use that record. Default `2`, the smallest number for which He and
#'   FIS are defined. Ignored when `complete_case = TRUE`.
#' @param complete_case Use a record only if every individual of every
#'   population is genotyped there (Schmidt et al. 2021, recommendation (b)).
#'   Default `FALSE`.
#' @param se_individuals Also compute the uncertainty from which individuals
#'   were sampled, and add three columns each for Ho, He and FIS: `_se_ind`
#'   (a jackknife over individuals) and `_se_combined` (locus and individual
#'   uncertainty together). **Report `_se_combined` for Ho, He and FIS.** Ar and
#'   privAr keep their locus-based `_se`, `_lo` and `_hi`, which already cover
#'   both sources. See "How standard errors are calculated" below. Default
#'   `FALSE`, because the jackknife over individuals takes longer on large
#'   datasets.
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
#'       population, each with its locus jackknife SE (`_se`) and 95% locus
#'       bootstrap interval (`_lo`, `_hi`). With `se_individuals = TRUE`, also
#'       `_se_ind` and `_se_combined` for Ho, He and FIS.}
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
#' @section How standard errors are calculated:
#' A standard error (SE) says how much a number would change if the study were
#' repeated; the estimate plus or minus 1.96 SE is a 95% confidence interval.
#' A population's He (or Ho, FIS, Ar) is uncertain for two reasons:
#' * you typed **some of the RAD loci** in the genome, and other loci would
#'   give a slightly different value;
#' * you caught **some of the individuals** in the population, and other
#'   individuals would give a slightly different value.
#'
#' The columns measure these as follows:
#' * `_se`: a jackknife over RAD loci. Each RAD locus is left out in turn and
#'   the statistic recomputed; how much those values spread gives the SE. It
#'   holds the individuals fixed.
#' * `_lo`, `_hi`: a bootstrap over RAD loci. Loci are drawn at random, with
#'   replacement, thousands of times; the middle 95% of the results is the
#'   interval. It also holds the individuals fixed.
#' * `_se_ind` (with `se_individuals = TRUE`): a jackknife over individuals.
#'   Each individual is left out in turn. It holds the loci fixed.
#' * `_se_combined`: both sources together, `sqrt(_se^2 + _se_ind^2)`.
#'
#' Whole RAD loci, not single SNPs, are resampled, because SNPs on one RAD
#' locus are inherited together.
#'
#' **Which to report.** For **Ho, He and FIS**, report `_se_combined`. For
#' **Ar and privAr**, report `_se` (or `_lo`, `_hi`). This rule comes from
#' simulated populations with a known true value, where each simulated study
#' was repeated 300 times with new loci and new individuals, and we counted how
#' often the 95% interval contained the truth (it should be about 95%). With
#' the locus SE alone, Ho and FIS intervals contained the truth only 30-65% of
#' the time in most settings once individuals differed in inbreeding: an inbred individual is
#' homozygous at many loci at once, so loci are not independent. With the
#' individual SE alone, He intervals contained it 47-86% of the time, because
#' He depends mostly on which loci were typed. The combined SE contained it
#' 90-99% of the time for all three, often a little more than 95% because some
#' noise is counted in both SEs. For Ar and privAr the locus SE alone did well
#' (88-96%). A bootstrap over individuals is not offered for reporting:
#' drawing an individual twice makes the sample look less diverse than it is,
#' and its He intervals contained the truth under 10% of the time. Details are
#' in `vignette("rationale")`, section 4, and `inst/sims/uncertainty_sources.R`.
#'
#' **What no SE covers.** Relatives in the sample (they bias the estimates
#' themselves; screen with [kinship_check()]), filtering choices and allele
#' dropout, and linkage between separate RAD loci.
#'
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
#' Nei, M. & Roychoudhury, A.K. (1974) Sampling variances of heterozygosity
#' and genetic distance. *Genetics* 76:379-390. (Loci and individuals as two
#' sources of uncertainty.)
#'
#' Shao, J. & Wu, C.F.J. (1989) A general theory for jackknife variance
#' estimation. *Annals of Statistics* 17:1176-1197.
#' \doi{10.1214/aos/1176347263}
#'
#' Owen, A.B. (2007) The pigeonhole bootstrap. *Annals of Applied Statistics*
#' 1:386-411. \doi{10.1214/07-AOAS122}
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
  sites_by_pop <- .resolve_sites(sites, pop_names, verbose)

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
    jack <- .jack_individuals(H, pops, rows, counts, n_typed, n_het, min_n_used)
    se_ind <- jack$se
    ## Built-in check that the two code paths agree.
    gap <- suppressWarnings(max(abs(jack$full - point[names(jack$full)]), na.rm = TRUE))
    if (is.finite(gap) && gap > 1e-8)
      warning("Individual jackknife: its full-sample estimates differ from the point ",
              "estimates by ", signif(gap, 3), ". Please report this.", call. = FALSE)
    ## Records with exactly 2 genotyped individuals drop out of a replicate
    ## when either is left out, which inflates these SEs (.jackknife_boundary()).
    jackknife_boundary <- .jackknife_boundary(n_typed, min_n_used)
    inflated <- pop_names[is.finite(jackknife_boundary) &
                            jackknife_boundary > .jackknife_boundary_share]
    if (length(inflated))
      warning("The individual SEs (_se_ind, and so _se_combined), mainly He's, are ",
              "probably too large for ",
              paste(sprintf("%s (%.0f%% of its records have only 2 genotyped individuals)",
                            inflated, 100 * jackknife_boundary[inflated]), collapse = ", "),
              ". Leaving one of those 2 out leaves too few to compute He, so the record drops ",
              "out of that jackknife replicate. A higher min_n (e.g. min_n = 3) avoids it; see ",
              "?diversity_stats, \"How standard errors are calculated\".", call. = FALSE)
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
  counts_autosomal <- .autosomal_counts(n_typed_all, sites_by_pop, H, is_haplotype)
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
                   prior_filters = prior_filters, filter_log = H$filter_log,
                   seed = seed, files = written)
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
      if (isTRUE(st$se_individuals)) "; _se_ind: jackknife over individuals; _se_combined: both (report it for Ho, He, Fis)",
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
    cat("  Ar and privAr: compare only populations analyzed together; which is richer is reliable, how much richer depends on SNPs per locus.\n")
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
  if (st$is_haplotype && any(st$sites > 0))
    out <- c(out, "`sites` was ignored: per-site values need the SNP VCF.")
  else if (!any(st$ok_sites) && any(st$sites > 0))
    out <- c(out, "`sites` was ignored: every value is smaller than the number of variant records.")
  if (isTRUE(st$se_individuals)) {
    share <- st$jackknife_boundary
    inflated <- names(share)[is.finite(share) & share > .jackknife_boundary_share]
    if (length(inflated))
      out <- c(out, sprintf("_se_ind and _se_combined (mainly He's) are probably too large for %s (many records with only 2 genotyped individuals); a higher min_n avoids it.",
                            paste(inflated, collapse = ", ")))
  } else {
    out <- c(out, "_se holds the individuals fixed. For Ho, He and Fis, report _se_combined: re-run with se_individuals = TRUE (see ?diversity_stats).")
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
  note(if (isTRUE(st$se_individuals)) .report_text$diversity_se_ind_columns
       else .report_text$diversity_se_loci_only)
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
  steps <- .format_filter_log(st$filter_log)
  if (length(steps)) {
    cat("\nFilters applied in R after reading the VCF\n")
    note(steps)
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
  if (st$is_haplotype && any(sites > 0)) {
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
.resolve_sites <- function(sites, pop_names = NULL, verbose = TRUE) {
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
  ## A repeated population would silently use its first value only.
  if (anyDuplicated(names(values)))
    stop("sites: population(s) given more than one value: ",
         paste(unique(names(values)[duplicated(names(values))]), collapse = ", "),
         ". Give each population exactly one sites value.", call. = FALSE)

  if (is.null(pop_names)) return(invisible(NULL))
  if (one_for_all) values <- stats::setNames(rep(as.numeric(values), length(pop_names)), pop_names)
  out <- .match_sites_to_pops(values, pop_names, verbose)
  if (!is.null(variant_sites)) attr(out, "variant_sites") <- variant_sites[pop_names]
  out
}

## Not exported. The SNP count behind each population's per-site values,
## He_autosomal = He * variant_records / sites (see autosomal_het()), and the
## warnings about it.
##
## WHAT THE FORMULA DOES. He_autosomal should be the population's He summed
## over its SNPs, divided by its sequenced sites. He (the mean over the records
## the population uses) times `variant_records` rebuilds that sum.
## `variant_records` is therefore the number of records in THIS data where the
## population has at least one genotyped individual. That matches how Stacks
## counts `Sites` (a site counts for a population only where it has a
## genotyped individual; SumStatsSummary::accumulate(), Stacks 2.68), and it
## handles Stacks' `-r`, which blanks a population at some sites and keeps
## the site for the others.
##
## WHY NOT STACKS' Variant_Sites. It also counts SNPs removed after Stacks.
## Multiplying by it treats every removed SNP as if it had the average He of
## the ones kept. After a MAC or MAF filter the removed SNPs had much lower He,
## so that inflated per-site He (+41% after filter_mac(3) in a simulation with
## 20 diploids). Using the data's own count, removed SNPs simply add nothing,
## which is the same loss Stacks' own --min-mac causes.
##
## The notes (warnings) cover the cases where `sites` may no longer match the
## data: records removed after reading, whole RAD loci removed after reading
## (told apart by H$filter_log: loci emptied by a MAF or MAC filter keep their
## sites and need no warning), and a VCF whose SNP count differs from the
## Stacks file's Variant_Sites.
##   n_typed_all   records x populations, genotyped individuals (all records)
##   sites_by_pop  from .resolve_sites(); attribute "variant_sites" when read
##                 from populations.sumstats_summary.tsv
##   H             the data (for n_records_read, n_loci_read, filter_log and
##                 locus_raw)
## Returns list(variant_records = named per-population counts, notes =
## warning texts for diversity_stats() to raise).
.autosomal_counts <- function(n_typed_all, sites_by_pop, H, is_haplotype) {
  variant_records <- stats::setNames(as.numeric(colSums(n_typed_all > 0)), names(sites_by_pop))
  notes <- character(0)
  if (!any(sites_by_pop > 0) || is_haplotype)
    return(list(variant_records = variant_records, notes = notes))

  n_rec <- nrow(n_typed_all)
  n_loci <- length(unique(H$locus_raw))
  records_removed <- if (is.null(H$n_records_read)) 0 else H$n_records_read - n_rec
  loci_removed <- if (is.null(H$n_loci_read)) 0 else H$n_loci_read - n_loci
  ## Which filters removed them (H$filter_log; see R/filter_loci.R). A RAD
  ## locus emptied by a minor allele count or frequency filter was still
  ## sequenced: its sites are now fixed sites, which is also how Stacks' own
  ## --min-mac treats them, so they rightly stay in `sites`. A locus removed by
  ## any other filter (call rate, excess heterozygosity) was removed as
  ## unreliable, so its sites should leave `sites`. Loci removed without a
  ## filter_*() function (or before the log existed) have no known reason.
  log <- H$filter_log
  if (is.null(log)) log <- .empty_filter_log()
  by_allele_frequency <- log$filter %in% c("filter_maf", "filter_mac")
  loci_quality <- sum(log$loci_removed[!by_allele_frequency])
  loci_unexplained <- max(0, loci_removed - sum(log$loci_removed))
  filters_named <- function(rows) paste(unique(log$filter[rows]), collapse = ", ")
  if (records_removed > 0) {
    removed_by <- log$records_removed > 0
    notes <- c(notes, paste0(
      "Per-site values: ", .big(records_removed), " of ", .big(H$n_records_read),
      " records were removed after reading",
      if (any(removed_by)) paste0(" (by ", filters_named(removed_by), ")") else "",
      ", so their He adds nothing to Ho_autosomal and He_autosomal. After a MAC or MAF filter ",
      "this is real rare-variant diversity being lost; report the threshold. See ",
      "?diversity_stats, `sites`.",
      if (any(log$filter == "filter_thin_one_snp" & removed_by))
        paste0(" filter_thin_one_snp() kept one SNP per RAD locus, so these per-site values ",
               "are far too low: compute them from the unthinned SNP data.")
      else ""))
  }
  if (loci_quality > 0) {
    quality_rows <- !by_allele_frequency & log$loci_removed > 0
    notes <- c(notes, paste0(
      "Per-site values: ", .big(loci_quality), " whole RAD loci were removed after reading by ",
      filters_named(quality_rows), ". Their sequenced sites should no longer count, but ",
      "Stacks' Sites still counts them, so if `sites` came from ",
      "populations.sumstats_summary.tsv, Ho_autosomal and He_autosomal are too low. Give `sites` ",
      "for the loci you kept, or use pi_allsites() on an all-sites VCF."))
  }
  if (loci_unexplained > 0)
    notes <- c(notes, paste0(
      "Per-site values: ", .big(loci_unexplained), " whole RAD loci were removed after reading, ",
      "but not by a filter_*() function, so the reason is not known. If they were removed as ",
      "unreliable (low call rate, paralogs), their sites should no longer count and Stacks' ",
      "Sites is too large, so Ho_autosomal and He_autosomal are too low: give `sites` for the ",
      "loci you kept, or use pi_allsites() on an all-sites VCF. If they only lost their SNPs to ",
      "a minor allele count or frequency filter, their sites were still sequenced and Stacks' ",
      "Sites is right."))

  from_file <- attr(sites_by_pop, "variant_sites")
  if (records_removed == 0 && !is.null(from_file) && all(is.finite(from_file))) {
    gap <- abs(variant_records - from_file) / pmax(from_file, 1)
    if (any(gap > 0.005)) {
      thinned <- n_rec / n_loci < 1.05 && any(from_file > 1.2 * variant_records)
      notes <- c(notes, paste0(
        "Per-site values: this VCF's SNP count does not match Variant_Sites in ",
        "populations.sumstats_summary.tsv (",
        paste(sprintf("%s: %s in the VCF, %s in the file", names(gap)[gap > 0.005],
                      .big(variant_records[gap > 0.005]), .big(from_file[gap > 0.005])),
              collapse = "; "),
        "). Either the VCF was filtered or thinned before it was read, or the two files come ",
        "from different populations runs. Ho_autosomal and He_autosomal use the VCF's own SNPs, so ",
        "SNPs missing from it add nothing, and if whole loci are missing, `sites` is too large; ",
        "both make the per-site values too low.",
        if (thinned) paste0(" The VCF has about one SNP per RAD locus, so it looks thinned ",
                            "(e.g. --write-single-snp): use the unthinned SNP VCF for per-site values.")
        else ""))
    }
  }
  list(variant_records = variant_records, notes = notes)
}

## Not exported. Matches named sites values to `pop_names` in both
## directions: every population needs a value (error naming the missing
## ones), and a name that matches no population is reported, since it is
## most likely a typo.
.match_sites_to_pops <- function(values, pop_names, verbose = TRUE) {
  missing_pops <- setdiff(pop_names, names(values))
  if (length(missing_pops))
    stop("sites: no value given for population(s): ", paste(missing_pops, collapse = ", "),
         ". Every population in the popmap needs a sites value.", call. = FALSE)
  extra <- setdiff(names(values), pop_names)
  if (length(extra))
    .inform(verbose, "sites: ", length(extra), " name(s) not among this run's populations, ignored: ",
            paste(extra, collapse = ", "))
  values[pop_names]
}
