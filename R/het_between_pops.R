###############################################################################
#
#  R/het_between_pops.R -- does mean heterozygosity differ between populations?
#  Each individual contributes ONE number, and the test counts both which
#  individuals and which loci were sampled.
#
#  WHY NOT THE USUAL TEST. Bootstrapping over loci, or a paired Wilcoxon on
#  per-locus values, treats LOCI as the only replicate. Loci are repeated
#  measures on the same individuals: when individuals differ in inbreeding,
#  that uncertainty does not shrink as loci are added, and a locus-only
#  interval becomes far too narrow (Van Dongen 1995).
#
#  WHY NOT WELCH'S t ALONE. The opposite mistake: Welch's t on one value per
#  individual treats the LOCI as fixed. Two populations with the same
#  genome-wide heterozygosity still differ a little at any given set of loci,
#  because drift acts on each locus separately.
#
#  So the primary test combines the two (.locus_variance_of_difference(),
#  .two_sample()). In het_between_pops_selftest()'s simulations of two
#  populations with identical true heterozygosity and their own allele
#  frequencies, a locus bootstrap and Welch's t each went wrong in the
#  situation the other handles, and the combined test held its nominal rate
#  in both (inst/sims/het_test_null.R has the full grid). Full discussion in
#  vignette("rationale"), section 5.
#
#  WHAT TO READ IN THE REPORT, IN ORDER
#    1. the missingness confound check -- if call rate correlates with
#       heterozygosity, the test measures library quality, not biology
#    2. g2 -- whether individuals differ in inbreeding
#    3. the tests themselves: the combined test, with Welch's t and Wilcoxon
#       (individuals only) for comparison
#
#  Needs the VCF: per-individual heterozygosity cannot be recovered from
#  populations.sumstats.tsv, which holds only per-population counts.
#
#  STEPS (functions below the exported one):
#    .warn_failed_individuals() names possible failed libraries (keeps them)
#    .check_population_loci() stops if a population has no usable locus
#    .individual_het()        pooled locus filter, per-individual table
#    .population_locus_sets() which loci each population uses (R/vcf_io.R)
#    .population_summary()    per-population mean heterozygosity and F
#    .missingness_confound()  does call rate track heterozygosity?
#    .g2_by_population()      identity disequilibrium
#    .pairwise_tests()        combined test, Welch, Wilcoxon, Hedges' g per pair
#    .omnibus_tests()         Welch's ANOVA and Kruskal-Wallis, 3+ populations
#  Checked by tests/testthat/test-diversity-stats.R, test-regression.R
#  (golden values), test-individual-inbreeding.R and
#  het_between_pops_selftest().
#
###############################################################################

#' Test whether mean heterozygosity differs between populations
#'
#' Gives each individual one value and compares populations with a test that
#' counts both sources of uncertainty: which individuals were caught, and which
#' loci were typed. It also adds a missingness-confound check and identity
#' disequilibrium (g2) per population, which shows whether individuals differ
#' in inbreeding. Printing the result shows the main tables. `summary()` is
#' short: the two tests, the populations with their g2, what to report and a
#' list of checks. `summary(details = TRUE)` is the full report, with every test
#' column and the notes on interpretation. `outdir` also writes the tables to
#' TSV files.
#'
#' @details
#' **The test.** `diff` is the difference in mean individual value. Its
#' standard error, `se_combined`, adds two parts:
#' * how much individuals vary (the variance Welch's t uses), and
#' * how much the difference would change with other loci: a jackknife over
#'   RAD loci, without counting genotype noise a second time.
#'
#' `ci_lo`, `ci_hi` and `p_combined` come from it (a t distribution, with
#' Welch-Satterthwaite degrees of freedom `df`). Report these. Why both parts:
#' loci alone go wrong when individuals differ in inbreeding, because an
#' inbred individual is homozygous at many loci at once. Individuals alone
#' (Welch's t, `p_welch`) go wrong when the populations are differentiated,
#' because each population's own allele frequencies make its heterozygosity at
#' the typed loci differ a little from its genome-wide value. In simulations
#' with identical true heterozygosity, Welch's t alone rejected up to 14% of
#' the time at FST 0.05 and up to 30% at FST 0.2, at a nominal 5%, when
#' individuals were alike in inbreeding; the combined test stayed at 3-7.5% in
#' every setting
#' (`inst/sims/het_test_null.R`). `p_welch` and the Wilcoxon rank-sum
#' `p_wilcox` (individuals only) are kept for comparison.
#'
#' **Three or more populations.** Start with `omnibus`, an overall test of
#' whether any population differs (Welch's one-way ANOVA, with Kruskal-Wallis
#' as the distribution-free check). Interpret the pairwise tests only if the
#' overall test is significant; they then say which pairs differ. The overall
#' tests use individuals only, like Welch's t, so with differentiated
#' populations whose individuals are alike in inbreeding they reject too
#' readily; a pair is then confirmed by its `p_combined_BH`.
#'
#' Two questions, same machinery. `pairwise_tests` compares mean individual
#' heterozygosity, which is observed heterozygosity, He x (1 - F): do the
#' populations differ in heterozygosity? That reflects diversity and
#' inbreeding together; to compare He itself, use [diversity_stats()].
#' `pairwise_F_tests` compares the individual inbreeding coefficient `F`
#' ([individual_inbreeding()], expected heterozygosity from each individual's
#' own population): do they differ in inbreeding? Report the one that matches
#' your question.
#'
#' **Which loci are used.** An individual's heterozygosity in
#' `individual_heterozygosity` and the missingness check use the loci
#' genotyped in at least `min_call` of ALL individuals pooled. Everything that
#' describes or compares populations uses loci that each population clears
#' on its own: `population_summary` uses each population's own loci, and a
#' pairwise test uses the loci BOTH populations of that pair clear. A pooled
#' rule would let a well-covered population carry a poorly covered one. An
#' individual not genotyped at any of the loci a table or test uses has no
#' value there and is left out of it (and only it); `n` counts the individuals
#' with a value. A population with no locus at `min_call` stops the function.
#'
#' **Failed individuals.** Nothing is removed. A failed library should be
#' removed during assembly or filtering, and one still in the popmap was
#' probably kept for a reason. Instead, an individual genotyped at fewer than
#' `min_loci` records (when the rest of its population has that many), or at
#' under half of the records where the rest of its population is genotyped (the dDocent tutorial's cutoff; judged within its
#' own population, as in Cerca et al. 2021), is named in a warning and in
#' `settings$flagged_individuals`. Keeping it means: its heterozygosity rests
#' on few loci yet counts as a full individual in the means and tests; heavy
#' missing data usually comes with allele dropout, which biases its
#' heterozygosity low and its F high; and it lowers its population's call
#' rate, so fewer loci clear `min_call`. Remove it from the popmap (or with
#' [filter_samples()]) if it was not meant to be analyzed.
#'
#' @inheritParams diversity_stats
#' @param vcf Path to a Stacks VCF (`populations.snps.vcf` or
#'   `populations.haps.vcf`, optionally gzip-compressed), or the object
#'   returned by [read_stacks_vcf()], optionally passed through `filter_*()`
#'   functions.
#' @param min_call Minimum fraction of individuals genotyped at a locus for the
#'   locus to be used, applied within each population (see Details). Default
#'   `0.9`.
#' @param min_loci Minimum number of loci. An individual genotyped at fewer
#'   records is named as a possible failed library (see Details; it is kept),
#'   a pair of populations sharing fewer loci is reported as `NA`, as is the
#'   overall test when fewer loci clear `min_call` in every population, and a
#'   message says so if fewer pass the pooled filter. Default `50`.
#' @param nboot_g2 Bootstrap replicates over individuals for the g2
#'   confidence interval. Default `1000`, as in [identity_disequilibrium()].
#' @param seed Random seed for the g2 bootstrap. Default `NULL`: use R's
#'   current random-number stream (call `set.seed()` first for a reproducible
#'   interval). A number makes the result reproducible on its own and leaves
#'   your session's random-number stream as it was.
#' @return An object of class `raddiv_het`, a list of:
#'   \describe{
#'     \item{individual_heterozygosity}{One row per individual: loci called
#'       and heterozygosity on the pooled locus set, and individual `F` on its
#'       population's own loci. `NA` where the individual is genotyped at
#'       none of those loci.}
#'     \item{population_summary}{Per population: individuals, loci, mean, SD,
#'       SE, range of individual heterozygosity, and mean `F`.}
#'     \item{omnibus}{With 3 or more populations: one overall test of whether
#'       any population differs, for heterozygosity and for `F` -- Welch's
#'       one-way ANOVA (`welch_F`, `welch_df1`, `welch_df2`, `p_welch`) and
#'       Kruskal-Wallis (`kruskal_chisq`, `kruskal_df`, `p_kruskal`) -- on the
#'       loci that clear `min_call` in every population (`n_loci`). `NULL`
#'       for 2 populations, where the pairwise test is the only test.}
#'     \item{pairwise_tests}{Heterozygosity, one row per pair of populations:
#'       means, the difference `diff` with `se_combined`, its 95% interval
#'       (`ci_lo`, `ci_hi`), `df` and `p_combined` from the combined test (see
#'       Details); `p_welch` and `p_wilcox` for comparison; Hedges' g; and
#'       Benjamini-Hochberg adjusted p-values (`_BH`).}
#'     \item{pairwise_F_tests}{The same tests on `F`.}
#'     \item{g2}{Identity disequilibrium per population, with bootstrap CI.}
#'     \item{missingness_confound}{Correlation between an individual's call
#'       rate and its heterozygosity, overall (`overall`) and within each
#'       population (`within`).}
#'     \item{settings}{The settings of this run.}
#'   }
#'   Values are stored at full precision; `print()`, `summary()` and the TSV
#'   files round them. With `outdir`, `individual_heterozygosity`,
#'   `pairwise_tests` and `pairwise_F_tests` are written to
#'   `individual_heterozygosity.<stem>.tsv`,
#'   `het_between_pops_tests.<stem>.tsv` and
#'   `het_between_pops_F_tests.<stem>.tsv`, and with 3 or more populations
#'   `omnibus` to `het_between_pops_omnibus.<stem>.tsv`.
#' @references
#' Cerca, J., Maurstad, M.F., Rochette, N.C., et al. (2021) Removing the bad
#' apples: a simple bioinformatic method to improve loci-recovery in de novo
#' RADseq data for non-model organisms. *Methods in Ecology and Evolution*
#' 12:805-817. (Failed individuals.)
#'
#' Van Dongen, S. (1995) How should we bootstrap allozyme data? *Heredity*
#' 74:445-447.
#'
#' Welch, B.L. (1947) The generalization of "Student's" problem when several
#' different population variances are involved. *Biometrika* 34:28-35.
#'
#' Welch, B.L. (1951) On the comparison of several mean values: an alternative
#' approach. *Biometrika* 38:330-336.
#'
#' Kruskal, W.H. & Wallis, W.A. (1952) Use of ranks in one-criterion variance
#' analysis. *Journal of the American Statistical Association* 47:583-621.
#'
#' Benjamini, Y. & Hochberg, Y. (1995) Controlling the false discovery rate: a
#' practical and powerful approach to multiple testing. *Journal of the Royal
#' Statistical Society B* 57:289-300.
#'
#' Hedges, L.V. (1981) Distribution theory for Glass's estimator of effect size
#' and related estimators. *Journal of Educational Statistics* 6:107-128.
#'
#' David, P., Pujol, B., Viard, F., Castella, V. & Goudet, J. (2007) Reliable
#' selfing rate estimates from imperfect population genetic data. *Molecular
#' Ecology* 16:2474-2487. (g2.)
#' @section Tables for a paper:
#' Put in the table exactly what the short `summary()` reports.
#' `summary(x)$tables` holds those tables as data frames: `heterozygosity` and
#' `inbreeding` (one row per pair of populations, every pair: the difference as
#' `difference (SE)`, its 95% interval, `p (combined)` and Hedges' g; the p-value
#' is Benjamini-Hochberg adjusted when there is more than one pair) and
#' `populations` (n, loci, mean individual heterozygosity, its spread among
#' individuals, mean F, and g2 with its 95% interval). Report the combined test.
#' The Welch and Wilcoxon p-values are for comparison, and a population's mean
#' heterozygosity has no SE here, because populations are compared with the
#' test, never by overlapping intervals; they are in `pairwise_tests`, for a
#' supplement.
#' @examples
#' # A toy dataset shipped with the package (4 and 3 individuals -- far too
#' # few for a real test, which needs individuals, not loci).
#' vcf    <- system.file("extdata", "small.haps.vcf", package = "RADdiversity")
#' popmap <- system.file("extdata", "small_popmap.tsv", package = "RADdiversity")
#' set.seed(1)
#' res <- het_between_pops(vcf, popmap, min_call = 0.5)
#' res                       # the main tables
#' res$individual_heterozygosity
#' summary(res)              # short: the two tests, the populations, the checks
#' summary(res, details = TRUE)   # the full report, with every column explained
#' summary(res)$tables$heterozygosity   # the table to put in a paper
#' @export
het_between_pops <- function(vcf, popmap, min_call = 0.9, min_loci = 50L, nboot_g2 = 1000L,
                             outdir = NULL, stem = NULL, seed = NULL, verbose = TRUE) {
  ## ---- 1. Check the arguments ----------------------------------------------
  if (!(is.numeric(min_call) && length(min_call) == 1L && !is.na(min_call) &&
        min_call >= 0 && min_call <= 1))
    stop("min_call must be in 0-1 (a single fraction of individuals, e.g. 0.9).", call. = FALSE)
  min_loci <- .check_count(min_loci, "min_loci", min = 1)
  nboot_g2 <- .check_count(nboot_g2, "nboot_g2")
  .check_flag(verbose, "verbose")
  .check_seed(seed)
  .check_run_inputs(vcf, popmap, stem, outdir)
  if (!is.null(seed)) {
    restore_rng <- .save_rng_state()
    on.exit(restore_rng(), add = TRUE)
    set.seed(seed)
  }

  ## ---- 2. Read the data -----------------------------------------------------
  H <- .resolve_H(vcf, verbose = verbose)
  pops <- .resolve_pops(popmap, H$samples, verbose = verbose)
  if (length(pops) < 2) stop("Need at least 2 populations.", call. = FALSE)
  .stop_tiny_pops(pops, "A test between populations needs at least 2 individuals in each.")

  ## ---- 3. Individuals: failed-library check, pooled loci --------------------
  ## Every individual in the popmap is kept; one that looks like a failed
  ## library is named in a warning (see "Failed individuals" in the help).
  flagged <- .warn_failed_individuals(H, pops, min_loci)
  individuals <- .individual_het(H, pops, min_call, min_loci, verbose)

  ## ---- 4. Populations -------------------------------------------------------
  locus_sets <- .population_locus_sets(H, pops, min_call)
  .check_population_loci(H, pops, locus_sets, min_call)
  per_pop <- .population_summary(H, pops, locus_sets)
  table_ind <- individuals$table
  table_ind$F <- unname(per_pop$F_by_individual[table_ind$sample])

  confound <- .missingness_confound(individuals$call_rate, individuals$het,
                                    table_ind$population, names(pops))
  g2 <- .g2_by_population(H, pops, locus_sets, nboot_g2)

  ## ---- 5. Tests: overall (3+ populations), then pairwise --------------------
  omnibus <- .omnibus_tests(H, pops, locus_sets, min_call, min_loci, verbose)
  tests <- .pairwise_tests(H, pops, locus_sets, min_call, min_loci, verbose)

  ## ---- 6. Files (only with outdir) ------------------------------------------
  written <- character(0)
  if (!is.null(outdir)) {
    stem <- .derive_stem(vcf, stem)
    written <- .write_tables(
      list(individual = .round_het_table(table_ind, "individual_heterozygosity"),
           omnibus = if (!is.null(omnibus)) .round_het_table(omnibus, "omnibus"),
           tests = .round_het_table(tests$heterozygosity, "pairwise_tests"),
           F_tests = .round_het_table(tests$F, "pairwise_tests")),
      outdir,
      c(individual = sprintf("individual_heterozygosity.%s.tsv", stem),
        omnibus = sprintf("het_between_pops_omnibus.%s.tsv", stem),
        tests = sprintf("het_between_pops_tests.%s.tsv", stem),
        F_tests = sprintf("het_between_pops_F_tests.%s.tsv", stem)))
  }

  settings <- list(n_loci_pooled = individuals$n_loci_pooled,
                   n_individuals = nrow(table_ind), n_pops = length(pops),
                   min_call = min_call, min_loci = min_loci, nboot_g2 = nboot_g2,
                   flagged_individuals = flagged, seed = seed, files = written)
  structure(list(individual_heterozygosity = table_ind,
                 population_summary = per_pop$table,
                 omnibus = omnibus,
                 pairwise_tests = tests$heterozygosity,
                 pairwise_F_tests = tests$F,
                 g2 = g2,
                 missingness_confound = confound,
                 settings = settings),
            class = "raddiv_het")
}

## ---------------------------------------------------------------------------
## Steps
## ---------------------------------------------------------------------------

## Not exported. FAILED INDIVIDUALS: a quick check, never a filter. A library
## that failed (too little DNA, a bad barcode, contamination) should have been
## removed during assembly or filtering, and one that is still in the popmap
## was probably kept on purpose, so nothing is removed here. Instead each
## individual that looks failed is named in a warning that says what keeping
## it does. An individual is flagged when it is
##   * genotyped at fewer than `min_loci` records, although the rest of its
##     population is genotyped at `min_loci` or more (a small dataset is not
##     a failed library), or
##   * missing at more than half (`.failed_call_rate`) of the records where
##     at least one OTHER individual of its own population is genotyped.
## The second rule uses the dDocent tutorial's 50% cutoff for individuals, and
## judges each individual against its own population's records, as Cerca et
## al. (2021) do: a record that a whole population lacks (Stacks' -r, a
## restriction-site mutation in a divergent population) says nothing about an
## individual's library.
## Returns the flagged sample names (character(0) when none), invisibly.
.warn_failed_individuals <- function(H, pops, min_loci = 50L) {
  per_pop <- lapply(names(pops), function(p) {
    typed <- !is.na(H$A1[, pops[[p]], drop = FALSE])       # records x individuals
    others_typed <- rowSums(typed) - typed > 0              # someone else typed there
    with_data <- colSums(others_typed)
    data.frame(sample = pops[[p]], population = p, called = colSums(typed),
               with_data = with_data,
               call_rate = ifelse(with_data > 0, colSums(typed & others_typed) / with_data,
                                  NA_real_),
               row.names = NULL)
  })
  tab <- do.call(rbind, per_pop)
  failed <- tab[(tab$called < min_loci & tab$with_data >= min_loci) |
                  (!is.na(tab$call_rate) & tab$call_rate < .failed_call_rate), , drop = FALSE]
  if (nrow(failed)) {
    shown <- utils::head(seq_len(nrow(failed)), 10)
    warning(nrow(failed), " individual(s) look like failed libraries (genotyped at fewer than ",
            min_loci, " records, or at under ", 100 * .failed_call_rate, "% of the records where ",
            "the rest of their population is genotyped): ",
            paste(sprintf("%s (%s: %s records, %s)", failed$sample[shown], failed$population[shown],
                          .big(failed$called[shown]),
                          ifelse(is.na(failed$call_rate[shown]), "no call rate",
                                 sprintf("%.0f%% call rate", 100 * failed$call_rate[shown]))),
                  collapse = ", "),
            if (nrow(failed) > 10) sprintf(", and %d more", nrow(failed) - 10) else "", ".\n",
            "  They are KEPT. What that does: (1) each one's heterozygosity rests on few loci and ",
            "is noisy, yet it counts as a full individual in population means and tests; (2) heavy ",
            "missing data usually comes with allele dropout, which biases its heterozygosity low ",
            "and its F high; (3) it lowers its population's call rate at every record, so fewer ",
            "loci clear min_call.\n  If they were not meant to be analyzed, remove them from the ",
            "popmap (or with filter_samples()) and re-run.", call. = FALSE)
  }
  invisible(failed$sample)
}

## Not exported. Largest share of missing genotypes, among the records where
## the rest of its population is genotyped, before an individual is flagged
## as a failed library (the dDocent tutorial's individual cutoff).
.failed_call_rate <- 0.5

## Not exported. Stops, naming them, if any population has no record
## genotyped in at least `min_call` of its own individuals: nothing about that
## population can be estimated at this min_call.
.check_population_loci <- function(H, pops, locus_sets, min_call) {
  empty <- names(pops)[colSums(locus_sets) == 0L]
  if (!length(empty)) return(invisible(NULL))
  call_rate <- sweep(.typed_by_pop(H, pops[empty]), 2L, lengths(pops[empty]), "/")
  best <- apply(call_rate, 2L, max)
  stop("No record is genotyped in >= ", 100 * min_call, "% of the individuals of population(s): ",
       paste(sprintf("%s (highest call rate %.0f%%)", empty, 100 * best), collapse = ", "),
       ".\n  Lower min_call to at most the highest call rate shown, or remove these ",
       "populations from the popmap.", call. = FALSE)
}

## Not exported. The individual-level part of het_between_pops():
##   1. keep loci genotyped in >= min_call of ALL individuals pooled (used
##      only for this table and the missingness check; see ?het_between_pops);
##   2. heterozygosity of each individual = heterozygous loci / called loci.
## Every individual is kept (see .warn_failed_individuals()). One not called
## at any pooled locus has heterozygosity NA here, and is simply left out of
## the missingness check.
## Returns list(table, call_rate, het, n_loci_pooled).
.individual_het <- function(H, pops, min_call, min_loci, verbose) {
  ids <- unlist(pops, use.names = FALSE)
  pop_of_id <- rep(names(pops), lengths(pops))
  a1 <- H$A1[, ids, drop = FALSE]
  a2 <- H$A2[, ids, drop = FALSE]

  pooled_call_rate <- rowMeans(!is.na(a1))
  pooled <- pooled_call_rate >= min_call - .threshold_tol
  .inform(verbose, sprintf("Loci genotyped in >= %.0f%% of individuals (pooled): %s of %s",
                           100 * min_call, .big(sum(pooled)), .big(nrow(a1))))
  if (sum(pooled) < min_loci)
    .inform(verbose, "  WARNING: fewer than ", min_loci, " loci pass the pooled filter. ",
            "The per-individual table and the missingness confound check use this ",
            "pooled set and may be unreliable. The population summary and the tests ",
            "use each population's own loci instead.")
  a1 <- a1[pooled, , drop = FALSE]
  a2 <- a2[pooled, , drop = FALSE]

  n_called <- colSums(!is.na(a1))
  het <- colMeans(a1 != a2, na.rm = TRUE)
  het[n_called == 0] <- NA_real_                 # 0/0: no pooled locus called
  call_rate <- if (sum(pooled)) unname(n_called) / sum(pooled) else rep(NA_real_, length(ids))

  table <- data.frame(sample = ids, population = pop_of_id,
                      loci_called = as.integer(n_called), heterozygosity = unname(het),
                      row.names = NULL)
  list(table = table, call_rate = call_rate, het = unname(het), n_loci_pooled = sum(pooled))
}

## Not exported. Per population, on its own locus set: every individual's
## heterozygosity and inbreeding F (.ind_F(), R/individual_inbreeding.R), and
## the summary table. An individual not called at any of these loci has no
## heterozygosity (NA) and is left out of the summary: `n` counts the
## individuals with a value.
.population_summary <- function(H, pops, locus_sets) {
  het_by_pop <- F_by_pop <- stats::setNames(vector("list", length(pops)), names(pops))
  for (p in names(pops)) {
    loci <- locus_sets[, p]
    a1 <- H$A1[loci, pops[[p]], drop = FALSE]
    a2 <- H$A2[loci, pops[[p]], drop = FALSE]
    het <- colMeans(a1 != a2, na.rm = TRUE)
    het[!is.finite(het)] <- NA_real_
    het_by_pop[[p]] <- het
    F_by_pop[[p]] <- stats::setNames(.ind_F(a1, a2)$F, pops[[p]])
  }
  table <- do.call(rbind, lapply(names(pops), function(p) {
    h <- het_by_pop[[p]][!is.na(het_by_pop[[p]])]
    of_h <- function(f) if (length(h)) f(h) else NA_real_
    data.frame(population = p, n = length(h), n_loci = sum(locus_sets[, p]),
               mean_het = of_h(mean), sd = of_h(stats::sd),
               se = of_h(stats::sd) / sqrt(length(h)), min = of_h(min), max = of_h(max),
               mean_F = if (any(is.finite(F_by_pop[[p]]))) mean(F_by_pop[[p]], na.rm = TRUE)
                        else NA_real_)
  }))
  list(table = table, het_by_pop = het_by_pop,
       F_by_individual = unlist(unname(F_by_pop)))
}

## Not exported. MISSINGNESS CONFOUND. An individual's heterozygosity is
## computed over the loci it was called at. If poorly sequenced individuals
## are both more often missing and (through allele dropout) more often called
## homozygous, a difference between populations can be a difference in
## library quality. Only individuals with a heterozygosity on the pooled loci
## take part. Returns list(no_variation, too_few, overall, within):
##   no_variation  call rate does not vary among the individuals that take part
##                 (all called at every pooled locus; an individual with no
##                 heterozygosity, such as a failed library, does not take part)
##   too_few       fewer than 3 individuals called at any pooled locus
##   overall  one row: correlation r of call rate with heterozygosity, its
##            p-value, n, and the gap in mean call rate between populations
##   within   the same correlation within each population of > 3 individuals
.missingness_confound <- function(call_rate, het, pop_of_individual, pop_names) {
  has_value <- is.finite(call_rate) & is.finite(het)
  if (sum(has_value) < 3L)
    return(list(no_variation = FALSE, too_few = TRUE, overall = NULL, within = NULL))
  call_rate <- call_rate[has_value]
  het <- het[has_value]
  pop_of_individual <- pop_of_individual[has_value]
  if (stats::sd(call_rate) < .zero_tol)
    return(list(no_variation = TRUE, too_few = FALSE, overall = NULL, within = NULL))
  test <- suppressWarnings(stats::cor.test(call_rate, het))
  within <- do.call(rbind, lapply(pop_names, function(p) {
    rows <- which(pop_of_individual == p)
    if (length(rows) > 3 && stats::sd(call_rate[rows]) > .zero_tol)
      data.frame(population = p, r = suppressWarnings(stats::cor(call_rate[rows], het[rows])),
                 min_call_rate = min(call_rate[rows]), max_call_rate = max(call_rate[rows]))
  }))
  overall <- data.frame(r = unname(test$estimate), p_value = test$p.value, n = length(het),
                        call_rate_gap = diff(range(tapply(call_rate, pop_of_individual, mean))))
  list(no_variation = FALSE, too_few = FALSE, overall = overall, within = within)
}

## Not exported. Identity disequilibrium g2 in each population, on its own
## locus set (.g2_summary(), R/identity_disequilibrium.R): the same excess
## variation in individual heterozygosity, in its standard, citable form.
.g2_by_population <- function(H, pops, locus_sets, nboot) {
  do.call(rbind, lapply(names(pops), function(p) {
    loci <- locus_sets[, p]
    g2 <- .g2_summary(H$A1[loci, pops[[p]], drop = FALSE],
                      H$A2[loci, pops[[p]], drop = FALSE], nboot = nboot, nperm = 0L)
    data.frame(population = p, n = g2$n_ind, g2 = g2$g2, g2_lo = g2$g2_lo, g2_hi = g2$g2_hi)
  }))
}

## Not exported. THE LOCUS PART OF THE UNCERTAINTY in a difference between
## two populations' mean individual values, when each individual's value is a
## ratio over loci: heterozygosity = heterozygous loci / typed loci, and
## F = 1 - observed / expected heterozygous loci.
##
## WHY IT IS NEEDED. Welch's t counts how much individuals vary. It treats
## the loci as fixed. But the loci are a sample of the genome, and two
## populations with the same genome-wide heterozygosity still differ a little
## at any particular set of loci, because drift acts on each locus
## separately. That chance difference does not shrink relative to Welch's
## variance as loci are added. In simulations where each population had its
## own allele frequencies (inst/sims/het_test_null.R), Welch's t alone gave
## up to 14% false positives at FST 0.05 and up to 30% at FST 0.2, at a
## nominal 5%, when individuals were alike in inbreeding.
##
## HOW. A delete-one-RAD-locus jackknife of the difference, in its linear
## form. Leaving out RAD locus b changes individual i's ratio by about
##   r_bi = (num_bi - ratio_i * den_bi) / den_i
## (num_bi, den_bi: i's numerator and denominator within locus b; den_i: its
## total), so the difference changes by u_b = mean_i r_bi (population 1)
## minus the same for population 2, and the jackknife variance is
## (B - 1)/B * sum_b u_b^2 over B RAD loci.
##
## NOT COUNTING GENOTYPE NOISE TWICE. u_b also contains the chance of which
## genotypes the sampled individuals happen to carry at locus b, and Welch's
## variance already counts that. Its size is estimated from the spread of r_bi
## among individuals, var_i(r_bi) / n, and subtracted (Owen 2007 discusses
## this double counting). What is left is the variance from which loci were
## typed. It is set to 0 if negative.
##
##   num1, den1, num2, den2  loci x individuals matrices (0 where untyped),
##                           one pair per population, only the individuals
##                           that have a value
##   block                   RAD locus of each row
## Returns list(var = locus variance of the difference, n_blocks = RAD loci).
.locus_variance_of_difference <- function(num1, den1, num2, den2, block) {
  block <- match(block, unique(block))
  n_blocks <- max(0L, block)
  if (n_blocks < 2L) return(list(var = NA_real_, n_blocks = n_blocks))
  influence <- function(num, den) {
    N <- rowsum(num, block, reorder = TRUE)            # RAD loci x individuals
    D <- rowsum(den, block, reorder = TRUE)
    total_den <- colSums(D)
    ratio <- colSums(N) / total_den
    R <- sweep(N - sweep(D, 2L, ratio, "*"), 2L, total_den, "/")
    n <- ncol(R)
    list(mean = rowMeans(R),
         noise = if (n >= 2L) apply(R, 1L, stats::var) / n else rep(NA_real_, nrow(R)))
  }
  one <- influence(num1, den1)
  two <- influence(num2, den2)
  u <- one$mean - two$mean
  v <- (n_blocks - 1) / n_blocks * sum(u^2 - one$noise - two$noise)
  list(var = if (is.finite(v)) max(0, v) else NA_real_, n_blocks = n_blocks)
}

## Not exported. Tests every pair of populations, on the loci where BOTH
## populations of the pair clear min_call on their own (not the pooled set,
## and not "every population in the run"). Each row carries its own n_loci,
## n1, mean1, n2, mean2. The combined test (individuals and loci; see
## .locus_variance_of_difference()) is primary; Welch's t and Wilcoxon are
## kept for comparison. With k populations there are k(k-1)/2 pairs, so
## p-values are also Benjamini-Hochberg adjusted across pairs.
## Returns list(heterozygosity = table, F = table).
.pairwise_tests <- function(H, pops, locus_sets, min_call, min_loci, verbose) {
  pairs <- utils::combn(names(pops), 2, simplify = FALSE)
  rows <- lapply(pairs, function(pair) {
    p1 <- pair[1]
    p2 <- pair[2]
    loci <- locus_sets[, p1] & locus_sets[, p2]
    n_loci <- sum(loci)
    if (n_loci < min_loci) {
      .inform(verbose, sprintf(
        "  %s vs %s: only %d loci clear %.0f%% call rate in BOTH populations ",
        p1, p2, n_loci, 100 * min_call),
        "(need >= ", min_loci, "). Reported as NA for this pair; other pairs are unaffected.")
      na_row <- .na_test_row(p1, p2, n_loci, length(pops[[p1]]), length(pops[[p2]]))
      return(list(het = na_row, F = na_row))
    }
    block <- H$locus_raw[loci]
    a1_1 <- H$A1[loci, pops[[p1]], drop = FALSE]
    a2_1 <- H$A2[loci, pops[[p1]], drop = FALSE]
    a1_2 <- H$A1[loci, pops[[p2]], drop = FALSE]
    a2_2 <- H$A2[loci, pops[[p2]], drop = FALSE]
    het_1 <- colMeans(a1_1 != a2_1, na.rm = TRUE)
    het_2 <- colMeans(a1_2 != a2_2, na.rm = TRUE)
    has_1 <- is.finite(het_1)
    has_2 <- is.finite(het_2)
    ## Heterozygosity: heterozygous loci / typed loci, per individual.
    het_cells <- function(a1, a2, keep) {
      typed <- !is.na(a1[, keep, drop = FALSE])
      list(num = (typed & a1[, keep, drop = FALSE] != a2[, keep, drop = FALSE]) * 1,
           den = typed * 1)
    }
    h1 <- het_cells(a1_1, a2_1, has_1)
    h2 <- het_cells(a1_2, a2_2, has_2)
    het_locus <- .locus_variance_of_difference(h1$num, h1$den, h2$num, h2$den, block)
    ## Individual F on the same loci and individuals, with expected
    ## heterozygosity from each population's own allele frequencies.
    f1 <- .ind_F_cells(a1_1, a2_1)
    f2 <- .ind_F_cells(a1_2, a2_2)
    F_1 <- .ind_F(a1_1, a2_1)$F[has_1]
    F_2 <- .ind_F(a1_2, a2_2)$F[has_2]
    ok_1 <- which(has_1)[is.finite(F_1)]
    ok_2 <- which(has_2)[is.finite(F_2)]
    F_locus <- .locus_variance_of_difference(f1$observed[, ok_1, drop = FALSE],
                                             f1$expected[, ok_1, drop = FALSE],
                                             f2$observed[, ok_2, drop = FALSE],
                                             f2$expected[, ok_2, drop = FALSE], block)
    list(het = .two_sample(het_1, het_2, p1, p2, n_loci, "heterozygosity", verbose, het_locus),
         F = .two_sample(F_1, F_2, p1, p2, n_loci, "F", verbose, F_locus))
  })
  bh_adjust <- function(tab) {
    tab$p_combined_BH <- stats::p.adjust(tab$p_combined, method = "BH")
    tab$p_welch_BH <- stats::p.adjust(tab$p_welch, method = "BH")
    tab$p_wilcox_BH <- stats::p.adjust(tab$p_wilcox, method = "BH")
    tab
  }
  list(heterozygosity = bh_adjust(do.call(rbind, lapply(rows, `[[`, "het"))),
       F = bh_adjust(do.call(rbind, lapply(rows, `[[`, "F"))))
}

## Not exported. OVERALL TESTS for 3 or more populations: is there any
## difference among them, before looking at pairs? With k populations there
## are k(k - 1)/2 pairwise tests; an overall test answers the question once.
## Both use one number per individual, as the pairwise tests do:
##   Welch's one-way ANOVA (stats::oneway.test(var.equal = FALSE); Welch
##   1951), which allows unequal variances and reduces to Welch's t for two
##   groups;
##   Kruskal-Wallis (stats::kruskal.test; Kruskal & Wallis 1952), the
##   distribution-free check.
## Loci: those that clear min_call in EVERY population (the k-population
## version of the pairwise rule). Returns NULL for 2 populations, otherwise a
## data frame with one row for heterozygosity and one for F.
.omnibus_tests <- function(H, pops, locus_sets, min_call, min_loci, verbose) {
  if (length(pops) < 3L) return(NULL)
  loci <- rowSums(locus_sets) == ncol(locus_sets)
  n_loci <- sum(loci)
  empty_row <- function(what) {
    data.frame(statistic = what, n_loci = n_loci, n_individuals = NA_integer_,
               welch_F = NA_real_, welch_df1 = NA_real_, welch_df2 = NA_real_,
               p_welch = NA_real_, kruskal_chisq = NA_real_, kruskal_df = NA_real_,
               p_kruskal = NA_real_)
  }
  if (n_loci < min_loci) {
    .inform(verbose, sprintf(
      "  Overall test: only %d loci clear %.0f%% call rate in ALL %d populations ",
      n_loci, 100 * min_call, length(pops)),
      "(need >= ", min_loci, "); reported as NA. The pairwise tests are unaffected.")
    return(rbind(empty_row("heterozygosity"), empty_row("F")))
  }
  values <- lapply(names(pops), function(p) {
    a1 <- H$A1[loci, pops[[p]], drop = FALSE]
    a2 <- H$A2[loci, pops[[p]], drop = FALSE]
    het <- colMeans(a1 != a2, na.rm = TRUE)
    list(heterozygosity = het, F = .ind_F(a1, a2)$F)
  })
  one_test <- function(what) {
    y <- unlist(lapply(values, `[[`, what), use.names = FALSE)
    group <- factor(rep(names(pops), vapply(values, function(v) length(v[[what]]), integer(1))),
                    levels = names(pops))
    keep <- is.finite(y)
    y <- y[keep]
    group <- droplevels(group[keep])
    row <- empty_row(what)
    row$n_individuals <- length(y)
    if (nlevels(group) < 3L || any(table(group) < 2L)) {
      .inform(verbose, "  Overall test (", what, "): fewer than 2 individuals with a value ",
              "in some population; reported as NA.")
      return(row)
    }
    welch <- tryCatch(suppressWarnings(stats::oneway.test(y ~ group, var.equal = FALSE)),
                      error = function(e) NULL)
    if (!is.null(welch) && is.finite(welch$p.value)) {
      row$welch_F <- unname(welch$statistic)
      row$welch_df1 <- unname(welch$parameter[1])
      row$welch_df2 <- unname(welch$parameter[2])
      row$p_welch <- welch$p.value
    } else {
      .inform(verbose, "  Overall Welch test (", what, "): not computable, usually because a ",
              "population has no variation among individuals; reported as NA.")
    }
    kw <- tryCatch(stats::kruskal.test(y, group), error = function(e) NULL)
    if (!is.null(kw) && is.finite(kw$p.value)) {
      row$kruskal_chisq <- unname(kw$statistic)
      row$kruskal_df <- unname(kw$parameter)
      row$p_kruskal <- kw$p.value
    }
    row
  }
  rbind(one_test("heterozygosity"), one_test("F"))
}

## Not exported. One all-NA row of the pairwise test table.
.na_test_row <- function(p1, p2, n_loci, n1, n2, mean1 = NA_real_, mean2 = NA_real_) {
  data.frame(pop1 = p1, pop2 = p2, n_loci = n_loci, n1 = n1, mean1 = mean1,
             n2 = n2, mean2 = mean2, diff = NA_real_, se_combined = NA_real_,
             ci_lo = NA_real_, ci_hi = NA_real_, df = NA_real_, p_combined = NA_real_,
             p_welch = NA_real_, p_wilcox = NA_real_, hedges_g = NA_real_)
}

## Not exported. The tests for one pair of populations, on one number per
## individual (`a`, `b`: heterozygosity or F). Returns one row of the pairwise
## table.
##
##   combined (primary)  se_combined^2 = var(a)/n_a + var(b)/n_b + locus part
##                       (`locus`, from .locus_variance_of_difference()); a t
##                       test with Welch-Satterthwaite degrees of freedom over
##                       the three parts (the locus part has RAD loci - 1);
##                       ci_lo, ci_hi and p_combined come from it
##   Welch's t           the individual part alone (p_welch), for comparison
##   Wilcoxon            distribution-free, individuals only (p_wilcox)
##   Hedges' g           the difference in among-individual SD units
.two_sample <- function(a, b, p1, p2, n_loci, what, verbose,
                        locus = list(var = NA_real_, n_blocks = 0L)) {
  ## An individual can pass the pooled check yet be called at none of this
  ## pair's loci, and so have no value here: drop it from this pair only.
  if (any(!is.finite(a)) || any(!is.finite(b))) {
    .inform(verbose, sprintf(
      "  %s vs %s (%s): %d individual(s) without a value on the %d shared loci ",
      p1, p2, what, sum(!is.finite(a)) + sum(!is.finite(b)), n_loci),
      "for this pair; excluded from this comparison only.")
    a <- a[is.finite(a)]
    b <- b[is.finite(b)]
  }
  if (length(a) < 2 || length(b) < 2) {
    .inform(verbose, sprintf(
      "  %s vs %s (%s): fewer than 2 individuals with a value in one group. ",
      p1, p2, what), "Reported as NA for this pair.")
    return(.na_test_row(p1, p2, n_loci, length(a), length(b),
                        if (length(a)) mean(a) else NA_real_,
                        if (length(b)) mean(b) else NA_real_))
  }
  ## "No variation among individuals, so no test is possible" shows up in two
  ## ways: t.test() stops with an error (both samples constant at the same
  ## nonzero value), or it returns NaN (both constant at exactly 0, where its
  ## zero-variance check compares against 0 and does not fire). Both become NA.
  no_test <- list(p.value = NA_real_, conf.int = c(NA_real_, NA_real_))
  welch <- tryCatch(stats::t.test(a, b), error = function(e) {
    .inform(verbose, sprintf("  %s vs %s (%s): t.test() failed (%s). ", p1, p2, what,
                             trimws(conditionMessage(e))),
            "This usually means no variation among individuals in one or both ",
            "populations, so no test is possible. Reported as NA.")
    no_test
  })
  if (is.nan(welch$p.value)) {
    .inform(verbose, sprintf("  %s vs %s (%s): t.test() returned NaN (both samples constant at exactly ",
                             p1, p2, what), "zero). No variation to test. Reported as NA.")
    welch <- no_test
  }
  ## `exact` is set explicitly so the p-value does not depend on the R
  ## version: R 4.6.0 changed the default with tied values from the normal
  ## approximation to exact conditional inference. This is the rule R used
  ## before 4.6.0: exact when there are no ties and both groups have fewer
  ## than 50 values, otherwise the normal approximation with continuity
  ## correction.
  exact <- !anyDuplicated(c(a, b)) && length(a) < 50 && length(b) < 50
  wilcox <- tryCatch(suppressWarnings(stats::wilcox.test(a, b, exact = exact, correct = TRUE)),
                     error = function(e) no_test)
  if (is.nan(wilcox$p.value)) wilcox <- no_test

  ## Hedges' g: the standardised difference with a small-sample correction J.
  n_a <- length(a)
  n_b <- length(b)
  pooled_sd <- sqrt(((n_a - 1) * stats::var(a) + (n_b - 1) * stats::var(b)) / (n_a + n_b - 2))
  cohen_d <- if (is.finite(pooled_sd) && pooled_sd > 0) (mean(a) - mean(b)) / pooled_sd else NA_real_
  small_sample_correction <- 1 - 3 / (4 * (n_a + n_b) - 9)

  ## The combined test: individuals (as in Welch's t) plus the loci.
  diff <- mean(a) - mean(b)
  v_a <- stats::var(a) / n_a
  v_b <- stats::var(b) / n_b
  v_locus <- locus$var
  se_combined <- df <- p_combined <- ci_lo <- ci_hi <- NA_real_
  if (is.finite(welch$p.value) && is.finite(v_locus)) {
    se_combined <- sqrt(v_a + v_b + v_locus)
    df_parts <- v_a^2 / (n_a - 1) + v_b^2 / (n_b - 1) +
      (if (v_locus > 0) v_locus^2 / (locus$n_blocks - 1) else 0)
    df <- (v_a + v_b + v_locus)^2 / df_parts
    p_combined <- 2 * stats::pt(-abs(diff / se_combined), df)
    ci_lo <- diff - stats::qt(0.975, df) * se_combined
    ci_hi <- diff + stats::qt(0.975, df) * se_combined
  } else if (is.finite(welch$p.value)) {
    .inform(verbose, sprintf("  %s vs %s (%s): fewer than 2 RAD loci, so the locus part of the ",
                             p1, p2, what), "combined test cannot be estimated. Reported as NA.")
  }

  data.frame(pop1 = p1, pop2 = p2, n_loci = n_loci,
             n1 = n_a, mean1 = mean(a), n2 = n_b, mean2 = mean(b),
             diff = diff, se_combined = se_combined, ci_lo = ci_lo, ci_hi = ci_hi, df = df,
             p_combined = p_combined, p_welch = welch$p.value, p_wilcox = wilcox$p.value,
             hedges_g = cohen_d * small_sample_correction)
}

## ---------------------------------------------------------------------------
## Rounding, print() and summary()
## ---------------------------------------------------------------------------

## Not exported. How het_between_pops() tables are rounded when printed or
## written (the returned object keeps full precision).
.het_rounding <- list(
  individual_heterozygosity = list(digits = 4, round_cols = c(heterozygosity = 5)),
  population_summary = list(digits = 4),
  omnibus = list(digits = 3, signif_cols = c(p_welch = 3, p_kruskal = 3)),
  pairwise_tests = list(digits = 4, round_cols = c(hedges_g = 2),
                        signif_cols = c(p_combined = 3, p_welch = 3, p_wilcox = 3,
                                        p_combined_BH = 3, p_welch_BH = 3, p_wilcox_BH = 3)),
  g2 = list(digits = 4)
)

## Not exported. Applies .het_rounding to one table.
.round_het_table <- function(df, table_name) {
  rules <- .het_rounding[[table_name]]
  .round_table(df, digits = rules$digits, round_cols = rules$round_cols,
               signif_cols = rules$signif_cols)
}

## Not exported. The pairwise table for display: at most `max_rows` rows,
## the most significant first when it has to be cut.
.show_pairwise <- function(tab, max_rows, full_name) {
  if (nrow(tab) <= max_rows) {
    .print_table(.round_het_table(tab, "pairwise_tests"))
  } else {
    cat(sprintf("  (%d most significant of %d pairs; the full table is %s)\n",
                max_rows, nrow(tab), full_name))
    .print_table(.round_het_table(tab[utils::head(order(tab$p_combined_BH), max_rows), ],
                                  "pairwise_tests"))
  }
}

#' @rdname het_between_pops
#' @param x,object A `raddiv_het` object, as returned by `het_between_pops()`.
#' @param ... Ignored.
#' @export
print.raddiv_het <- function(x, ...) {
  st <- x$settings
  cat(sprintf("Heterozygosity between populations (combined test: individuals and loci): %d populations, %d individuals\n",
              st$n_pops, st$n_individuals))
  cat(sprintf("  loci: each population's own (>= %.0f%% call rate within it); a pair uses loci both clear\n",
              100 * st$min_call))
  cat("\n$population_summary\n")
  .print_table(.round_het_table(x$population_summary, "population_summary"))
  if (!is.null(x$omnibus)) {
    cat("\n$omnibus (overall test across all populations: look at this before the pairs)\n")
    .print_table(.round_het_table(x$omnibus, "omnibus"))
  }
  cat("\n$pairwise_tests (heterozygosity)\n")
  .show_pairwise(x$pairwise_tests, 10, "x$pairwise_tests")
  cat("\n$pairwise_F_tests (individual inbreeding F)\n")
  .show_pairwise(x$pairwise_F_tests, 10, "x$pairwise_F_tests")
  cat("\n")
  for (w in .het_warnings(x)) cat("NOTE:", w, "\n")
  cat("Also in the result: $individual_heterozygosity, $g2, $missingness_confound\n")
  cat("summary() shows the two tests and what to check; summary(x, details = TRUE) has every column.\n")
  invisible(x)
}

## Not exported. The short warnings print.raddiv_het() shows.
.het_warnings <- function(x) {
  out <- character(0)
  flagged <- x$settings$flagged_individuals
  if (length(flagged))
    out <- c(out, sprintf("%d individual(s) look like failed libraries and were kept: %s%s (see the warning and ?het_between_pops).",
                          length(flagged), paste(utils::head(flagged, 5), collapse = ", "),
                          if (length(flagged) > 5) ", ..." else ""))
  cf <- x$missingness_confound
  if (!is.null(cf$overall) && is.finite(cf$overall$p_value) && cf$overall$p_value < 0.05 &&
      cf$overall$r > 0)
    out <- c(out, "call rate is positively correlated with heterozygosity (possible allele dropout).")
  if (!is.null(cf$overall) && cf$overall$call_rate_gap > 0.02)
    out <- c(out, "mean call rate differs between populations by more than 2 percentage points.")
  if (any(x$g2$g2_lo > 0, na.rm = TRUE))
    out <- c(out, "g2 > 0 in some population: individuals differ in inbreeding; report _se_combined from diversity_stats() (computed by default).")
  out
}

#' @rdname het_between_pops
#' @param details For `summary()`: `FALSE` (the default) prints a short view --
#'   the two tests (heterozygosity, and inbreeding F) as the combined test's
#'   difference with its SE, 95% interval and p-value, the populations with
#'   their g2, what to report, and a list of checks. `TRUE` prints the full
#'   report: every test column (including the Welch and Wilcoxon p-values), the
#'   missingness check in full and the notes on interpretation. Either way, the
#'   object returned holds the tables the short view prints, in `$tables` (see
#'   "Tables for a paper").
#' @export
summary.raddiv_het <- function(object, details = FALSE, ...) {
  .check_flag(details, "details")
  ## The tables the short view prints, with every pair (the short view shows at
  ## most 5 per question), so a table for a paper can be built from exactly what
  ## it shows.
  structure(list(result = object, details = details,
                 tables = list(
                   heterozygosity = .het_pair_table(object$pairwise_tests,
                                                    max_rows = nrow(object$pairwise_tests))$table,
                   inbreeding = .het_pair_table(object$pairwise_F_tests,
                                                max_rows = nrow(object$pairwise_F_tests))$table,
                   populations = .het_population_table(object))),
            class = "summary.raddiv_het")
}

#' @export
print.summary.raddiv_het <- function(x, ...) {
  if (isTRUE(x$details)) .het_report(x$result) else .het_brief(x$result)
  invisible(x)
}

## Not exported. One question's pairwise tests as a short table: the COMBINED
## test only (the one to report): the difference in means with its SE, its 95%
## interval, its p-value (Benjamini-Hochberg adjusted when there is more than
## one pair) and Hedges' g. Welch's and Wilcoxon's p-values are for comparison
## and stay in the full report. At most `max_rows` pairs (5 in the short view,
## so that it stays short with many populations), most significant first.
.het_pair_table <- function(tab, max_rows = 5L) {
  shown <- if (nrow(tab) <= max_rows) tab
           else tab[utils::head(order(tab$p_combined_BH), max_rows), ]
  multiple <- nrow(tab) > 1L
  d <- .se_decimals(shown$se_combined)
  out <- data.frame(pair = paste(shown$pop1, "-", shown$pop2), stringsAsFactors = FALSE,
                    check.names = FALSE)
  out[["difference (SE)"]] <- .est_se(shown$diff, shown$se_combined)
  out[["95% CI"]] <- sprintf("%.*f to %.*f", d, shown$ci_lo, d, shown$ci_hi)
  out[[if (multiple) "p (combined, BH)" else "p (combined)"]] <-
    .format_p(if (multiple) shown$p_combined_BH else shown$p_combined)
  out[["Hedges g"]] <- sprintf("%.2f", shown$hedges_g)
  list(table = out, n_pairs = nrow(tab), n_shown = nrow(shown))
}

## Not exported. The POPULATIONS table of the short summary: each population's
## mean individual heterozygosity (no SE: populations are compared with the
## tests, never by overlapping intervals), its spread among individuals, mean F,
## and g2 with its 95% interval.
.het_population_table <- function(x) {
  ps <- x$population_summary
  g2 <- x$g2[match(ps$population, x$g2$population), ]
  data.frame(population = ps$population, n = ps$n, n_loci = ps$n_loci,
             mean_het = sprintf("%.4f", ps$mean_het), sd = sprintf("%.4f", ps$sd),
             mean_F = sprintf("%.4f", ps$mean_F),
             `g2 (95% CI)` = sprintf("%.4f (%.4f to %.4f)", g2$g2, g2$g2_lo, g2$g2_hi),
             stringsAsFactors = FALSE, check.names = FALSE)
}

## Not exported. The overall (omnibus) test for one question, as one sentence
## (3 or more populations only).
.het_omnibus_line <- function(x, statistic, full = FALSE) {
  om <- x$omnibus
  if (is.null(om)) return(NULL)
  row <- om[om$statistic == statistic, , drop = FALSE]
  if (!nrow(row)) return(NULL)
  sprintf("Overall test across all %d populations (Welch's one-way ANOVA): p = %s. Read the pairs only if it is significant%s.",
          x$settings$n_pops, .format_p(row$p_welch),
          if (full) "; it counts individuals only, so confirm a pair by its combined p" else "")
}

## Not exported. The checks of a het_between_pops() result as .check() lists:
## failed libraries, the missingness confound, the call-rate gap, and g2.
.het_checks <- function(x) {
  st <- x$settings
  checks <- list()
  add <- function(tag, text) checks[[length(checks) + 1L]] <<- .check(tag, text)

  flagged <- st$flagged_individuals
  if (length(flagged))
    add("look", sprintf("KEPT, but they look like failed libraries: %s. Each rests on few loci yet counts as a full individual; remove them from the popmap if they were not meant to be analyzed",
                        paste(flagged, collapse = ", ")))
  else add("ok", "no individual looks like a failed library")

  cf <- x$missingness_confound
  if (isTRUE(cf$too_few)) {
    add("info", "fewer than 3 individuals were called at any pooled locus; the missingness check needs more")
  } else if (isTRUE(cf$no_variation)) {
    add("ok", "call rate does not vary among the individuals checked, so it cannot confound the test")
  } else if (!is.null(cf$overall)) {
    if (is.finite(cf$overall$p_value) && cf$overall$p_value < 0.05 && cf$overall$r > 0)
      add("look", "individuals with more missing data look LESS heterozygous (the allele-dropout signature); re-run with min_call = 1.0 before quoting the test")
    else
      add("ok", sprintf("heterozygosity does not track call rate (r = %+.2f, p = %.2g)",
                        cf$overall$r, cf$overall$p_value))
    if (cf$overall$call_rate_gap > 0.02)
      add("look", "mean call rate differs between populations by more than 2 percentage points; that asymmetry is itself a candidate explanation")
  }

  g2 <- x$g2
  above <- g2$population[which(g2$g2_lo > 0)]
  if (length(above))
    add("look", sprintf("individuals differ in inbreeding in %s (g2 above 0), so SEs over loci alone are too small for Ho and Fis: report _se_combined from diversity_stats() (computed by default)",
                        paste(above, collapse = " and ")))
  else
    add("ok", "g2 is not above 0 in any population: no sign that individuals differ in inbreeding")
  checks
}

## Not exported. What to report, for the short and the full view.
.het_what_to_report <- function(x) {
  c("Report the combined test: the difference, its SE and 95% CI, and p.",
    "Heterozygosity reflects diversity and inbreeding together. Use the F test",
    "for a question about inbreeding, and diversity_stats() to compare He itself.",
    if (!is.null(x$omnibus))
      "With 3+ populations, read the pairs only if the overall test is significant.")
}

## Not exported. The short summary(): the two tests, the populations, what to
## report, and the checks. `summary(x, details = TRUE)` prints the full report.
.het_brief <- function(x) {
  st <- x$settings
  cat(sprintf("HETEROZYGOSITY BETWEEN POPULATIONS: %d individuals, %d populations\n",
              st$n_individuals, st$n_pops))
  .legend(sprintf("Loci: each population's own (at least %.0f%% call rate within it);", 100 * st$min_call),
          "a pair uses the loci both populations clear.")

  question <- function(title, tab, statistic, after = NULL) {
    .section(title, "   combined test")
    .legend(strwrap(.het_omnibus_line(x, statistic), width = 76))
    r <- .het_pair_table(tab)
    .print_compact(r$table)
    if (r$n_shown < r$n_pairs)
      .legend(sprintf("Showing the %d most significant of %d pairs;", r$n_shown, r$n_pairs),
              sprintf("the full table is x$%s.",
                      if (statistic == "F") "pairwise_F_tests" else "pairwise_tests"))
    .legend(after)
  }
  question("DO THE POPULATIONS DIFFER IN HETEROZYGOSITY?", x$pairwise_tests, "heterozygosity",
           strwrap(sprintf("The combined test counts individuals and loci. In simulations with no true difference it wrongly rejected %s of the time (Welch's t alone: up to 30%%).",
                           .coverage$het_test), width = 76))
  question("DO THEY DIFFER IN INBREEDING (F)?", x$pairwise_F_tests, "F",
           c("F = 1 - observed / expected heterozygosity, with the expected value taken",
             "from each individual's own population."))

  .section("POPULATIONS   Identity disequilibrium g2: 0 when individuals do not differ")
  .print_compact(.het_population_table(x))
  .legend("mean_het has no SE here. Compare populations with the tests above, never",
          "by overlapping intervals. sd is the spread among individuals.")

  .section("WHAT TO REPORT")
  .legend(.het_what_to_report(x))
  .section("CHECKS")
  .check_lines(.het_checks(x))
  cat("\nNot shown here: Welch's and Wilcoxon's p-values, the missingness check in\n",
      "full, and the notes on interpretation.\n", sep = "")
  .legend("summary(x, details = TRUE)   explains every column and check",
          "x$pairwise_tests, x$pairwise_F_tests   hold every test column",
          "x$population_summary, x$g2             hold the population tables")
  invisible(NULL)
}

## Not exported. The pairwise tests as small labelled tables: the means, the
## combined test (the one to report), and the tests kept for comparison.
.het_pair_blocks <- function(tab, max_rows, full_name) {
  shown <- if (nrow(tab) <= max_rows) tab
           else tab[utils::head(order(tab$p_combined_BH), max_rows), ]
  if (nrow(shown) < nrow(tab))
    cat(sprintf("  (%d most significant of %d pairs; the full table is %s)\n",
                max_rows, nrow(tab), full_name))
  r <- .round_het_table(shown, "pairwise_tests")
  r <- cbind(pair = paste(r$pop1, "-", r$pop2), r[, setdiff(names(r), c("pop1", "pop2")), drop = FALSE])
  cat("Means\n")
  .print_compact(r[, c("pair", "n_loci", "n1", "mean1", "n2", "mean2", "diff")])
  cat("Combined test (the one to report)\n")
  .print_compact(r[, c("pair", "diff", "se_combined", "ci_lo", "ci_hi", "df", "p_combined",
                       "p_combined_BH")])
  cat("For comparison\n")
  .print_compact(r[, c("pair", "hedges_g", "p_welch", "p_welch_BH", "p_wilcox", "p_wilcox_BH")])
  invisible(NULL)
}

## Not exported. The full heterozygosity report (summary(x, details = TRUE), and
## what the command-line script prints with --details). The same two tests as
## the short summary come first, then every column.
.het_report <- function(x) {
  st <- x$settings
  note <- function(...) cat(paste0("  ", c(...), "\n"), sep = "")
  tests <- x$pairwise_tests
  rule <- strrep("=", 69)

  cat("\n", rule, "\n", sep = "")
  cat("  HETEROZYGOSITY BETWEEN POPULATIONS -- one value per individual; tests count\n")
  cat("  both which individuals and which loci were sampled\n")
  cat("  ", .big(st$n_loci_pooled), " loci pass the pooled filter (used for the\n",
      "  individual-level table and the confound check below);\n",
      "  ", st$n_individuals, " individuals, ", st$n_pops, " populations\n", sep = "")
  cat(rule, "\n", sep = "")
  if (length(st$flagged_individuals))
    note("", sprintf("KEPT, but they look like failed libraries: %s.",
                     paste(st$flagged_individuals, collapse = ", ")),
         .report_text$het_failed_individuals)

  ## The two tests, as in the short summary.
  for (q in list(list("DO THE POPULATIONS DIFFER IN HETEROZYGOSITY?", tests, "heterozygosity"),
                 list("DO THEY DIFFER IN INBREEDING (F)?", x$pairwise_F_tests, "F"))) {
    .section(q[[1]], "   combined test")
    .legend(strwrap(.het_omnibus_line(x, q[[3]], full = TRUE), width = 76))
    r <- .het_pair_table(q[[2]], max_rows = 10L)
    .print_compact(r$table)
  }
  cat("\n", strrep("-", 69), "\n", sep = "")
  note(.het_what_to_report(x))
  cat(strrep("-", 69), "\n", sep = "")

  cat("\nFULL TABLES\n")
  cat("\nPer-population summary of individual heterozygosity\n")
  cat("  (each population's OWN locus set: loci at >= ", 100 * st$min_call,
      "% call rate\n   WITHIN that population, see 'n_loci'.",
      " A pair's test below uses only loci BOTH\n",
      "   members of that pair clear, so its loci can differ from these; that is\n",
      "   expected, not an inconsistency.)\n", sep = "")
  .print_compact(.round_het_table(x$population_summary, "population_summary"))

  cat("\nMissingness confound check\n")
  cf <- x$missingness_confound
  if (isTRUE(cf$too_few)) {
    note("Fewer than 3 individuals were called at any pooled locus; the check needs more.")
  } else if (cf$no_variation) {
    note("Call rate does not vary among the individuals checked; no confound possible.")
  } else {
    cat(sprintf("  cor(per-individual call rate, heterozygosity) = %+.3f  (p = %.3g, n = %d)\n",
                cf$overall$r, cf$overall$p_value, cf$overall$n))
    for (k in seq_len(NROW(cf$within)))
      cat(sprintf("    within %-12s r = %+.3f   call rate %.3f-%.3f\n", cf$within$population[k],
                  cf$within$r[k], cf$within$min_call_rate[k], cf$within$max_call_rate[k]))
    if (is.finite(cf$overall$p_value) && cf$overall$p_value < 0.05 && cf$overall$r > 0)
      note(.report_text$het_confound_positive)
    else
      note("No positive association, so heterozygosity is not tracking coverage.")
    if (cf$overall$call_rate_gap > 0.02)
      note(.report_text$het_call_rate_gap)
  }
  cat("\n")
  note(.report_text$het_sd_among_individuals)

  cat("\nIdentity disequilibrium g2 (David et al. 2007) -- the standard measure of\n")
  cat("variance in inbreeding among individuals; 0 when individuals do not differ.\n")
  .print_compact(.round_het_table(x$g2, "g2"))
  note(.report_text$het_g2)
  note(.report_text$het_van_dongen)

  if (!is.null(x$omnibus)) {
    cat(sprintf("\nOverall test across all %d populations (loci that clear min_call in every one)\n",
                st$n_pops))
    .print_compact(.round_het_table(x$omnibus, "omnibus"))
    note(.report_text$het_omnibus)
  }

  n_pairs <- nrow(tests)
  cat(sprintf("\nPairwise tests, heterozygosity -- %d comparison%s\n", n_pairs,
              if (n_pairs == 1) "" else "s"))
  note("The combined test (p_combined, ci_lo, ci_hi) is primary: it counts both which",
       "individuals and which loci were sampled. p_welch (individuals only) and",
       "p_wilcox (distribution-free) are for comparison.",
       if (n_pairs > 1)
         sprintf("p_*_BH are Benjamini-Hochberg adjusted across all %d comparisons.", n_pairs))
  .het_pair_blocks(tests, 25, "x$pairwise_tests")
  n_significant <- function(tab) {
    c(combined = sum(tab$p_combined_BH < 0.05, na.rm = TRUE),
      welch = sum(tab$p_welch_BH < 0.05, na.rm = TRUE),
      wilcox = sum(tab$p_wilcox_BH < 0.05, na.rm = TRUE))
  }
  sig <- n_significant(tests)
  cat(sprintf("\n  pairs significant at BH < 0.05: combined %d (Welch %d, Wilcoxon %d), of %d\n",
              sig[["combined"]], sig[["welch"]], sig[["wilcox"]], n_pairs))
  note("Hedges' g is the standardised difference: ~0.2 small, ~0.5 medium, ~0.8 large.")

  F_tests <- x$pairwise_F_tests
  cat("\n")
  cat(paste0(.report_text$het_F_tests, "\n"), sep = "")
  .het_pair_blocks(F_tests, 25, "x$pairwise_F_tests")
  sig_F <- n_significant(F_tests)
  cat(sprintf("\n  pairs significant at BH < 0.05: combined %d (Welch %d, Wilcoxon %d), of %d\n",
              sig_F[["combined"]], sig_F[["welch"]], sig_F[["wilcox"]], nrow(F_tests)))

  cat("\nInterpretation notes\n")
  note(.report_text$het_interpretation)
  if (length(st$files)) cat("\nWrote ", paste(st$files, collapse = " and "), "\n", sep = "")
  cat("\n")
  invisible(NULL)
}
