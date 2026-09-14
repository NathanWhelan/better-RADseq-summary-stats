###############################################################################
#
#  R/het_between_pops.R -- does mean heterozygosity differ between populations?
#  Each individual contributes ONE number, so individuals are the replicate.
#
#  WHY NOT THE USUAL TEST. Bootstrapping over loci, or a paired Wilcoxon on
#  per-locus values, treats LOCI as the replicate. Loci are repeated measures
#  on the same individuals, so the interval shrinks as 1/sqrt(n_loci) while
#  the real uncertainty -- which individuals were caught -- does not shrink at
#  all. In het_between_pops_selftest()'s simulations of two populations with
#  identical true heterozygosity, a locus bootstrap rejected about 8% of the
#  time when individuals did not differ in inbreeding, but about half the
#  time when the SD of F among individuals was 0.10 and about 80% at 0.25.
#  Welch's t on per-individual heterozygosity stayed near 5% throughout. The
#  objection is old (Van Dongen 1995). Full discussion in
#  vignette("rationale"), section 5.
#
#  WHAT TO READ IN THE REPORT, IN ORDER
#    1. the missingness confound check -- if call rate correlates with
#       heterozygosity, the test measures library quality, not biology
#    2. the overdispersion factor -- how badly a locus bootstrap would have
#       misled on THESE data
#    3. the tests themselves: Welch's t, with Wilcoxon as a distribution-free
#       check
#
#  Needs the VCF: per-individual heterozygosity cannot be recovered from
#  populations.sumstats.tsv, which holds only per-population counts.
#
#  STEPS (functions below the exported one):
#    .individual_het()        pooled locus filter, per-individual table
#    .population_locus_sets() which loci each population uses (R/vcf_io.R)
#    .population_summary()    per-population mean heterozygosity and F
#    .missingness_confound()  does call rate track heterozygosity?
#    .overdispersion()        excess variance among individuals
#    .g2_by_population()      identity disequilibrium
#    .pairwise_tests()        Welch / Wilcoxon / Hedges' g per pair
#  Checked by tests/testthat/test-diversity-stats.R, test-regression.R
#  (golden values), test-individual-inbreeding.R and
#  het_between_pops_selftest().
#
###############################################################################

#' Test whether mean heterozygosity differs between populations
#'
#' Uses the INDIVIDUAL, not the locus, as the unit of replication (Welch's t,
#' with Wilcoxon as a distribution-free check), and adds a missingness-confound
#' check, an overdispersion check showing how badly a locus bootstrap would
#' have misled on this dataset, and identity disequilibrium (g2) per
#' population. Printing the result shows the main tables; `summary()` shows the
#' full report; `outdir` also writes the tables to TSV files.
#'
#' @details
#' **Three or more populations.** Start with `omnibus`, an overall test of
#' whether any population differs (Welch's one-way ANOVA, with Kruskal-Wallis
#' as the distribution-free check). Interpret the pairwise tests only if the
#' overall test is significant; they then say which pairs differ.
#'
#' Two questions, same machinery. `pairwise_tests` compares mean individual
#' heterozygosity: do the populations differ in diversity?
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
#' rule would let a well-covered population carry a poorly covered one.
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
#'   loci is excluded (a failed library cannot give a heterozygosity), a pair
#'   of populations sharing fewer is reported as `NA`, and a warning is given
#'   if fewer pass the pooled filter. Default `50`.
#' @param nboot_g2 Bootstrap replicates over individuals for the g2
#'   confidence interval. Default `1000`, as in [identity_disequilibrium()].
#' @param seed Random seed for the g2 bootstrap. Default `NULL`: use R's
#'   current random-number stream (call `set.seed()` first for a reproducible
#'   interval). A number makes the result reproducible on its own and leaves
#'   your session's random-number stream as it was.
#' @return An object of class `raddiv_het`, a list of:
#'   \describe{
#'     \item{individual_heterozygosity}{One row per individual: loci called,
#'       heterozygosity (on the pooled locus set) and individual `F`.}
#'     \item{population_summary}{Per population: individuals, loci, mean, SD,
#'       SE, range of individual heterozygosity, and mean `F`.}
#'     \item{omnibus}{With 3 or more populations: one overall test of whether
#'       any population differs, for heterozygosity and for `F` -- Welch's
#'       one-way ANOVA (`welch_F`, `welch_df1`, `welch_df2`, `p_welch`) and
#'       Kruskal-Wallis (`kruskal_chisq`, `kruskal_df`, `p_kruskal`) -- on the
#'       loci that clear `min_call` in every population (`n_loci`). `NULL`
#'       for 2 populations, where the pairwise test is the only test.}
#'     \item{pairwise_tests}{Heterozygosity, one row per pair of populations:
#'       means, difference with Welch 95% CI, Welch and Wilcoxon p-values
#'       (with Benjamini-Hochberg adjusted versions), Hedges' g.}
#'     \item{pairwise_F_tests}{The same tests on `F`.}
#'     \item{overdispersion}{Observed vs binomial SD of individual
#'       heterozygosity per population, their variance ratio, and roughly how
#'       much a locus bootstrap would understate the SE.}
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
#' @examples
#' # A toy dataset shipped with the package (4 and 3 individuals -- far too
#' # few for a real test, which needs individuals, not loci).
#' vcf    <- system.file("extdata", "small.haps.vcf", package = "RADdiversity")
#' popmap <- system.file("extdata", "small_popmap.tsv", package = "RADdiversity")
#' set.seed(1)
#' res <- het_between_pops(vcf, popmap, min_call = 0.5)
#' res                       # the main tables
#' res$individual_heterozygosity
#' summary(res)              # the full report, with the confound and overdispersion checks
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

  ## ---- 3. Individuals: pooled loci, heterozygosity, exclusions -------------
  individuals <- .individual_het(H, pops, min_call, min_loci, verbose)
  pops <- individuals$pops                  # without excluded individuals

  ## ---- 4. Populations -------------------------------------------------------
  locus_sets <- .population_locus_sets(H, pops, min_call)
  per_pop <- .population_summary(H, pops, locus_sets)
  table_ind <- individuals$table
  table_ind$F <- unname(per_pop$F_by_individual[table_ind$sample])

  confound <- .missingness_confound(individuals$call_rate, individuals$het,
                                    table_ind$population, names(pops))
  overdispersion <- .overdispersion(H, pops, locus_sets, per_pop$het_by_pop)
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
                   seed = seed, files = written)
  structure(list(individual_heterozygosity = table_ind,
                 population_summary = per_pop$table,
                 omnibus = omnibus,
                 pairwise_tests = tests$heterozygosity,
                 pairwise_F_tests = tests$F,
                 overdispersion = overdispersion,
                 g2 = g2,
                 missingness_confound = confound,
                 settings = settings),
            class = "raddiv_het")
}

## ---------------------------------------------------------------------------
## Steps
## ---------------------------------------------------------------------------

## Not exported. The individual-level part of het_between_pops():
##   1. keep loci genotyped in >= min_call of ALL individuals pooled (used
##      only for this table and the missingness check; see ?het_between_pops);
##   2. heterozygosity of each individual = heterozygous loci / called loci;
##   3. exclude individuals called at fewer than `min_loci` loci. An
##      individual called nowhere would otherwise give NaN, which would
##      propagate silently through every test.
## Returns list(table, pops (without excluded individuals), call_rate, het,
## n_loci_pooled).
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
            "pooled set and may be unreliable; the population summary and the ",
            "pairwise tests use each population's own loci and are not affected.")
  a1 <- a1[pooled, , drop = FALSE]
  a2 <- a2[pooled, , drop = FALSE]

  het <- colMeans(a1 != a2, na.rm = TRUE)
  n_called <- colSums(!is.na(a1))

  too_few <- n_called < min_loci
  if (any(too_few)) {
    .inform(verbose, sprintf("  %d individual(s) called at fewer than %d loci -- EXCLUDED: %s",
                             sum(too_few), min_loci, paste(ids[too_few], collapse = ", ")))
    .inform(verbose, "    (a library this poor cannot contribute a heterozygosity ",
            "estimate; check it before re-running)")
    ids <- ids[!too_few]
    pop_of_id <- pop_of_id[!too_few]
    het <- het[!too_few]
    n_called <- n_called[!too_few]
    pops <- lapply(pops, function(v) v[v %in% ids])
    small <- names(pops)[lengths(pops) < 2]
    if (length(small))
      stop("After excluding those individuals, population(s) ",
           paste(small, collapse = ", "), " have fewer than 2 individuals.\n",
           "  A two-sample test needs at least 2 per group.", call. = FALSE)
  }
  if (!all(is.finite(het)))
    stop("Non-finite per-individual heterozygosity remains after filtering. ",
         "This should not happen; please report it.", call. = FALSE)

  table <- data.frame(sample = ids, population = pop_of_id,
                      loci_called = as.integer(n_called), heterozygosity = unname(het),
                      row.names = NULL)
  list(table = table, pops = pops, call_rate = unname(n_called) / sum(pooled),
       het = unname(het), n_loci_pooled = sum(pooled))
}

## Not exported. Per population, on its own locus set: every individual's
## heterozygosity and inbreeding F (.ind_F(), R/individual_inbreeding.R), and
## the summary table.
.population_summary <- function(H, pops, locus_sets) {
  het_by_pop <- F_by_pop <- stats::setNames(vector("list", length(pops)), names(pops))
  for (p in names(pops)) {
    loci <- locus_sets[, p]
    a1 <- H$A1[loci, pops[[p]], drop = FALSE]
    a2 <- H$A2[loci, pops[[p]], drop = FALSE]
    het_by_pop[[p]] <- colMeans(a1 != a2, na.rm = TRUE)
    F_by_pop[[p]] <- stats::setNames(.ind_F(a1, a2)$F, pops[[p]])
  }
  table <- do.call(rbind, lapply(names(pops), function(p) {
    h <- het_by_pop[[p]]
    data.frame(population = p, n = length(h), n_loci = sum(locus_sets[, p]),
               mean_het = mean(h), sd = stats::sd(h), se = stats::sd(h) / sqrt(length(h)),
               min = min(h), max = max(h), mean_F = mean(F_by_pop[[p]], na.rm = TRUE))
  }))
  list(table = table, het_by_pop = het_by_pop,
       F_by_individual = unlist(unname(F_by_pop)))
}

## Not exported. MISSINGNESS CONFOUND. An individual's heterozygosity is
## computed over the loci it was called at. If poorly sequenced individuals
## are both more often missing and (through allele dropout) more often called
## homozygous, a difference between populations can be a difference in
## library quality. Returns list(no_variation, overall, within):
##   overall  one row: correlation r of call rate with heterozygosity, its
##            p-value, n, and the gap in mean call rate between populations
##   within   the same correlation within each population of > 3 individuals
.missingness_confound <- function(call_rate, het, pop_of_individual, pop_names) {
  if (stats::sd(call_rate) < .zero_tol)
    return(list(no_variation = TRUE, overall = NULL, within = NULL))
  test <- suppressWarnings(stats::cor.test(call_rate, het))
  within <- do.call(rbind, lapply(pop_names, function(p) {
    rows <- which(pop_of_individual == p)
    if (length(rows) > 3 && stats::sd(call_rate[rows]) > .zero_tol)
      data.frame(population = p, r = suppressWarnings(stats::cor(call_rate[rows], het[rows])),
                 min_call_rate = min(call_rate[rows]), max_call_rate = max(call_rate[rows]))
  }))
  overall <- data.frame(r = unname(test$estimate), p_value = test$p.value, n = length(het),
                        call_rate_gap = diff(range(tapply(call_rate, pop_of_individual, mean))))
  list(no_variation = FALSE, overall = overall, within = within)
}

## Not exported. IS THE LOCUS BOOTSTRAP INVALID FOR *THIS* DATASET? That
## depends on whether individuals vary in heterozygosity by MORE than
## coin-flip (binomial) noise. If every individual had the same true
## heterozygosity, individual j's heterozygosity would be a mean of L
## Bernoulli draws, with variance sum_l h_l(1 - h_l) / L^2 (h_l = the locus's
## heterozygote frequency in the population). The ratio of the observed
## variance to that expectation is the overdispersion factor: 1 means no
## excess variation, and a locus bootstrap understates the standard error of a
## population mean by roughly sqrt(overdispersion).
##
## h_l is estimated from the same n_l individuals typed at the locus, and the
## plug-in h_l(1 - h_l) is then too small by a factor (n_l - 1)/n_l, which
## would make the ratio read about n/(n - 1) with no excess variation at all
## (1.21 at n = 5, 1.07 at n = 10 in simulation). Multiplying each locus's
## term by n_l/(n_l - 1) removes that; loci typed in fewer than 2 individuals
## are left out.
.overdispersion <- function(H, pops, locus_sets, het_by_pop) {
  do.call(rbind, lapply(names(pops), function(p) {
    loci <- locus_sets[, p]
    a1 <- H$A1[loci, pops[[p]], drop = FALSE]
    a2 <- H$A2[loci, pops[[p]], drop = FALSE]
    n_typed <- rowSums(!is.na(a1))
    locus_het <- rowSums(a1 != a2, na.rm = TRUE) / n_typed
    keep <- n_typed >= 2
    n_typed <- n_typed[keep]
    locus_het <- locus_het[keep]
    expected_var <- sum(locus_het * (1 - locus_het) * n_typed / (n_typed - 1)) /
      length(locus_het)^2
    observed_var <- stats::var(het_by_pop[[p]])
    ## Monomorphic at every locus: 0/0, which is "undefined", not "no excess".
    ratio <- if (expected_var > 0) observed_var / expected_var else NA_real_
    data.frame(population = p, n = length(pops[[p]]),
               obs_sd = sqrt(observed_var), binomial_sd = sqrt(expected_var),
               overdispersion = ratio, se_understated_by = sqrt(ratio))
  }))
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

## Not exported. Tests every pair of populations, on the loci where BOTH
## populations of the pair clear min_call on their own (not the pooled set,
## and not "every population in the run"). Each row carries its own n_loci,
## n1, mean1, n2, mean2. Welch's t is primary, Wilcoxon the distribution-free
## check; with k populations there are k(k-1)/2 pairs, so p-values are also
## Benjamini-Hochberg adjusted across pairs.
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
    a1_1 <- H$A1[loci, pops[[p1]], drop = FALSE]
    a2_1 <- H$A2[loci, pops[[p1]], drop = FALSE]
    a1_2 <- H$A1[loci, pops[[p2]], drop = FALSE]
    a2_2 <- H$A2[loci, pops[[p2]], drop = FALSE]
    het_1 <- colMeans(a1_1 != a2_1, na.rm = TRUE)
    het_2 <- colMeans(a1_2 != a2_2, na.rm = TRUE)
    ## Individual F on the same loci and individuals, with expected
    ## heterozygosity from each population's own allele frequencies.
    F_1 <- .ind_F(a1_1, a2_1)$F[is.finite(het_1)]
    F_2 <- .ind_F(a1_2, a2_2)$F[is.finite(het_2)]
    list(het = .two_sample(het_1, het_2, p1, p2, n_loci, "heterozygosity", verbose),
         F = .two_sample(F_1, F_2, p1, p2, n_loci, "F", verbose))
  })
  bh_adjust <- function(tab) {
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
             n2 = n2, mean2 = mean2, diff = NA_real_, ci_lo = NA_real_,
             ci_hi = NA_real_, p_welch = NA_real_, p_wilcox = NA_real_,
             hedges_g = NA_real_)
}

## Not exported. Welch's t (primary), Wilcoxon (distribution-free check) and
## Hedges' g for one pair of populations, on one number per individual
## (`a`, `b`: heterozygosity or F). Returns one row of the pairwise table.
.two_sample <- function(a, b, p1, p2, n_loci, what, verbose) {
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

  data.frame(pop1 = p1, pop2 = p2, n_loci = n_loci,
             n1 = n_a, mean1 = mean(a), n2 = n_b, mean2 = mean(b),
             diff = mean(a) - mean(b), ci_lo = welch$conf.int[1], ci_hi = welch$conf.int[2],
             p_welch = welch$p.value, p_wilcox = wilcox$p.value,
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
                        signif_cols = c(p_welch = 3, p_wilcox = 3, p_welch_BH = 3, p_wilcox_BH = 3)),
  overdispersion = list(digits = 4, round_cols = c(obs_sd = 5, binomial_sd = 5,
                                                   overdispersion = 1, se_understated_by = 1)),
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
    .print_table(.round_het_table(tab[utils::head(order(tab$p_welch_BH), max_rows), ],
                                  "pairwise_tests"))
  }
}

#' @rdname het_between_pops
#' @param x,object A `raddiv_het` object, as returned by `het_between_pops()`.
#' @param ... Ignored.
#' @export
print.raddiv_het <- function(x, ...) {
  st <- x$settings
  cat(sprintf("Heterozygosity between populations, individuals as replicates: %d populations, %d individuals\n",
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
  cat("Also in the result: $individual_heterozygosity, $overdispersion, $g2, $missingness_confound\n")
  cat("summary() prints the full report, with the confound and overdispersion checks explained.\n")
  invisible(x)
}

## Not exported. The short warnings print.raddiv_het() shows.
.het_warnings <- function(x) {
  out <- character(0)
  cf <- x$missingness_confound
  if (!cf$no_variation && is.finite(cf$overall$p_value) && cf$overall$p_value < 0.05 &&
      cf$overall$r > 0)
    out <- c(out, "call rate is positively correlated with heterozygosity (possible allele dropout).")
  if (!cf$no_variation && cf$overall$call_rate_gap > 0.02)
    out <- c(out, "mean call rate differs between populations by more than 2 percentage points.")
  od <- x$overdispersion$overdispersion
  if (any(is.finite(od)) && max(od, na.rm = TRUE) >= 2)
    out <- c(out, "individuals vary more than binomial noise; locus-based intervals would be too narrow.")
  if (any(x$g2$g2_lo > 0, na.rm = TRUE))
    out <- c(out, "g2 > 0 in some population: report diversity_stats(se_individuals = TRUE).")
  out
}

#' @rdname het_between_pops
#' @export
summary.raddiv_het <- function(object, ...) {
  structure(list(result = object), class = "summary.raddiv_het")
}

#' @export
print.summary.raddiv_het <- function(x, ...) {
  .het_report(x$result)
  invisible(x)
}

## Not exported. The full heterozygosity report (summary() of a result, and
## what the command-line script prints).
.het_report <- function(x) {
  st <- x$settings
  note <- function(...) cat(paste0("  ", c(...), "\n"), sep = "")
  tests <- x$pairwise_tests
  rule <- strrep("=", 69)

  cat("\n", rule, "\n", sep = "")
  cat("  HETEROZYGOSITY BETWEEN POPULATIONS -- individual as the unit\n")
  cat("  ", .big(st$n_loci_pooled), " loci pass the pooled filter (used for the\n",
      "  individual-level table and confound check below); ", st$n_individuals,
      " individuals, ", st$n_pops, " populations\n", sep = "")
  cat(rule, "\n", sep = "")

  cat("\nPer-population summary of individual heterozygosity\n")
  cat("  (each population's OWN locus set -- loci at >= ", 100 * st$min_call,
      "% call rate WITHIN that population; see 'n_loci'. This can differ from\n",
      "   a specific pair's test below, which uses only loci BOTH members of\n",
      "   that pair clear -- that is expected, not an inconsistency.)\n", sep = "")
  .print_table(.round_het_table(x$population_summary, "population_summary"))

  cat("\nMissingness confound check\n")
  cf <- x$missingness_confound
  if (cf$no_variation) {
    note("Every individual was called at every locus; no confound possible.")
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

  cat("\nOverdispersion check -- does the locus bootstrap fail on THIS dataset?\n")
  note("(each population's own locus set, as in the summary above)")
  od <- x$overdispersion
  .print_table(.round_het_table(od, "overdispersion"))
  cat("\n")
  if (all(is.na(od$overdispersion))) {
    note(.report_text$het_overdispersion_undefined)
  } else {
    worst <- max(od$overdispersion, na.rm = TRUE)
    note(if (worst < 2) .report_text$het_overdispersion_small
         else if (worst < 5) .report_text$het_overdispersion_moderate
         else .report_text$het_overdispersion_large)
  }
  note(.report_text$het_van_dongen)

  cat("\nIdentity disequilibrium g2 (David et al. 2007) -- the standard measure of\n")
  cat("variance in inbreeding among individuals; 0 when individuals do not differ.\n")
  .print_table(.round_het_table(x$g2, "g2"))
  note(.report_text$het_g2)

  if (!is.null(x$omnibus)) {
    cat(sprintf("\nOverall test across all %d populations (loci that clear min_call in every one)\n",
                st$n_pops))
    .print_table(.round_het_table(x$omnibus, "omnibus"))
    note(.report_text$het_omnibus)
  }

  n_pairs <- nrow(tests)
  cat(sprintf("\nPairwise tests -- %d comparison%s. Welch's t is the primary test;\n",
              n_pairs, if (n_pairs == 1) "" else "s"))
  cat("Wilcoxon is the distribution-free check. 'diff' and its CI are from Welch.\n")
  if (n_pairs > 1)
    cat(sprintf("p_*_BH are Benjamini-Hochberg adjusted across all %d comparisons.\n", n_pairs))
  .show_pairwise(tests, 25, "x$pairwise_tests")
  n_significant <- function(tab) {
    c(welch = sum(tab$p_welch_BH < 0.05, na.rm = TRUE),
      wilcox = sum(tab$p_wilcox_BH < 0.05, na.rm = TRUE))
  }
  sig <- n_significant(tests)
  cat(sprintf("\n  pairs significant at BH < 0.05: Welch %d, Wilcoxon %d, of %d\n",
              sig[["welch"]], sig[["wilcox"]], n_pairs))
  note("Hedges' g is the standardised difference: ~0.2 small, ~0.5 medium, ~0.8 large.")

  F_tests <- x$pairwise_F_tests
  cat("\n")
  cat(paste0(.report_text$het_F_tests, "\n"), sep = "")
  .show_pairwise(F_tests, 25, "x$pairwise_F_tests")
  sig_F <- n_significant(F_tests)
  cat(sprintf("\n  pairs significant at BH < 0.05: Welch %d, Wilcoxon %d, of %d\n",
              sig_F[["welch"]], sig_F[["wilcox"]], nrow(F_tests)))

  cat("\nInterpretation notes\n")
  note(.report_text$het_interpretation)
  if (length(st$files)) cat("\nWrote ", paste(st$files, collapse = " and "), "\n", sep = "")
  cat("\n")
  invisible(NULL)
}
