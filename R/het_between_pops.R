###############################################################################
#
#  R/het_between_pops.R -- does mean heterozygosity differ between populations?
#  Each individual contributes ONE number, so individuals are the replicate.
#
#  WHY NOT THE USUAL TEST. Bootstrapping over loci, or a paired Wilcoxon on
#  per-locus values, treats LOCI as the replicate. Loci are repeated measures on
#  the same individuals, so the interval shrinks as 1/sqrt(n_loci) while the
#  real uncertainty -- which individuals you caught -- does not shrink at all.
#  Measured type I error against a nominal 5%, two populations with identical
#  true heterozygosity: 5% when individuals are identical, but 55% and 76% as
#  soon as they differ, which is exactly what inbreeding produces. Welch's t on
#  per-individual heterozygosity held 5% throughout. het_between_pops_selftest()
#  reproduces that table. The objection is old (Van Dongen 1995); most papers
#  ignore it.
#  Full discussion in README.md, "Step 3".
#
#  WHAT TO READ, IN ORDER
#    1. the missingness confound check -- if call rate correlates with
#       heterozygosity you are measuring library quality, not biology
#    2. the overdispersion factor -- how badly a locus bootstrap would have
#       misled you ON THESE DATA; quote it if a reviewer asks
#    3. the test itself, Welch t with Wilcoxon as a distribution-free check
#
#  Power is set by the number of INDIVIDUALS, not loci: ~52% for a moderate
#  difference at n = 15 vs 10. Report "no difference was detected", not "there
#  is no difference". Relatives are not independent units; check first.
#
#  Needs the VCF: per-individual heterozygosity cannot be recovered from
#  sumstats.tsv, which holds only per-population counts.
#
###############################################################################

#' Test whether mean heterozygosity differs between populations
#'
#' Uses the INDIVIDUAL, not the locus, as the unit of replication (Welch's t
#' with Wilcoxon as a distribution-free check), plus a missingness-confound
#' diagnostic and an overdispersion check quantifying how badly a locus
#' bootstrap would have misled on this dataset. Prints the same report as the
#' `het_between_pops.R` command-line script, writes the same two TSV files,
#' and additionally returns the result tables as data frames.
#'
#' @param vcf_file Path to a Stacks VCF (`populations.snps.vcf` or
#'   `populations.haps.vcf`, optionally gzip-compressed).
#' @param popmap_f Path to a two-column, no-header popmap TSV (`sample_id
#'   <TAB> population`).
#' @param min_call Minimum per-individual genotyping rate for a locus to be
#'   used. Default `0.9`.
#' @param outdir Directory to write the two output TSVs into. Default `"."`
#'   (the current directory), created if it does not exist.
#' @param seed Random seed set before any randomness, for reproducibility.
#'   Default `2024`. The caller's own RNG state (e.g. from their own prior
#'   `set.seed()`) is restored when this function returns, so calling it
#'   interactively does not affect subsequent random draws in the caller's
#'   session.
#' @param verbose Print the VCF/popmap parsing summary. Default `TRUE`. The
#'   main report is always printed regardless of this setting, matching the
#'   command-line script.
#' @return Invisibly, a list with elements `individual_heterozygosity` and
#'   `pairwise_tests` (the same two data frames written to
#'   `individual_heterozygosity.tsv` and `het_between_pops_tests.tsv`).
#' @export
het_between_pops <- function(vcf_file, popmap_f, min_call = 0.9, outdir = ".",
                              seed = 2024, verbose = TRUE) {

  if (!file.exists(vcf_file))
    stop("VCF file not found: ", vcf_file, "\n  Check the path and try again.")
  if (!file.exists(popmap_f))
    stop("Popmap file not found: ", popmap_f, "\n  Check the path and try again.")
  min_call <- as.numeric(min_call)
  if (is.na(min_call) || min_call < 0 || min_call > 1) stop("min_call must be in 0-1.")
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  ## Restore the caller's RNG state on exit -- see the matching comment in
  ## diversity_stats().
  restore_rng <- .save_rng_state()
  on.exit(restore_rng(), add = TRUE)
  set.seed(seed)

  message("Reading ", vcf_file, " ...")
  H    <- read_haps_vcf(vcf_file, verbose = verbose)
  pops <- read_popmap(popmap_f, H$samples, verbose = verbose)
  r <- length(pops)
  if (r < 2) stop("Need at least 2 populations.")

  ids  <- unlist(pops, use.names = FALSE)
  plab <- rep(names(pops), lengths(pops))
  A1 <- H$A1[, ids, drop = FALSE]; A2 <- H$A2[, ids, drop = FALSE]

  ## keep loci genotyped in enough individuals, POOLED across all populations.
  ## This pooled set is used ONLY for individual-level diagnostics below (the
  ## raw per-individual table and the missingness confound check) -- not for
  ## anything that compares populations to each other. A pooled threshold lets
  ## a locus pass because the larger population's coverage carries the smaller
  ## one's (see README.md's warning against -R for this same reason),
  ## which is exactly the asymmetry a BETWEEN-population estimate must avoid.
  ## The population summary, overdispersion check, and pairwise tests further
  ## down use a genuinely per-population/per-pair criterion instead (see below).
  cr <- rowMeans(!is.na(A1))
  keep <- cr >= min_call - 1e-9
  message(sprintf("Loci genotyped in >= %.0f%% of individuals (pooled): %s of %s",
                  100 * min_call, format(sum(keep), big.mark = ","),
                  format(nrow(A1), big.mark = ",")))
  if (sum(keep) < 50)
    message("  WARNING: fewer than 50 loci pass the pooled filter. The raw ",
            "per-individual table and the missingness confound check below ",
            "use this pooled set and may be unreliable as a result; the ",
            "per-population summary and the pairwise tests further down use a ",
            "different, per-population basis and are not affected by this.")
  A1 <- A1[keep, , drop = FALSE]; A2 <- A2[keep, , drop = FALSE]
  L <- nrow(A1)

  ## ---------------------------------------------------------------------------
  ## Per-individual heterozygosity: the proportion of that individual's CALLED
  ## loci at which it is heterozygous. Each individual contributes ONE number.
  ## ---------------------------------------------------------------------------
  het   <- colMeans(A1 != A2, na.rm = TRUE)
  ncall <- colSums(!is.na(A1))

  ## An individual called at NO loci gives colMeans(na.rm = TRUE) = NaN, which
  ## then propagates silently: mean_het reads NaN, every test reads NA, and the
  ## function would still return normally. A failed library is exactly the case
  ## this function is meant to detect, so it must not be the case it fails
  ## quietly on. Drop such individuals loudly and refuse if that leaves a
  ## population too small.
  .thin <- ncall < 50
  if (any(.thin)) {
    message(sprintf("  %d individual(s) called at fewer than 50 loci -- EXCLUDED: %s",
                    sum(.thin), paste(ids[.thin], collapse = ", ")))
    message("    (a library this poor cannot contribute a heterozygosity estimate;",
            " check it before re-running)")
    ids <- ids[!.thin]; plab <- plab[!.thin]
    A1 <- A1[, !.thin, drop = FALSE]; A2 <- A2[, !.thin, drop = FALSE]
    het <- het[!.thin]; ncall <- ncall[!.thin]
    pops <- lapply(pops, function(v) v[v %in% ids])
    small <- names(pops)[lengths(pops) < 2]
    if (length(small))
      stop("After excluding those individuals, population(s) ",
           paste(small, collapse = ", "), " have fewer than 2 individuals.\n",
           "  A two-sample test needs at least 2 per group.")
    pops <- pops[lengths(pops) >= 2]
    r <- length(pops)
    if (r < 2) stop("Fewer than 2 usable populations remain.")
  }
  if (!all(is.finite(het)))
    stop("Non-finite per-individual heterozygosity remains after filtering. ",
         "This should not happen; please report it.")
  ind <- data.frame(sample = ids, population = plab,
                    loci_called = ncall, heterozygosity = round(het, 5),
                    row.names = NULL)

  ## ---------------------------------------------------------------------------
  ## Per-locus, per-population call rate, and per-population locus sets.
  ##
  ## The pooled `keep` set above answers "well genotyped across everyone" and is
  ## right for individual-level diagnostics (see the comment where it's built).
  ## Everything below that describes or compares POPULATIONS uses a genuinely
  ## per-population criterion instead: a locus belongs to population p's own
  ## set iff population p ITSELF clears min_call there, independent of how any
  ## other population happens to be covered. This is computed from the full,
  ## unfiltered locus set (H$A1/H$A2), not from the pooled `keep` set, because
  ## the pooled filter can drop a locus that is perfectly well covered in one
  ## population just because some unrelated population drags the pooled rate
  ## down -- the asymmetry this whole redesign exists to avoid.
  ## ---------------------------------------------------------------------------
  n_rec_all <- nrow(H$A1)
  cr_pop <- matrix(NA_real_, n_rec_all, r, dimnames = list(NULL, names(pops)))
  for (p in names(pops))
    cr_pop[, p] <- rowMeans(!is.na(H$A1[, pops[[p]], drop = FALSE]))
  keep_pop <- cr_pop >= min_call - 1e-9

  ## Per-population individual heterozygosity, each on ITS OWN locus set.
  pop_het <- stats::setNames(vector("list", r), names(pops))
  for (p in names(pops)) {
    lp <- keep_pop[, p]
    a1p <- H$A1[lp, pops[[p]], drop = FALSE]; a2p <- H$A2[lp, pops[[p]], drop = FALSE]
    pop_het[[p]] <- colMeans(a1p != a2p, na.rm = TRUE)
  }

  cat("\n=====================================================================\n")
  cat("  HETEROZYGOSITY BETWEEN POPULATIONS -- individual as the unit\n")
  cat("  ", format(L, big.mark = ","), " loci pass the pooled filter (used for the\n",
      "  individual-level table and confound check below); ", length(ids),
      " individuals, ", r, " populations\n", sep = "")
  cat("=====================================================================\n")

  summ <- do.call(rbind, lapply(names(pops), function(pn) {
    v <- pop_het[[pn]]
    data.frame(population = pn, n = length(v), n_loci = sum(keep_pop[, pn]),
               mean_het = round(mean(v), 4), sd = round(stats::sd(v), 4),
               se = round(stats::sd(v) / sqrt(length(v)), 4),
               min = round(min(v), 4), max = round(max(v), 4)) }))
  cat("\nPer-population summary of individual heterozygosity\n")
  cat("  (each population's OWN locus set -- loci at >= ", 100 * min_call,
      "% call rate WITHIN that population; see 'n_loci'. This can differ from\n",
      "   a specific pair's test below, which uses only loci BOTH members of\n",
      "   that pair clear -- that is expected, not an inconsistency.)\n", sep = "")
  print(summ, row.names = FALSE)
  ## ---------------------------------------------------------------------------
  ## MISSINGNESS CONFOUND. An individual's heterozygosity is computed over the
  ## loci it was CALLED at. If poorly sequenced individuals are both more missing
  ## and (through allele dropout) more apparently homozygous, then a difference
  ## in heterozygosity between populations can be a difference in library
  ## quality. This is the same mechanism the rest of the package exists to catch,
  ## and it operates on the individual axis here rather than the site axis.
  cat("\nMissingness confound check\n")
  cr_ind <- ncall / L
  if (stats::sd(cr_ind) < 1e-12) {
    cat("  Every individual was called at every locus; no confound possible.\n")
  } else {
    ct <- suppressWarnings(stats::cor.test(cr_ind, het))
    cat(sprintf("  cor(per-individual call rate, heterozygosity) = %+.3f  (p = %.3g, n = %d)\n",
                unname(ct$estimate), ct$p.value, length(het)))
    ## and per population, so a between-population difference in call rate is
    ## not read as a within-population correlation and vice versa
    for (pn in names(pops)) {
      jj <- which(plab == pn)
      if (length(jj) > 3 && stats::sd(cr_ind[jj]) > 1e-12)
        cat(sprintf("    within %-12s r = %+.3f   call rate %.3f-%.3f\n", pn,
                    suppressWarnings(stats::cor(cr_ind[jj], het[jj])),
                    min(cr_ind[jj]), max(cr_ind[jj])))
    }
    if (is.finite(ct$p.value) && ct$p.value < 0.05 && ct$estimate > 0) {
      cat("  POSITIVE and significant: individuals with more missing data look LESS\n")
      cat("  heterozygous, which is the allele-dropout signature. Re-run on a\n")
      cat("  complete-data locus set (min_call = 1.0) before quoting the test\n")
      cat("  below; if the difference survives that, it is not a coverage artefact.\n")
    } else {
      cat("  No positive association, so heterozygosity is not tracking coverage.\n")
    }
    if (diff(range(tapply(cr_ind, plab, mean))) > 0.02)
      cat("  WARNING: mean call rate differs between populations by more than 2\n          percentage points. That asymmetry is itself a candidate explanation.\n")
  }

  cat("\n  'sd' is the spread AMONG INDIVIDUALS. That spread, divided by sqrt(n), is\n")
  cat("  the real uncertainty in a population mean -- and it does not shrink at all\n")
  cat("  as you add loci, which is why locus-based tests can be invalid here.\n")

  ## ---------------------------------------------------------------------------
  ## IS THE LOCUS BOOTSTRAP ACTUALLY INVALID FOR *THIS* DATASET?
  ##
  ## It depends on one empirical question: do individuals vary in heterozygosity
  ## by MORE than coin-flip (binomial) noise?
  ##
  ## If every individual had the same true heterozygosity, the spread among
  ## individuals would be exactly the binomial expectation, and locus-based
  ## resampling would be correctly calibrated (measured type I error 5.0%). Any
  ## EXCESS spread is variation that averaging over loci cannot remove, and it is
  ## what breaks locus-based tests.
  ##
  ## Under the null of no individual variation, individual j's heterozygosity is
  ## a mean of L Bernoulli draws, so Var = sum over loci of h(1-h) / L^2, where h
  ## is the locus's heterozygote frequency in that population. The ratio of the
  ## observed variance to that expectation is an overdispersion factor: 1 means no
  ## excess individual variation, and the standard error of a population mean is
  ## understated by roughly sqrt(overdispersion) if you resample loci instead.
  ## ---------------------------------------------------------------------------
  cat("\nOverdispersion check -- does the locus bootstrap fail on THIS dataset?\n")
  cat("  (each population's own locus set, as in the summary above)\n")
  od <- do.call(rbind, lapply(names(pops), function(pn) {
    lp <- keep_pop[, pn]; nj <- length(pops[[pn]])
    a1 <- H$A1[lp, pops[[pn]], drop = FALSE]; a2 <- H$A2[lp, pops[[pn]], drop = FALSE]
    hl <- rowMeans(a1 != a2, na.rm = TRUE)          # per-locus het freq in this pop
    hl <- hl[is.finite(hl)]
    exp_var <- sum(hl * (1 - hl)) / length(hl)^2    # binomial expectation
    obs_var <- stats::var(pop_het[[pn]])
    ## A population monomorphic at every locus used gives exp_var = 0 (and
    ## usually obs_var = 0 too, since every individual is 0% heterozygous),
    ## so the ratio is 0/0 -- undefined, not "no excess variation". Report NA
    ## rather than the literal text "NaN"/"NaNx".
    od_ratio <- if (exp_var > 0) obs_var / exp_var else NA_real_
    data.frame(population = pn, n = nj,
               obs_sd = round(sqrt(obs_var), 5),
               binomial_sd = round(sqrt(exp_var), 5),
               overdispersion = round(od_ratio, 1),
               SE_understated_by = if (is.na(od_ratio)) "n/a (no variation)"
                                    else paste0(round(sqrt(od_ratio), 1), "x"))
  }))
  print(od, row.names = FALSE)
  if (all(is.na(od$overdispersion))) {
    cat("\n  Overdispersion is undefined for every population (no variation in either\n")
    cat("  heterozygosity or per-locus heterozygote frequency) -- too degenerate a\n")
    cat("  dataset to say whether a locus bootstrap would fail here.\n")
  } else {
    worst <- max(od$overdispersion, na.rm = TRUE)
    cat("\n")
    if (worst < 2) {
      cat("  Overdispersion is near 1: individuals differ about as much as coin-flip\n")
      cat("  noise alone predicts. For this dataset a locus bootstrap would be roughly\n")
      cat("  correctly calibrated, and the two approaches should broadly agree. The\n")
      cat("  individual-level test below remains the safer default.\n")
    } else if (worst < 5) {
      cat("  Overdispersion is moderate. Individuals vary by more than coin-flip noise,\n")
      cat("  so a locus bootstrap understates the standard error by roughly the factor\n")
      cat("  in the last column and its p-values are too small. Use the test below.\n")
    } else {
      cat("  Overdispersion is LARGE. Individuals differ far more than coin-flip noise,\n")
      cat("  so a locus bootstrap would understate the standard error by the factor in\n")
      cat("  the last column and reject far too often. Use the test below, and inspect\n")
      cat("  the per-individual table for outliers (a bad library looks like one\n")
      cat("  individual with unusually LOW heterozygosity).\n")
    }
  }
  cat("  This is the point made by Van Dongen (1995, Heredity 74:445-447): the unit\n")
  cat("  of resampling changes what the bootstrap means, and loci are usually the\n")
  cat("  wrong unit because they are all measured on the same individuals.\n")

  ## ---------------------------------------------------------------------------
  ## Pairwise tests. Each pair uses ONLY the loci where BOTH members of that
  ## pair individually clear min_call -- the intersection of their own
  ## per-population locus sets above (S_A n S_B). Deliberately NOT the pooled
  ## set (would let one population's coverage carry the other's, the bug this
  ## redesign fixes) and NOT "every population in the run must pass" (would
  ## penalize every pair for whichever population anywhere in the dataset
  ## happens to be worst covered -- needlessly strict once r > 2). Each row
  ## carries its own n_loci/n1/mean1/n2/mean2, so the table is self-contained
  ## and does not need to be cross-checked against the summary above, which is
  ## computed on a different (per-population, not per-pair) basis by design.
  ## Welch's t is primary; Wilcoxon is the distribution-free check. With k
  ## populations there are k(k-1)/2 pairs, so p-values are BH-adjusted across
  ## all of them.
  ## ---------------------------------------------------------------------------
  rows <- list()
  for (i in seq_len(r - 1)) for (j in (i + 1):r) {
    pi <- names(pops)[i]; pj <- names(pops)[j]
    lij <- keep_pop[, pi] & keep_pop[, pj]
    nloc <- sum(lij)
    if (nloc < 50) {
      message(sprintf(
        "  %s vs %s: only %d loci clear %.0f%% call rate in BOTH populations ",
        pi, pj, nloc, 100 * min_call),
        "(need >= 50). Reported as NA for this pair; other pairs are unaffected.")
      rows[[length(rows) + 1]] <- data.frame(
        pop1 = pi, pop2 = pj, n_loci = nloc,
        n1 = length(pops[[pi]]), mean1 = NA_real_,
        n2 = length(pops[[pj]]), mean2 = NA_real_,
        diff = NA_real_, ci_lo = NA_real_, ci_hi = NA_real_,
        p_welch = NA_real_, p_wilcox = NA_real_, hedges_g = NA_real_)
      next
    }
    a1i <- H$A1[lij, pops[[pi]], drop = FALSE]; a2i <- H$A2[lij, pops[[pi]], drop = FALSE]
    a1j <- H$A1[lij, pops[[pj]], drop = FALSE]; a2j <- H$A2[lij, pops[[pj]], drop = FALSE]
    a <- colMeans(a1i != a2i, na.rm = TRUE)
    b <- colMeans(a1j != a2j, na.rm = TRUE)
    ## An individual can pass the earlier GLOBAL thinning check (>= 50 calls on
    ## the pooled locus set) yet still be called at NONE of THIS pair's smaller,
    ## pair-specific locus set, giving colMeans(na.rm = TRUE) = NaN here. Drop
    ## such individuals from this comparison only -- same reasoning as the
    ## global thinning step above, just re-applied on the pair's own basis.
    if (any(!is.finite(a)) || any(!is.finite(b))) {
      n_bad <- sum(!is.finite(a)) + sum(!is.finite(b))
      message(sprintf(
        "  %s vs %s: %d individual(s) called at none of the %d shared loci for ",
        pi, pj, n_bad, nloc), "this pair; excluded from this comparison only.")
      a <- a[is.finite(a)]; b <- b[is.finite(b)]
    }
    if (length(a) < 2 || length(b) < 2) {
      message(sprintf(
        "  %s vs %s: fewer than 2 individuals with a valid value in one group ",
        pi, pj), "after excluding those above. Reported as NA for this pair.")
      rows[[length(rows) + 1]] <- data.frame(
        pop1 = pi, pop2 = pj, n_loci = nloc,
        n1 = length(a), mean1 = if (length(a)) round(mean(a), 4) else NA_real_,
        n2 = length(b), mean2 = if (length(b)) round(mean(b), 4) else NA_real_,
        diff = NA_real_, ci_lo = NA_real_, ci_hi = NA_real_,
        p_welch = NA_real_, p_wilcox = NA_real_, hedges_g = NA_real_)
      next
    }
    ## Two distinct degenerate cases, both meaning "no variation among
    ## individuals in one or both populations, so no test is possible":
    ##   (a) t.test()/wilcox.test() ERROR outright -- happens when both samples
    ##       are constant at the SAME nonzero value (their pooled variance is 0
    ##       and R's zero-stderr guard, a RELATIVE check against that shared
    ##       value, correctly fires).
    ##   (b) they return SUCCESSFULLY but with p.value/conf.int = NaN -- happens
    ##       when both samples are constant at exactly ZERO (e.g. a population
    ##       monomorphic at every retained locus). t.test.default's zero-stderr
    ##       guard is `stderr < 10*eps*max(abs(mx), abs(my))`; with both means at
    ##       0 the right-hand side is also 0, so the guard never fires and R
    ##       hands back a "valid" object full of NaN instead of erroring.
    ##       wilcox.test() has the same blind spot. Catch both explicitly rather
    ##       than assuming an error is the only way this fails.
    tt <- try(stats::t.test(a, b), silent = TRUE)
    if (inherits(tt, "try-error")) {
      message(sprintf(
        "  %s vs %s: t.test() failed (%s). This usually means no variation ",
        pi, pj, trimws(conditionMessage(attr(tt, "condition")))),
        "among individuals in one or both populations, so no test is ",
        "possible. Reported as NA.")
      tt <- list(p.value = NA_real_, conf.int = c(NA_real_, NA_real_))
    } else if (is.nan(tt$p.value)) {
      message(sprintf(
        "  %s vs %s: t.test() returned NaN (both samples constant at exactly ",
        pi, pj),
        "zero -- e.g. monomorphic at every retained locus in both populations). ",
        "No variation to test. Reported as NA.")
      tt <- list(p.value = NA_real_, conf.int = c(NA_real_, NA_real_))
    }
    wt <- suppressWarnings(try(stats::wilcox.test(a, b), silent = TRUE))
    if (inherits(wt, "try-error") || is.nan(wt$p.value)) wt <- list(p.value = NA_real_)
    ## Hedges' g, a sample-size-corrected standardised difference
    s <- sqrt(((length(a)-1)*stats::var(a) + (length(b)-1)*stats::var(b)) / (length(a)+length(b)-2))
    d <- if (is.finite(s) && s > 0) (mean(a) - mean(b)) / s else NA_real_
    J <- 1 - 3 / (4 * (length(a) + length(b)) - 9)
    rows[[length(rows) + 1]] <- data.frame(
      pop1 = pi, pop2 = pj, n_loci = nloc,
      n1 = length(a), mean1 = round(mean(a), 4),
      n2 = length(b), mean2 = round(mean(b), 4),
      diff = round(mean(a) - mean(b), 4),
      ci_lo = round(tt$conf.int[1], 4), ci_hi = round(tt$conf.int[2], 4),
      p_welch = tt$p.value, p_wilcox = wt$p.value,
      hedges_g = round(d * J, 2))
  }
  tab <- do.call(rbind, rows)
  tab$p_welch_BH  <- stats::p.adjust(tab$p_welch,  method = "BH")
  tab$p_wilcox_BH <- stats::p.adjust(tab$p_wilcox, method = "BH")
  out <- tab
  out$p_welch <- signif(out$p_welch, 3); out$p_wilcox <- signif(out$p_wilcox, 3)
  out$p_welch_BH <- signif(out$p_welch_BH, 3); out$p_wilcox_BH <- signif(out$p_wilcox_BH, 3)

  npair <- nrow(tab)
  cat(sprintf("\nPairwise tests -- %d comparison%s. Welch's t is the primary test;\n",
              npair, if (npair == 1) "" else "s"))
  cat("Wilcoxon is the distribution-free check. 'diff' and its CI are from Welch.\n")
  if (npair > 1)
    cat(sprintf("p_*_BH are Benjamini-Hochberg adjusted across all %d comparisons.\n", npair))
  if (npair <= 25) print(out, row.names = FALSE) else {
    o <- order(out$p_welch_BH)
    cat("  (25 most significant pairs; full table written to file)\n")
    print(out[utils::head(o, 25), ], row.names = FALSE)
  }
  cat(sprintf("\n  pairs significant at BH < 0.05: Welch %d, Wilcoxon %d, of %d\n",
              sum(tab$p_welch_BH < 0.05, na.rm = TRUE),
              sum(tab$p_wilcox_BH < 0.05, na.rm = TRUE), npair))
  cat("  Hedges' g is the standardised difference: ~0.2 small, ~0.5 medium, ~0.8 large.\n")

  cat("\nInterpretation notes\n")
  cat("  * Power is limited by the NUMBER OF INDIVIDUALS, not the number of loci.\n")
  cat("    Simulation at n = 15 vs 10 with realistic individual variation gave 52%\n")
  cat("    power for a moderate difference, so a non-significant result here is\n")
  cat("    weak evidence of no difference, not evidence of no difference.\n")
  cat("  * If some individuals are relatives, they are not independent units and\n")
  cat("    even this test is anti-conservative. Check relatedness first.\n")
  cat("  * A single poor-quality library shows up as one individual with unusually\n")
  cat("    LOW heterozygosity. Inspect the per-individual table before trusting a\n")
  cat("    result, especially at small n where one outlier can drive it.\n")

  utils::write.table(ind, file.path(outdir, "individual_heterozygosity.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  utils::write.table(out, file.path(outdir, "het_between_pops_tests.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  cat("\nWrote individual_heterozygosity.tsv and het_between_pops_tests.tsv to ",
      outdir, "\n\n", sep = "")

  invisible(list(individual_heterozygosity = ind, pairwise_tests = out))
}
