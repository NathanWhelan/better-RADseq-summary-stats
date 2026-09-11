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
#' bootstrap would have misled on this dataset. Returns the result tables;
#' printing the result shows the same report as the `het_between_pops.R`
#' command-line script, and `outdir` writes the tables to TSV files.
#'
#' @param vcf_file Path to a Stacks VCF (`populations.snps.vcf` or
#'   `populations.haps.vcf`, optionally gzip-compressed) -- OR, an already-
#'   parsed (and optionally filtered) `H` list, i.e. the object returned by
#'   [read_stacks_vcf()], on its own or passed through one or more `filter_*()`
#'   functions first (see `R/filter_loci.R`). Passing a list skips re-reading
#'   and re-parsing the VCF file.
#' @param popmap_f Path to a two-column, no-header popmap TSV (`sample_id
#'   <TAB> population`).
#' @param min_call Minimum fraction of individuals that must be genotyped at
#'   a locus for that locus to be used, applied within each population.
#'   Default `0.9`.
#' @param outdir Directory to write the two result tables into as TSV files
#'   (created if needed). Default `NULL`: write no files -- the tables are in
#'   the returned object. The command-line script writes to the current
#'   directory.
#' @param seed Random seed set before any randomness, for reproducibility.
#'   Default `2024`. The caller's own RNG state (e.g. from their own prior
#'   `set.seed()`) is restored when this function returns, so calling it
#'   interactively does not affect subsequent random draws in the caller's
#'   session.
#' @param verbose Print the VCF/popmap parsing summary. Default `TRUE`.
#'   Progress messages can be silenced with `suppressMessages()`; the report
#'   itself appears only when the result is printed.
#' @param stem Text used in the two output filenames
#'   (`individual_heterozygosity.<stem>.tsv` and
#'   `het_between_pops_tests.<stem>.tsv`). Derived from `vcf_file`'s name when
#'   it is a path (`populations.snps.vcf` gives `"snps"`), so a run on the SNP
#'   VCF and a run on the haplotype VCF into the same `outdir` do not
#'   overwrite each other. Required when `vcf_file` is an already-parsed
#'   list, which has no filename to derive it from.
#' @details
#' Two questions, two tests, same machinery. `pairwise_tests` compares mean
#' individual heterozygosity: do the populations differ in diversity?
#' `pairwise_F_tests` compares the individual inbreeding coefficient `F`
#' ([individual_inbreeding()], expected heterozygosity from each individual's
#' own population): do they differ in inbreeding? Report the one that matches
#' your question -- a difference in FIS is a difference in inbreeding.
#'
#' @return An object of class `raddiv_het`: a list with elements
#'   `individual_heterozygosity` (one row per individual, with its
#'   heterozygosity and its `F`), `pairwise_tests` (heterozygosity, one row per
#'   pair of populations) and `pairwise_F_tests` (the same tests on `F`).
#'   Printing it shows the full report, including the missingness-confound
#'   and overdispersion checks. With `outdir`, the tables are also written to
#'   `individual_heterozygosity.<stem>.tsv`,
#'   `het_between_pops_tests.<stem>.tsv` and
#'   `het_between_pops_F_tests.<stem>.tsv`.
#' @examples
#' # A toy dataset shipped with the package (4 and 3 individuals -- far too
#' # few for a real test, which needs individuals, not loci).
#' vcf    <- system.file("extdata", "small.haps.vcf", package = "RADdiversity")
#' popmap <- system.file("extdata", "small_popmap.tsv", package = "RADdiversity")
#' res <- het_between_pops(vcf, popmap, min_call = 0.5)
#' res$individual_heterozygosity
#' res$pairwise_tests
#' print(res)   # the full report, with the confound and overdispersion checks
#' @export
het_between_pops <- function(vcf_file, popmap_f, min_call = 0.9, outdir = NULL,
                              seed = 2024, verbose = TRUE, stem = NULL) {

  ## vcf_file may be a path (checked with file.exists() below) or an
  ## already-parsed H list (see .resolve_H() in R/vcf_io.R) -- only a path
  ## needs this existence check before we try to read it.
  .check_run_inputs(vcf_file, popmap_f, stem)
  min_call <- as.numeric(min_call)
  if (is.na(min_call) || min_call < 0 || min_call > 1) stop("min_call must be in 0-1.")
  ## Restore the caller's RNG state on exit -- see the matching comment in
  ## diversity_stats().
  restore_rng <- .save_rng_state()
  on.exit(restore_rng(), add = TRUE)
  set.seed(seed)

  if (is.character(vcf_file)) message("Reading ", vcf_file, " ...")
  H    <- .resolve_H(vcf_file, verbose = verbose)
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
  cr_pop <- sweep(.typed_by_pop(H, pops), 2L, lengths(pops), "/")
  keep_pop <- cr_pop >= min_call - 1e-9

  ## Per-population individual heterozygosity and individual inbreeding F
  ## (see individual_inbreeding()), each on ITS OWN locus set.
  pop_het <- pop_F <- stats::setNames(vector("list", r), names(pops))
  for (p in names(pops)) {
    lp <- keep_pop[, p]
    a1p <- H$A1[lp, pops[[p]], drop = FALSE]; a2p <- H$A2[lp, pops[[p]], drop = FALSE]
    pop_het[[p]] <- colMeans(a1p != a2p, na.rm = TRUE)
    pop_F[[p]]   <- stats::setNames(.ind_F(a1p, a2p)$F, pops[[p]])
  }
  F_all <- unlist(unname(pop_F))
  ind$F <- round(unname(F_all[ind$sample]), 4)

  summ <- do.call(rbind, lapply(names(pops), function(pn) {
    v <- pop_het[[pn]]
    data.frame(population = pn, n = length(v), n_loci = sum(keep_pop[, pn]),
               mean_het = round(mean(v), 4), sd = round(stats::sd(v), 4),
               se = round(stats::sd(v) / sqrt(length(v)), 4),
               min = round(min(v), 4), max = round(max(v), 4),
               mean_F = round(mean(pop_F[[pn]], na.rm = TRUE), 4)) }))

  ## ---------------------------------------------------------------------------
  ## MISSINGNESS CONFOUND. An individual's heterozygosity is computed over the
  ## loci it was CALLED at. If poorly sequenced individuals are both more missing
  ## and (through allele dropout) more apparently homozygous, then a difference
  ## in heterozygosity between populations can be a difference in library
  ## quality. This is the same mechanism the rest of the package exists to catch,
  ## and it operates on the individual axis here rather than the site axis.
  ## Computed here; print.raddiv_het() reports it.
  ## ---------------------------------------------------------------------------
  cr_ind <- ncall / L
  confound <- list(no_variation = stats::sd(cr_ind) < 1e-12)
  if (!confound$no_variation) {
    ct <- suppressWarnings(stats::cor.test(cr_ind, het))
    ## and per population, so a between-population difference in call rate is
    ## not read as a within-population correlation and vice versa
    within <- do.call(rbind, lapply(names(pops), function(pn) {
      jj <- which(plab == pn)
      if (length(jj) > 3 && stats::sd(cr_ind[jj]) > 1e-12)
        data.frame(population = pn, r = suppressWarnings(stats::cor(cr_ind[jj], het[jj])),
                   min = min(cr_ind[jj]), max = max(cr_ind[jj]))
    }))
    confound <- c(confound, list(r = unname(ct$estimate), p = ct$p.value, n = length(het),
                                 within = within,
                                 pop_gap = diff(range(tapply(cr_ind, plab, mean)))))
  }

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

  ## Identity disequilibrium g2 in each population (identity_disequilibrium()):
  ## the same excess variation in individual heterozygosity, in its standard,
  ## citable form, with a bootstrap CI over individuals.
  g2 <- do.call(rbind, lapply(names(pops), function(pn) {
    lp <- keep_pop[, pn]
    cbind(population = pn,
          .g2_summary(H$A1[lp, pops[[pn]], drop = FALSE],
                      H$A2[lp, pops[[pn]], drop = FALSE], nboot = 200L, nperm = 0L))
  }))
  g2 <- data.frame(population = g2$population, n = g2$n_ind, g2 = round(g2$g2, 4),
                   g2_lo = round(g2$g2_lo, 4), g2_hi = round(g2$g2_hi, 4))

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
  rows <- rows_F <- list()
  for (i in seq_len(r - 1)) for (j in (i + 1):r) {
    p1 <- names(pops)[i]; p2 <- names(pops)[j]
    lij <- keep_pop[, p1] & keep_pop[, p2]
    nloc <- sum(lij)
    if (nloc < 50) {
      message(sprintf(
        "  %s vs %s: only %d loci clear %.0f%% call rate in BOTH populations ",
        p1, p2, nloc, 100 * min_call),
        "(need >= 50). Reported as NA for this pair; other pairs are unaffected.")
      na_row <- .na_test_row(p1, p2, nloc, length(pops[[p1]]), length(pops[[p2]]))
      rows[[length(rows) + 1]] <- na_row; rows_F[[length(rows_F) + 1]] <- na_row
      next
    }
    a1i <- H$A1[lij, pops[[p1]], drop = FALSE]; a2i <- H$A2[lij, pops[[p1]], drop = FALSE]
    a1j <- H$A1[lij, pops[[p2]], drop = FALSE]; a2j <- H$A2[lij, pops[[p2]], drop = FALSE]
    a <- colMeans(a1i != a2i, na.rm = TRUE)
    b <- colMeans(a1j != a2j, na.rm = TRUE)
    rows[[length(rows) + 1]] <- .two_sample(a, b, p1, p2, nloc, "heterozygosity")
    ## The same test on individual F (see individual_inbreeding()), over the
    ## same loci and the same individuals, with the expected heterozygosity
    ## taken from each population's own allele frequencies.
    fa <- .ind_F(a1i, a2i)$F[is.finite(a)]; fb <- .ind_F(a1j, a2j)$F[is.finite(b)]
    rows_F[[length(rows_F) + 1]] <- .two_sample(fa, fb, p1, p2, nloc, "F")
  }
  ## Benjamini-Hochberg across all k(k-1)/2 pairs, on the unrounded p-values;
  ## the returned tables round them to 3 significant digits.
  adjust <- function(t) {
    t$p_welch_BH  <- stats::p.adjust(t$p_welch,  method = "BH")
    t$p_wilcox_BH <- stats::p.adjust(t$p_wilcox, method = "BH")
    t
  }
  rounded <- function(t) {
    for (cl in c("p_welch", "p_wilcox", "p_welch_BH", "p_wilcox_BH")) t[[cl]] <- signif(t[[cl]], 3)
    t
  }
  tab <- adjust(do.call(rbind, rows)); tabF <- adjust(do.call(rbind, rows_F))
  out <- rounded(tab); outF <- rounded(tabF)

  ## Files, only when asked for (`outdir`), named by the input's stem so the
  ## SNP-VCF and haplotype-VCF runs can share a directory.
  written <- character(0)
  if (!is.null(outdir)) {
    dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
    stem <- .derive_stem(vcf_file, stem)                # populations.snps.vcf -> "snps"
    written <- file.path(outdir, sprintf(c("individual_heterozygosity.%s.tsv",
                                           "het_between_pops_tests.%s.tsv",
                                           "het_between_pops_F_tests.%s.tsv"), stem))
    utils::write.table(ind,  written[1], sep = "\t", quote = FALSE, row.names = FALSE)
    utils::write.table(out,  written[2], sep = "\t", quote = FALSE, row.names = FALSE)
    utils::write.table(outF, written[3], sep = "\t", quote = FALSE, row.names = FALSE)
  }

  structure(
    list(individual_heterozygosity = ind, pairwise_tests = out, pairwise_F_tests = outF),
    class = "raddiv_het",
    report = list(n_loci_pooled = L, n_ind = length(ids), n_pop = r,
                  min_call = min_call, summary = summ, confound = confound,
                  overdispersion = od, g2 = g2,
                  ## significance counted on the unrounded BH p-values
                  n_sig = c(welch  = sum(tab$p_welch_BH  < 0.05, na.rm = TRUE),
                            wilcox = sum(tab$p_wilcox_BH < 0.05, na.rm = TRUE)),
                  n_sig_F = c(welch  = sum(tabF$p_welch_BH  < 0.05, na.rm = TRUE),
                              wilcox = sum(tabF$p_wilcox_BH < 0.05, na.rm = TRUE)),
                  files = written))
}

## Not exported. One all-NA row of the pairwise test table.
.na_test_row <- function(p1, p2, nloc, n1, n2, mean1 = NA_real_, mean2 = NA_real_)
  data.frame(pop1 = p1, pop2 = p2, n_loci = nloc, n1 = n1, mean1 = mean1,
             n2 = n2, mean2 = mean2, diff = NA_real_, ci_lo = NA_real_,
             ci_hi = NA_real_, p_welch = NA_real_, p_wilcox = NA_real_,
             hedges_g = NA_real_)

## Not exported. Welch's t (primary), Wilcoxon (the distribution-free check)
## and Hedges' g for one pair of populations, on one number per individual
## (`a`, `b`: heterozygosity, or F). Returns one row of the pairwise table.
.two_sample <- function(a, b, p1, p2, nloc, what) {
  ## An individual can pass the earlier GLOBAL thinning check (>= 50 calls on
  ## the pooled locus set) yet still be called at NONE of THIS pair's smaller,
  ## pair-specific locus set, giving no value here. Drop such individuals from
  ## this comparison only.
  if (any(!is.finite(a)) || any(!is.finite(b))) {
    n_bad <- sum(!is.finite(a)) + sum(!is.finite(b))
    message(sprintf(
      "  %s vs %s (%s): %d individual(s) without a value on the %d shared loci ",
      p1, p2, what, n_bad, nloc), "for this pair; excluded from this comparison only.")
    a <- a[is.finite(a)]; b <- b[is.finite(b)]
  }
  if (length(a) < 2 || length(b) < 2) {
    message(sprintf(
      "  %s vs %s (%s): fewer than 2 individuals with a value in one group. ",
      p1, p2, what), "Reported as NA for this pair.")
    return(.na_test_row(p1, p2, nloc, length(a), length(b),
                        if (length(a)) round(mean(a), 4) else NA_real_,
                        if (length(b)) round(mean(b), 4) else NA_real_))
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
  ##       wilcox.test() has the same blind spot. Catch both explicitly.
  tt <- try(stats::t.test(a, b), silent = TRUE)
  if (inherits(tt, "try-error")) {
    message(sprintf(
      "  %s vs %s (%s): t.test() failed (%s). This usually means no variation ",
      p1, p2, what, trimws(conditionMessage(attr(tt, "condition")))),
      "among individuals in one or both populations, so no test is ",
      "possible. Reported as NA.")
    tt <- list(p.value = NA_real_, conf.int = c(NA_real_, NA_real_))
  } else if (is.nan(tt$p.value)) {
    message(sprintf(
      "  %s vs %s (%s): t.test() returned NaN (both samples constant at exactly ",
      p1, p2, what),
      "zero). No variation to test. Reported as NA.")
    tt <- list(p.value = NA_real_, conf.int = c(NA_real_, NA_real_))
  }
  wt <- suppressWarnings(try(stats::wilcox.test(a, b), silent = TRUE))
  if (inherits(wt, "try-error") || is.nan(wt$p.value)) wt <- list(p.value = NA_real_)
  ## Hedges' g, a sample-size-corrected standardised difference
  s <- sqrt(((length(a)-1)*stats::var(a) + (length(b)-1)*stats::var(b)) / (length(a)+length(b)-2))
  d <- if (is.finite(s) && s > 0) (mean(a) - mean(b)) / s else NA_real_
  J <- 1 - 3 / (4 * (length(a) + length(b)) - 9)
  data.frame(
    pop1 = p1, pop2 = p2, n_loci = nloc,
    n1 = length(a), mean1 = round(mean(a), 4),
    n2 = length(b), mean2 = round(mean(b), 4),
    diff = round(mean(a) - mean(b), 4),
    ci_lo = round(tt$conf.int[1], 4), ci_hi = round(tt$conf.int[2], 4),
    p_welch = tt$p.value, p_wilcox = wt$p.value,
    hedges_g = round(d * J, 2))
}

#' @rdname het_between_pops
#' @param x A `raddiv_het` object, as returned by `het_between_pops()`.
#' @param ... Ignored.
#' @export
print.raddiv_het <- function(x, ...) {
  rp <- attr(x, "report"); out <- x$pairwise_tests
  old <- options(width = max(200L, getOption("width")))
  on.exit(options(old), add = TRUE)

  cat("\n=====================================================================\n")
  cat("  HETEROZYGOSITY BETWEEN POPULATIONS -- individual as the unit\n")
  cat("  ", format(rp$n_loci_pooled, big.mark = ","), " loci pass the pooled filter (used for the\n",
      "  individual-level table and confound check below); ", rp$n_ind,
      " individuals, ", rp$n_pop, " populations\n", sep = "")
  cat("=====================================================================\n")

  cat("\nPer-population summary of individual heterozygosity\n")
  cat("  (each population's OWN locus set -- loci at >= ", 100 * rp$min_call,
      "% call rate WITHIN that population; see 'n_loci'. This can differ from\n",
      "   a specific pair's test below, which uses only loci BOTH members of\n",
      "   that pair clear -- that is expected, not an inconsistency.)\n", sep = "")
  print(rp$summary, row.names = FALSE)

  cat("\nMissingness confound check\n")
  cf <- rp$confound
  if (cf$no_variation) {
    cat("  Every individual was called at every locus; no confound possible.\n")
  } else {
    cat(sprintf("  cor(per-individual call rate, heterozygosity) = %+.3f  (p = %.3g, n = %d)\n",
                cf$r, cf$p, cf$n))
    for (k in seq_len(NROW(cf$within)))
      cat(sprintf("    within %-12s r = %+.3f   call rate %.3f-%.3f\n", cf$within$population[k],
                  cf$within$r[k], cf$within$min[k], cf$within$max[k]))
    if (is.finite(cf$p) && cf$p < 0.05 && cf$r > 0) {
      cat("  POSITIVE and significant: individuals with more missing data look LESS\n")
      cat("  heterozygous, which is the allele-dropout signature. Re-run on a\n")
      cat("  complete-data locus set (min_call = 1.0) before quoting the test\n")
      cat("  below; if the difference survives that, it is not a coverage artefact.\n")
    } else {
      cat("  No positive association, so heterozygosity is not tracking coverage.\n")
    }
    if (cf$pop_gap > 0.02)
      cat("  WARNING: mean call rate differs between populations by more than 2\n          percentage points. That asymmetry is itself a candidate explanation.\n")
  }

  cat("\n  'sd' is the spread AMONG INDIVIDUALS. That spread, divided by sqrt(n), is\n")
  cat("  the real uncertainty in a population mean -- and it does not shrink at all\n")
  cat("  as you add loci, which is why locus-based tests can be invalid here.\n")

  cat("\nOverdispersion check -- does the locus bootstrap fail on THIS dataset?\n")
  cat("  (each population's own locus set, as in the summary above)\n")
  od <- rp$overdispersion
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

  cat("\nIdentity disequilibrium g2 (David et al. 2007) -- the standard measure of\n")
  cat("variance in inbreeding among individuals; 0 when individuals do not differ.\n")
  print(rp$g2, row.names = FALSE)
  cat("  95% CI: bootstrap over individuals (see ?identity_disequilibrium). A CI above\n")
  cat("  0 means individuals differ in inbreeding, so diversity_stats()'s locus-based\n")
  cat("  intervals are too narrow for population-level inference: report its\n")
  cat("  individual-jackknife SEs (se_individuals = TRUE) instead.\n")

  npair <- nrow(out)
  cat(sprintf("\nPairwise tests -- %d comparison%s. Welch's t is the primary test;\n",
              npair, if (npair == 1) "" else "s"))
  cat("Wilcoxon is the distribution-free check. 'diff' and its CI are from Welch.\n")
  if (npair > 1)
    cat(sprintf("p_*_BH are Benjamini-Hochberg adjusted across all %d comparisons.\n", npair))
  if (npair <= 25) print(out, row.names = FALSE) else {
    o <- order(out$p_welch_BH)
    cat("  (25 most significant pairs; the full table is x$pairwise_tests)\n")
    print(out[utils::head(o, 25), ], row.names = FALSE)
  }
  cat(sprintf("\n  pairs significant at BH < 0.05: Welch %d, Wilcoxon %d, of %d\n",
              rp$n_sig[["welch"]], rp$n_sig[["wilcox"]], npair))
  cat("  Hedges' g is the standardised difference: ~0.2 small, ~0.5 medium, ~0.8 large.\n")

  outF <- x$pairwise_F_tests
  cat("\nThe same tests on individual inbreeding, F = 1 - observed/expected\n")
  cat("heterozygosity (expected from each individual's own population; see\n")
  cat("?individual_inbreeding). Heterozygosity asks whether the populations differ\n")
  cat("in DIVERSITY; F asks whether they differ in INBREEDING.\n")
  if (nrow(outF) <= 25) print(outF, row.names = FALSE) else {
    o <- order(outF$p_welch_BH)
    cat("  (25 most significant pairs; the full table is x$pairwise_F_tests)\n")
    print(outF[utils::head(o, 25), ], row.names = FALSE)
  }
  cat(sprintf("\n  pairs significant at BH < 0.05: Welch %d, Wilcoxon %d, of %d\n",
              rp$n_sig_F[["welch"]], rp$n_sig_F[["wilcox"]], nrow(outF)))

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
  if (length(rp$files)) cat("\nWrote ", paste(rp$files, collapse = " and "), "\n", sep = "")
  cat("\n")
  invisible(x)
}
