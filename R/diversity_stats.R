###############################################################################
#
#  R/diversity_stats.R -- Ho, He, FIS, rarefied allelic richness and rarefied
#  private allelic richness from a Stacks VCF, in one run.
#
#  WHY THE NUMBERS DIFFER FROM STACKS -- one line each. The full explanation,
#  and the formulas with every Stacks column matched to its counterpart here,
#  are in README.md. This header deliberately does not restate them.
#
#    He     Nei & Chesser (1983), not Stacks' `Pi`, whose correction assumes
#           FIS = 0. Both printed side by side.            [-> "Formulas"]
#    FIS    ratio of sums, not Stacks' mean of per-locus ratios, and not
#           divided by Pi.                                 [-> "Formulas"]
#    CIs    bootstrap over RAD LOCI (boot="loci", the default), not SNP rows:
#           resampling linked SNPs as independent gave intervals ~1.8x too
#           narrow. boot="individuals"/boot="both" exist for comparison
#           only -- they badly undercover He/Fis/Ar/privAr (see
#           README.md, "Bootstrap mode"). Do NOT use any of these
#           to compare populations -- that is het_between_pops(). [-> "Step 3"]
#    sites  DEFAULT: a locus is used by a population once that population has
#           >= min_n typed individuals there (available-data, per population,
#           per locus). Pass complete_case = TRUE for the old rule: kept only
#           where EVERY individual of EVERY population is genotyped (Schmidt
#           et al. 2021). Watch the retention lines either way.
#
#  Engine: this package's own estimators (R/estimators.R). With
#  hierfstat_check = TRUE they are cross-checked against hierfstat when it is
#  installed. Validation provenance is in README.md, Step 0.
#
###############################################################################

#' Diversity statistics from a Stacks VCF
#'
#' Computes Ho, He, FIS (Nei & Chesser 1983, ratio of sums), rarefied allelic
#' richness and rarefied private allelic richness from a Stacks 2
#' `populations` VCF and popmap, with a block bootstrap over RAD loci and a
#' delete-one-block jackknife standard error. Returns the result tables;
#' printing the result (at the console, or with `print()`) shows the same
#' report as the `diversity_stats.R` command-line script, and `outdir` writes
#' the tables to TSV files.
#'
#' Run it twice: `populations.snps.vcf` gives Ho, He, `pct_poly` and the
#' per-sequenced-site values; `populations.haps.vcf` gives FIS, Ar and
#' privAr. Each run prints which of its numbers to take. See `README.md`
#' ("Which file for which statistic") for the full reasoning.
#'
#' @param vcf_file Path to `populations.snps.vcf` or `populations.haps.vcf`
#'   (optionally gzip-compressed) -- OR, an already-parsed (and optionally
#'   filtered) `H` list, i.e. the object returned by [read_stacks_vcf()], on
#'   its own or passed through one or more `filter_*()` functions first (see
#'   `R/filter_loci.R`). When passing a list, `stem` must also be given (see
#'   below), since the usual output-filename logic needs a real filename to
#'   work from.
#' @param popmap_f Path to a two-column, no-header popmap TSV (`sample_id
#'   <TAB> population`).
#' @param g Rarefaction size in GENE COPIES (10 diploids = 20). Required, no
#'   default -- must be `<= 2x` the smallest population.
#' @param nboot Bootstrap replicates. Default `10000L`; `0` = none.
#' @param boot Which axis the bootstrap resamples: one of `"loci"` (RAD loci,
#'   the default -- matches hierfstat), `"individuals"` (within each
#'   population, matches diveRsity), or `"both"` (Owen & Eckles 2012
#'   crossed-factor scheme) -- matched exactly, not partially. `"individuals"`/
#'   `"both"` are comparison modes only -- an empirical coverage simulation
#'   found they badly undercover the true value for every statistic with a
#'   finite-sample correction; see `README.md`, "Bootstrap mode".
#' @param sites Total sequenced sites, for the autosomal Ho/He (nucleotide
#'   diversity) conversion -- adds `Ho_autosomal`/`He_autosomal` alongside
#'   Ho/He (SNP VCF only; ignored, with a message, on a haplotype VCF).
#'   Accepts four shapes: a single non-negative number, broadcast to every
#'   population (the default `0` means "off"); a named numeric vector, one
#'   value per population (names must match the popmap's population names);
#'   a data frame or matrix with a population column (`population`/`pop`)
#'   and a sites column (`sites`); or a path to a Stacks
#'   `populations.sumstats_summary.tsv` file, from which each population's
#'   own `Sites` (the "All positions (variant and fixed)" block) is read
#'   automatically via [read_sumstats_summary()]. The vector/data-frame/path
#'   forms let a population with more missing loci -- which is exactly what
#'   makes its own `Sites` count smaller in Stacks' own accounting -- use
#'   its own denominator instead of another population's.
#' @param min_n Available-data mode only: minimum typed individuals a
#'   population needs at a locus to use that locus for that population.
#'   Default `2L`, the mathematical floor below which Hs/FIS are undefined.
#'   Ignored when `complete_case = TRUE`.
#' @param complete_case Use the old, stricter rule: a record is used only if
#'   every individual of every population is genotyped there (Schmidt et al.
#'   2021 recommendation (b)). Default `FALSE`.
#' @param outdir Directory to write the result tables into as TSV files
#'   (created if needed). Default `NULL`: write no files -- the tables are in
#'   the returned object. The command-line script writes to the current
#'   directory.
#' @param seed Random seed set before bootstrapping, for reproducibility.
#'   Default `2024`. The caller's own RNG state (e.g. from their own prior
#'   `set.seed()`) is restored when this function returns, so calling it
#'   interactively does not affect subsequent random draws in the caller's
#'   session.
#' @param verbose Print the VCF/popmap parsing summary. Default `TRUE`.
#'   Progress messages can be silenced with `suppressMessages()`; the report
#'   itself appears only when the result is printed.
#' @param stem Text used to build the two output filenames (e.g.
#'   `diversity_per_population.<stem>.tsv`) -- ordinarily derived
#'   automatically from `vcf_file`'s own filename (e.g.
#'   `populations.haps.vcf` gives `stem = "haps"`), which is why this
#'   package's documented workflow runs this function once on the haplotype
#'   VCF and once on the SNP VCF without the second run silently overwriting
#'   the first run's files. That automatic derivation has no filename to
#'   work from when `vcf_file` is an already-parsed list rather than a path,
#'   so `stem` must be supplied explicitly in that case (e.g.
#'   `stem = "snps"`) -- this function stops with an explanatory error
#'   rather than guessing. Default `NULL` (derive automatically from
#'   `vcf_file` when it is a path; required otherwise). Supplying `stem`
#'   alongside a path overrides the automatic derivation.
#' @param hierfstat_check If `TRUE` and the `hierfstat` package is installed,
#'   also compute per-locus Ho/Hs with `hierfstat::basic.stats()` and
#'   rarefied allelic richness with `hierfstat::allelic.richness()`, and
#'   report the largest difference from this package's own values (with a
#'   warning if they disagree). Default `FALSE`: every number is computed by
#'   this package either way, the package tests already make this comparison,
#'   and on large datasets the hierfstat calls were most of the run time. For
#'   Weir & Cockerham FST, use [differentiation_stats()].
#' @param se_individuals Also report delete-one-INDIVIDUAL jackknife standard
#'   errors (`Ho_se_ind`, `He_se_ind`, `Fis_se_ind`, `Ar_se_ind`,
#'   `privAr_se_ind`): the uncertainty from which individuals were sampled,
#'   which the locus-based `_se`, `_lo` and `_hi` leave out. They matter when
#'   individuals differ in inbreeding or include relatives -- check with
#'   [identity_disequilibrium()] (g2 > 0) -- and are then the ones to report
#'   for population-level inference. The Ar/privAr ones need
#'   `g <= 2 * (n - 1)` for every population and are `NA` otherwise. Default
#'   `FALSE`.
#' @return An object of class `raddiv_diversity`: a list with elements
#'   `per_population` (Ho, He, FIS and % polymorphic, with jackknife SEs and
#'   bootstrap CIs), `richness` (Ar, privAr and priv_total) and `autosomal`
#'   (per-sequenced-site Ho/He; `NULL` unless `sites` gave at least one
#'   population a plausible value, and always `NULL` on a haplotype VCF).
#'   Printing it shows the full report. With `outdir`, the tables are also
#'   written to `diversity_per_population.<stem>.tsv`,
#'   `diversity_richness.<stem>.tsv` and (when present)
#'   `diversity_autosomal.<stem>.tsv`.
#' @examples
#' # A toy dataset shipped with the package: 80 RAD loci, 2 populations of 4
#' # and 3 individuals (real studies need far more individuals).
#' haps   <- system.file("extdata", "small.haps.vcf", package = "RADdiversity")
#' snps   <- system.file("extdata", "small.snps.vcf", package = "RADdiversity")
#' popmap <- system.file("extdata", "small_popmap.tsv", package = "RADdiversity")
#'
#' # Haplotype VCF: FIS, allelic richness, private allelic richness.
#' res <- diversity_stats(haps, popmap, g = 4, nboot = 100)
#' res$per_population
#' res$richness
#'
#' # SNP VCF with each population's own sequenced-site count, read from Stacks'
#' # populations.sumstats_summary.tsv: Ho and He per sequenced site.
#' sumstats <- system.file("extdata", "populations.sumstats_summary.tsv",
#'                         package = "RADdiversity")
#' res_snps <- diversity_stats(snps, popmap, g = 4, nboot = 0, sites = sumstats)
#' res_snps$autosomal
#'
#' # The full report, as the command-line script prints it:
#' print(res_snps)
#' @export
diversity_stats <- function(vcf_file, popmap_f, g, nboot = 10000L,
                             boot = "loci",
                             sites = 0, min_n = 2L, complete_case = FALSE,
                             outdir = NULL, seed = 2024, verbose = TRUE,
                             stem = NULL, hierfstat_check = FALSE,
                             se_individuals = FALSE) {

  g <- as.integer(g)
  if (is.na(g)) stop("g must be an integer (gene copies).")
  nboot <- as.integer(nboot)
  if (is.na(nboot) || nboot < 0) stop("nboot must be a non-negative integer.")
  ## Shape/type check only -- population names aren't known yet (the popmap
  ## hasn't been read), so a bad `sites` (wrong type, a negative number, a
  ## nonexistent file path) still fails here, before the potentially slow
  ## VCF read below. Per-population name-matching happens later, once `pops`
  ## exists (see .resolve_sites(), just below this function).
  .resolve_sites(sites)
  min_n <- as.integer(min_n)
  if (is.na(min_n) || min_n < 2)
    stop("min_n must be an integer >= 2 (Hs/FIS are undefined below n = 2).")
  ## Deliberately an exact match, not match.arg()'s partial-prefix matching:
  ## "individuals"/"both" are comparison-only modes this package's own
  ## documentation says badly undercover the true value and must not be
  ## treated as publication-ready CIs, so a typo like `boot="individual"`
  ## must be rejected loudly rather than silently resolving to a real,
  ## different (and non-default) mode.
  if (!(boot %in% c("loci", "individuals", "both")))
    stop("boot must be one of: loci, individuals, both (got: ", boot, ")")
  boot_mode <- boot
  ## vcf_file may be a path (checked with file.exists() below) or an
  ## already-parsed H list (see .resolve_H() in R/vcf_io.R) -- only a path
  ## needs this existence check before we try to read it. A list also needs
  ## `stem` supplied up front (checked here, before any of the potentially
  ## slow work below runs), since the output filenames later in this
  ## function have no other way to know what to call themselves.
  .check_run_inputs(vcf_file, popmap_f, stem)
  ## Restore the caller's RNG state on exit -- this function is meant to be
  ## called interactively (not just as a fresh Rscript process), so
  ## set.seed() below must not silently overwrite the caller's own
  ## random-number stream for the rest of their session.
  restore_rng <- .save_rng_state()
  on.exit(restore_rng(), add = TRUE)
  set.seed(seed)

  have_hf <- isTRUE(hierfstat_check) && .hierfstat_available()
  if (isTRUE(hierfstat_check) && !have_hf)
    message("hierfstat_check = TRUE but hierfstat is not installed; skipping the cross-check.")

  if (is.character(vcf_file)) message("Reading ", vcf_file, " ...")
  H <- .resolve_H(vcf_file, verbose = verbose)
  pops <- read_popmap(popmap_f, H$samples, verbose = verbose)
  r <- length(pops)
  if (r < 2) stop("Need at least 2 populations.")
  ## Full resolution, now that population names are known: a named vector,
  ## one sequenced-site count per population, in `names(pops)` order. See
  ## .resolve_sites() for the four accepted `sites` shapes.
  sites_vec <- .resolve_sites(sites, names(pops))
  nmax <- lengths(pops)
  ## A population of one individual has no defined Hs, and every statistic
  ## downstream then silently reads as 0 rather than as missing, because
  ## colSums(na.rm = TRUE) drops the NA from the numerator while the record is
  ## still counted in the denominator. Refuse rather than print zeros.
  tiny <- names(pops)[nmax < 2]
  if (length(tiny))
    stop("Population(s) with fewer than 2 individuals: ", paste(tiny, collapse = ", "),
         "\n  Expected heterozygosity is undefined at n = 1, and rarefaction to ",
         "2 gene\n  copies is meaningless. Drop these populations or merge them.")

  ## Check g against the population sizes now, right after they're known,
  ## rather than after the (potentially slow) per-locus work below -- a bad
  ## g value should fail in seconds, not minutes.
  if (g > 2L * min(nmax))
    stop("g = ", g, " gene copies exceeds the smallest population's ",
         2L * min(nmax), ". g is in GENE COPIES: 10 diploids = 20.")
  if (g < 2) stop("g must be at least 2 gene copies.")

  ## Is this a haplotype VCF? Two independent signals: Stacks writes multi-
  ## nucleotide alleles ("AC", "CA") for haplotype records and single bases for
  ## SNP records, and haplotype records are frequently multi-allelic. Either one
  ## is decisive. This matters because the autosomal conversion below is only
  ## nucleotide diversity when one record is one SITE.
  is_hapvcf <- .is_haplotype_H(H)
  message(sprintf("  record type: %s",
                  if (is_hapvcf) "HAPLOTYPE (one multi-allelic locus per RAD tag)"
                  else "SNP (one site per record)"))

  ## RAD locus of each record: several SNPs share one in a SNP VCF, one per
  ## record in a haplotype VCF. This is the bootstrap resampling unit.
  rad <- H$locus_raw
  n_rec <- nrow(H$A1)
  message(sprintf("  %s records on %s RAD loci (%.2f per locus)",
                  format(n_rec, big.mark = ","),
                  format(length(unique(rad)), big.mark = ","),
                  n_rec / length(unique(rad))))

  ## ---------------------------------------------------------------------------
  ## Locus inclusion.
  ##
  ## `n_typed[j, p]` = individuals of population p genotyped at record j. Kept
  ## for every record (not just the ones ultimately used) because it drives the
  ## per-cell threshold below, the reporting, and (for complete_case) the old
  ## joint filter.
  ## ---------------------------------------------------------------------------
  n_typed <- .typed_by_pop(H, pops)

  cells_used <- cells_total <- NA_real_
  if (complete_case) {
    ## OLD rule: a record is used only if EVERY individual of EVERY population
    ## is genotyped there (Schmidt et al. 2021 recommendation (b)) -- a JOINT
    ## condition across all populations on the SAME record, so one population's
    ## missing individual costs every other population that record too.
    keep <- rep(TRUE, n_rec)
    for (p in names(pops)) keep <- keep & (n_typed[, p] == nmax[[p]])
    message(sprintf("Complete-data locus set (complete_case): %s of %s records (%.1f%%)",
                    format(sum(keep), big.mark = ","), format(n_rec, big.mark = ","),
                    100 * mean(keep)))
    if (sum(keep) < 0.05 * n_rec)
      message("  WARNING: under 5% survive. Report a looser threshold alongside this.")
    if (!sum(keep)) stop("No record is genotyped in every individual. With many ",
                         "populations this is common; drop complete_case or drop the worst individuals.")
  } else {
    ## DEFAULT rule: available data, decided PER POPULATION, independently,
    ## locus by locus. A record is worth keeping at all once >= 1 population
    ## clears min_n there; which population's cells actually get used is
    ## decided later, per (locus, population) cell (see the Ho/Hs loop below).
    keep <- rowSums(n_typed >= min_n) > 0
    cells_used  <- sum(n_typed >= min_n)
    cells_total <- n_rec * r
    message(sprintf("Available-data mode (min_n = %d): %s of %s records usable by >= 1 population",
                    min_n, format(sum(keep), big.mark = ","), format(n_rec, big.mark = ",")))
    message(sprintf("  %s of %s locus-by-population cells meet min_n (%.1f%%)",
                    format(cells_used, big.mark = ","), format(cells_total, big.mark = ","),
                    100 * cells_used / cells_total))
    for (p in names(pops)) {
      ok_p <- n_typed[, p] >= min_n
      ## Typed-n summary only where the population has any usable locus:
      ## min() of nothing is Inf, which sprintf("%d") cannot print, and the
      ## all-missing case must reach the clear error just below instead.
      message(sprintf("    %-22s %s of %s loci usable%s",
                      p, format(sum(ok_p), big.mark = ","), format(n_rec, big.mark = ","),
                      if (any(ok_p))
                        sprintf(" (typed-n there: mean %.1f, min %d, max %d)",
                                mean(n_typed[ok_p, p]), min(n_typed[ok_p, p]),
                                max(n_typed[ok_p, p]))
                      else ""))
    }
    if (!sum(keep)) stop("No record has >= min_n typed individuals in any population. ",
                         "Lower min_n (currently ", min_n, ") or check that genotypes were parsed.")
  }
  idx <- which(keep); L <- length(idx)
  rad_k <- rad[idx]; uloc <- unique(rad_k); li <- match(rad_k, uloc); nL <- length(uloc)

  ## ---------------------------------------------------------------------------
  ## Per-locus Ho, Hs and Hp, and the allele-count table everything later is
  ## built from:
  ##   Ho[j, i]  fraction of population i's typed individuals heterozygous at j
  ##   Hs[j, i]  Nei & Chesser (1983) gene diversity, hs_from_counts()
  ##   Hp[j, i]  the SAME quantity with the estimator behind Stacks' `Pi`
  ##             (gene_div_2n_counts()), so the two can be printed side by side
  ##             and the reader can see exactly where they diverge; not used
  ##             for anything downstream
  ## `lown` marks cells with fewer than min_n typed individuals: below 2, Hs is
  ## undefined; above that, min_n is the caller's own floor. Under
  ## complete_case every used cell is fully typed, so the floor is just 2 and
  ## the mask is empty.
  ## ---------------------------------------------------------------------------
  message("Engine: internal (Nei & Chesser Hs, R/estimators.R)",
          if (have_hf) "; hierfstat cross-check requested" else "")
  min_n_eff <- if (complete_case) 2L else min_n
  lown <- n_typed[idx, , drop = FALSE] < min_n_eff
  Ho <- Hs <- Hp <- matrix(NA_real_, L, r, dimnames = list(NULL, names(pops)))
  for (i in seq_len(r))    # NaN where nobody is typed; masked by `lown` below
    Ho[, i] <- rowMeans(H$A1[idx, pops[[i]], drop = FALSE] !=
                        H$A2[idx, pops[[i]], drop = FALSE], na.rm = TRUE)
  n_all <- matrix(NA_integer_, L, r, dimnames = list(NULL, names(pops)))
  cmats <- vector("list", L)
  for (j in seq_len(L)) {
    k <- H$n_alleles[idx[j]]
    m <- matrix(0L, r, k, dimnames = list(names(pops), NULL))
    for (i in seq_len(r)) {
      v <- c(H$A1[idx[j], pops[[i]]], H$A2[idx[j], pops[[i]]])
      m[i, ] <- tabulate(v[!is.na(v)], nbins = k)
      Hp[j, i] <- gene_div_2n_counts(m[i, ])
      if (!lown[j, i]) Hs[j, i] <- hs_from_counts(m[i, ], Ho[j, i], n_typed[idx[j], i])
    }
    cmats[[j]] <- m; n_all[j, ] <- rowSums(m > 0)
  }
  Ho[lown] <- NA_real_

  ## ---------------------------------------------------------------------------
  ## Rarefied allelic and private allelic richness, over GENE COPIES
  ## (g checked and set right after nmax was computed, above)
  ## ---------------------------------------------------------------------------
  message(sprintf("Rarefying to g = %d gene copies (smallest population has %d)",
                  g, 2L * min(nmax)))
  Ar <- Pr <- matrix(NA_real_, L, r, dimnames = list(NULL, names(pops)))
  for (j in seq_len(L)) {
    for (i in seq_len(r)) Ar[j, i] <- rare_richness(cmats[[j]][i, ], g)
    Pr[j, ] <- rare_private_all(cmats[[j]], g)
  }
  ## Enforce `min_n` on Ar/Pr too, mirroring the Ho/Hs mask above. Without
  ## this, min_n changes Ho/He/Fis but silently has no effect on allelic
  ## richness: `cmats` above is built from every typed individual regardless of
  ## min_n, so a locus/population cell with, say, only 1 typed individual
  ## would otherwise still contribute an Ar/privAr value even after raising
  ## min_n past that.
  Ar[lown] <- NA_real_; Pr[lown] <- NA_real_

  ## Optional cross-check against hierfstat (hierfstat_check = TRUE). Every
  ## number reported comes from this package; this only prints how far
  ## hierfstat's own values are from them.
  if (have_hf) {
    dat <- .to_hierfstat_df(H, pops, idx)
    bs <- try(hierfstat::basic.stats(dat, diploid = TRUE, digits = 12), silent = TRUE)
    if (!inherits(bs, "try-error")) {
      d_ho <- max(abs(as.matrix(bs$Ho)[!lown] - Ho[!lown]), na.rm = TRUE)
      d_hs <- max(abs(as.matrix(bs$Hs)[!lown] - Hs[!lown]), na.rm = TRUE)
      message(sprintf("  cross-check vs hierfstat::basic.stats(): max diff Ho %.2e, Hs %.2e",
                      d_ho, d_hs))
      if (max(d_ho, d_hs) > 1e-8) warning("Ho/Hs disagree with hierfstat::basic.stats().")
    }
    hfar <- try(hierfstat::allelic.richness(dat, min.n = g)$Ar, silent = TRUE)
    ## hierfstat pads a "dummy.loc" row onto its result when given exactly one
    ## locus, so hfar and Ar disagree in shape and cannot be diffed. That only
    ## happens when very few records survived the locus filter -- a dataset too
    ## thin to trust anyway -- so warn instead of crashing.
    if (inherits(hfar, "try-error") ||
        !identical(dim(as.matrix(hfar)), dim(Ar))) {
      warning("Skipped the hierfstat::allelic.richness() cross-check: its ",
              "result has a different shape than this run's (", L,
              " locus/loci). This usually means too few loci survived the ",
              "locus filter to trust these numbers.")
    } else {
      d <- max(abs(as.matrix(hfar) - Ar), na.rm = TRUE)
      message(sprintf("  cross-check vs hierfstat::allelic.richness(): max diff %.2e", d))
      if (d > 1e-8) warning("Allelic richness disagrees with hierfstat.")
    }
  }

  ## ---------------------------------------------------------------------------
  ## Point estimates and a BLOCK bootstrap over RAD loci
  ## ---------------------------------------------------------------------------
  poly <- n_all > 1
  ## Numerator and denominator must count the SAME records. Summing with
  ## na.rm = TRUE while dividing by length(sel) silently turns an undefined
  ## locus into a zero, which then looks like a measurement.
  usable <- is.finite(Ho) & is.finite(Hs)
  n_drop <- sum(!usable)
  if (n_drop)
    message(sprintf("  %s locus-by-population cells have no defined Ho/Hs and are ",
                    format(n_drop, big.mark = ",")),
            "excluded from that population's means (not counted as zero).")
  ## ---------------------------------------------------------------------------
  ## Every statistic reported below is a ratio (or a plain sum) of COLUMN SUMS
  ## over records -- so the calculation needs only these per-record pieces:
  ##   n    1 if Ho and Hs are both defined        ho, hs  those Ho, Hs (else 0)
  ##   n2   1 if Hp (Stacks-Pi form) is too        ho2, hp Ho, Hp on that mask
  ##   nP   1 if privAr is defined                 P       privAr (else 0)
  ##   nA   1 if Ar is defined                     A       Ar (else 0)
  ##   pol  1 if polymorphic here (among the n)
  ## An undefined value contributes 0 to its sum AND 0 to its count, so a
  ## numerator and its denominator always count the same records -- summing
  ## with na.rm while dividing by the number of records would silently turn an
  ## undefined locus into a measured zero. The Stacks-Pi comparison (ho2, hp)
  ## gets its own mask for the same reason.
  ##
  ## `S` sums those pieces once per RAD LOCUS (one row per locus). The point
  ## estimate, the jackknife and the loci bootstrap are all computed from S --
  ## see R/resampling.R. stat_mat() turns summed pieces into statistics.
  ## ---------------------------------------------------------------------------
  usp <- usable & is.finite(Hp)
  z <- function(x, ok) ifelse(ok, x, 0)
  rec <- cbind(usable, z(Ho, usable), z(Hs, usable),        # n, ho, hs
               usp, z(Ho, usp), z(Hp, usp),                 # n2, ho2, hp
               !is.na(Pr), z(Pr, !is.na(Pr)),               # nP, P
               !is.na(Ar), z(Ar, !is.na(Ar)),               # nA, A
               z(poly, usable))                             # pol
  storage.mode(rec) <- "double"
  S <- rowsum(rec, li, reorder = TRUE)          # nL x (11 * r); row b = RAD locus b
  stat_names <- paste0(rep(c("Ho_", "He_", "Fis_", "He2n_", "Fis2n_", "Ar_", "Pr_",
                             "poly_", "PrTot_"), each = r), names(pops))
  ## Summed pieces (one row per estimate or replicate) -> statistics. A mean or
  ## ratio whose count is 0 is NA ("no information"), never the NaN that 0/0
  ## would print. priv_total is the plain sum of privAr over the loci where it
  ## is defined; a population with NO defined locus gets NA, which is a
  ## different claim from a population with defined loci and zero private
  ## alleles. The FIS guards use > 1e-12 rather than > 0 only to absorb the
  ## last-bit rounding a jackknife leaves when it subtracts one locus from the
  ## total (an exact 0 can come back as 1e-17).
  stat_mat <- function(tot) {
    tot <- matrix(tot, ncol = 11L * r)
    blk <- function(k) tot[, (k - 1L) * r + seq_len(r), drop = FALSE]
    n  <- blk(1); ho  <- blk(2); hs <- blk(3)
    n2 <- blk(4); ho2 <- blk(5); hp <- blk(6)
    nP <- blk(7); P   <- blk(8); nA <- blk(9); A <- blk(10); pol <- blk(11)
    per <- function(x, cnt) { out <- x / cnt; out[cnt < 0.5] <- NA_real_; out }
    PrTot <- P; PrTot[nP < 0.5] <- NA_real_
    out <- cbind(per(ho, n), per(hs, n),
                 ifelse(hs > 1e-12, 1 - ho / hs, NA_real_),
                 per(hp, n2),
                 ifelse(hp > 1e-12, 1 - ho2 / hp, NA_real_),
                 per(A, nA), per(P, nP), 100 * per(pol, n), PrTot)
    colnames(out) <- stat_names
    out
  }


  ## The boot="individuals"/boot="both" engine: recomputes every statistic
  ## from a fresh per-population individual reweighting, for a given locus
  ## multiset `sel` (either every locus once, or a resampled multiset with
  ## repeats, as one loci-bootstrap draw would give). Reuses hs_from_counts(),
  ## gene_div_2n_counts(), rare_richness(), rare_private_all() unmodified --
  ## only how the input counts are built differs.
  ##
  ## Mechanism (Owen & Eckles 2012 product-weight bootstrap, r=2 crossed
  ## factors: locus x individual -- see README.md, "Bootstrap mode"):
  ## draw ONE fresh individual weight vector per population (step 1), reuse it
  ## for every locus in `sel` including repeats (step 2), which is equivalent
  ## to duplicating each resampled individual's row that many times and
  ## recomputing directly (Rubin 1981 exchangeable weighting) -- cheaper here
  ## because a weight vector is drawn once per replicate, not once per locus.
  ##
  ## Uses the same internal formulas as the point estimate, fed the
  ## resampled (weighted) counts.
  ##
  ## COVERAGE CAVEAT (why this is not the default): a duplicated individual
  ## (weight > 1) makes the resample's He/Fis/Ar/privAr come out biased low,
  ## because those statistics carry a finite-sample correction (Nei-Chesser's
  ## n/(n-1) for He/Fis; the hypergeometric rarefaction formula's dependence
  ## on total gene copies for Ar/privAr) that assumes n DISTINCT draws -- a
  ## with-replacement resample violates that. An empirical coverage simulation
  ## (known-truth He/Fis/Ar, n=10-15, nboot=500, 150-300 simulated datasets)
  ## found this drives coverage of the TRUE value as low as 0%, and it WORSENS
  ## as locus count grows (the bias is fixed; the bootstrap spread shrinks).
  ## Ho (no correction factor) was unaffected (88-92% coverage). See
  ## README.md, "Bootstrap mode" for the full numbers. boot="loci" is
  ## therefore the default for every statistic; this function only runs under
  ## an explicit boot="individuals"/"both" comparison request.
  stat_from_resampled <- function(sel) {
    w <- lapply(seq_len(r), function(i)
      tabulate(sample.int(nmax[i], nmax[i], replace = TRUE), nbins = nmax[i]))

    ## Compute once per DISTINCT locus row in `sel`, then weight by how many
    ## times that row appears in `sel` when aggregating below -- mathematically
    ## identical to summing over `sel` with repeats (each repeat sees the same
    ## locus data and the same replicate's weight vectors), cheaper to compute.
    urow <- unique(sel)
    mult <- tabulate(match(sel, urow), nbins = length(urow))
    nU <- length(urow)

    Ho_w <- Hs_w <- Hp_w <- Ar_w <- Pr_w <- matrix(NA_real_, nU, r,
                                                    dimnames = list(NULL, names(pops)))
    poly_w <- matrix(FALSE, nU, r, dimnames = list(NULL, names(pops)))

    for (k in seq_len(nU)) {
      j <- urow[k]
      na <- H$n_alleles[idx[j]]
      wm <- matrix(0, r, na, dimnames = list(names(pops), NULL))
      for (i in seq_len(r)) {
        a <- H$A1[idx[j], pops[[i]]]; b <- H$A2[idx[j], pops[[i]]]
        ok <- !is.na(a)
        wi <- w[[i]]
        typed_w <- sum(wi[ok])
        if (typed_w < min_n_eff) next
        het <- a[ok] != b[ok]
        ho_w <- sum(wi[ok][het]) / typed_w
        cnt <- tabulate(rep.int(a[ok], wi[ok]), nbins = na) +
               tabulate(rep.int(b[ok], wi[ok]), nbins = na)
        wm[i, ] <- cnt
        Ho_w[k, i] <- ho_w
        Hs_w[k, i] <- hs_from_counts(cnt, ho_w, typed_w)
        Hp_w[k, i] <- gene_div_2n_counts(cnt)
        poly_w[k, i] <- sum(cnt > 0) > 1
      }
      for (i in seq_len(r)) Ar_w[k, i] <- rare_richness(wm[i, ], g)
      Pr_w[k, ] <- rare_private_all(wm, g)
    }

    usable_w <- is.finite(Ho_w) & is.finite(Hs_w)
    nsp <- colSums(usable_w * mult); nsp[nsp == 0] <- NA
    ho <- colSums(ifelse(usable_w, Ho_w, 0) * mult)
    hs <- colSums(ifelse(usable_w, Hs_w, 0) * mult)
    usp  <- usable_w & is.finite(Hp_w)
    nspp <- colSums(usp * mult); nspp[nspp == 0] <- NA
    hop <- colSums(ifelse(usp, Ho_w, 0) * mult)
    hpp <- colSums(ifelse(usp, Hp_w, 0) * mult)

    ar_ok <- is.finite(Ar_w)
    ar_cnt <- colSums(ar_ok * mult); ar_cnt[ar_cnt == 0] <- NA
    Ar_ <- colSums(ifelse(ar_ok, Ar_w, 0) * mult) / ar_cnt

    pr_ok <- is.finite(Pr_w)
    pr_cnt <- colSums(pr_ok * mult)
    PrSum_ <- colSums(ifelse(pr_ok, Pr_w, 0) * mult)
    pr_cnt_na <- pr_cnt; pr_cnt_na[pr_cnt_na == 0] <- NA
    Pr_ <- PrSum_ / pr_cnt_na
    PrTot_ <- PrSum_; PrTot_[pr_cnt == 0] <- NA_real_

    c(stats::setNames(ho / nsp, paste0("Ho_", names(pops))),
      stats::setNames(hs / nsp, paste0("He_", names(pops))),
      stats::setNames(ifelse(hs > 0, 1 - ho / hs, NA_real_), paste0("Fis_", names(pops))),
      stats::setNames(hpp / nspp, paste0("He2n_", names(pops))),
      stats::setNames(ifelse(hpp > 0, 1 - hop / hpp, NA_real_), paste0("Fis2n_", names(pops))),
      stats::setNames(Ar_, paste0("Ar_", names(pops))),
      stats::setNames(Pr_, paste0("Pr_", names(pops))),
      stats::setNames(100 * colSums(ifelse(usable_w, poly_w, 0) * mult) / nsp,
               paste0("poly_", names(pops))),
      stats::setNames(PrTot_, paste0("PrTot_", names(pops))))
  }

  if (!any(poly, na.rm = TRUE))
    stop("No locus is polymorphic in any population, so He, FIS and allelic\n",
         "  richness are all degenerate. Check that genotypes were parsed: the\n",
         "  missing-genotype rate reported above should not be near 100%.")
  point <- stat_mat(colSums(S))[1, ]

  ## ---------------------------------------------------------------------------
  ## Delete-one-BLOCK jackknife over RAD loci -- a standard error alongside the
  ## bootstrap CI, computed the way population-genetics software has long done
  ## it (Weir 1996, "Genetic Data Analysis II"; the same method behind
  ## GENEPOP/FSTAT's Fst/Fis standard errors), on the SAME loci axis this
  ## project's own coverage simulation validated. Always uses the per-locus
  ## block sums (the boot="loci" mechanism) regardless of which boot mode is
  ## active: jackknife-over-loci is a cross-check for that one validated axis,
  ## not a per-mode option.
  ##
  ## Unlike a with-replacement bootstrap, a jackknife replicate never
  ## duplicates anything -- it just drops one whole RAD-locus block -- so it
  ## does NOT have the compositional-duplication bias that made individual-
  ## resampling undercover He/Fis/Ar (see "Bootstrap mode" in
  ## README.md). It needs no random numbers and no nboot: with nL
  ## RAD loci there are exactly nL delete-one replicates, deterministic and
  ## reproducible run to run.
  ##
  ##   theta_bar         = mean over the nL delete-one estimates
  ##   SE_jackknife       = sqrt( (nL-1)/nL * sum_b (theta_(-b) - theta_bar)^2 )
  ## All nL delete-one estimates come from the block sums at once (see
  ## .jack_block_sums() in R/resampling.R).
  se_jk <- .jack_block_sums(S, stat_mat)

  ## Optional delete-one-INDIVIDUAL jackknife (see .jack_individuals() in
  ## R/resampling.R). Its full-sample recomputation must reproduce the point
  ## estimates exactly -- a built-in check that the two code paths agree.
  se_ind <- NULL
  if (isTRUE(se_individuals)) {
    message("Delete-one-individual jackknife over ", format(sum(nmax), big.mark = ","),
            " individuals ...")
    ji <- .jack_individuals(H, pops, idx, cmats, n_typed[idx, , drop = FALSE], min_n_eff, g)
    se_ind <- ji$se
    gap <- suppressWarnings(max(abs(ji$full - point[names(ji$full)]), na.rm = TRUE))
    if (is.finite(gap) && gap > 1e-8)
      warning("Individual jackknife: its full-sample estimates differ from the point ",
              "estimates by ", signif(gap, 3), ". Please report this.")
    if (g > 2L * (min(nmax) - 1L))
      message("  Ar/privAr individual SEs are NA for populations with fewer than ",
              "g/2 + 1 individuals: removing one leaves fewer than g = ", g,
              " gene copies. Use g <= ", 2L * (min(nmax) - 1L), " to get them.")
  }

  ## Ar/Pr coverage: how many loci actually contributed a defined value per
  ## population. rare_richness() is single-population (a locus counts for
  ## population i once i itself has >= g gene copies there), but
  ## rare_private_all() is a JOINT statistic across ALL populations at once --
  ## "private" only means something when every population's presence/absence
  ## at that locus is defined at the SAME rarefaction depth g, so one
  ## population short of g gene copies at a locus makes privAr undefined for
  ## every population there, not just the short one. Under complete_case this
  ## was unreachable (every used locus had exactly nmax gene copies everywhere,
  ## and g <= 2*min(nmax) by the check above); under available data, per-locus
  ## coverage varies, so report it rather than let it silently shrink privAr's
  ## effective locus set.
  ar_n <- colSums(!is.na(Ar[seq_len(L), , drop = FALSE]))
  pr_n <- colSums(!is.na(Pr[seq_len(L), , drop = FALSE]))
  if (nboot > 0) {
    boot_desc <- switch(boot_mode,
      loci        = sprintf("%s RAD loci", format(nL, big.mark = ",")),
      individuals = "individuals within each population",
      both        = sprintf("%s RAD loci AND individuals within each population",
                            format(nL, big.mark = ",")))
    message(sprintf("Bootstrapping %s replicates over %s (boot=\"%s\") ...",
                    format(nboot, big.mark = ","), boot_desc, boot_mode))
    rows_of <- split(seq_len(L), li)
    resample_loci <- function() unlist(rows_of[sample.int(nL, nL, replace = TRUE)], use.names = FALSE)
    bt <- switch(boot_mode,
      loci        = t(.boot_block_sums(S, nboot, stat_mat)),
      individuals = replicate(nboot, stat_from_resampled(seq_len(L))),
      both        = replicate(nboot, stat_from_resampled(resample_loci())))
    ci <- t(apply(bt, 1, stats::quantile, c(0.025, 0.975), na.rm = TRUE))
  } else {
    bt <- NULL
    ci <- matrix(NA_real_, length(point), 2, dimnames = list(names(point), NULL))
  }

  ## ---------------------------------------------------------------------------
  ## Result tables. Nothing is printed here: printing the returned object
  ## (print.raddiv_diversity(), below) gives the full report -- which is also
  ## what the command-line script shows.
  ## ---------------------------------------------------------------------------
  gv  <- function(pfx) point[paste0(pfx, names(pops))]
  lo  <- function(pfx) ci[paste0(pfx, names(pops)), 1]
  hi  <- function(pfx) ci[paste0(pfx, names(pops)), 2]
  sej <- function(pfx) se_jk[paste0(pfx, names(pops))]
  tab <- data.frame(population = names(pops), n = as.integer(nmax),
                    Ho = round(gv("Ho_"), 4), Ho_se = round(sej("Ho_"), 4),
                    Ho_lo = round(lo("Ho_"), 4), Ho_hi = round(hi("Ho_"), 4),
                    He = round(gv("He_"), 4), He_se = round(sej("He_"), 4),
                    He_lo = round(lo("He_"), 4), He_hi = round(hi("He_"), 4),
                    Fis = round(gv("Fis_"), 4), Fis_se = round(sej("Fis_"), 4),
                    Fis_lo = round(lo("Fis_"), 4), Fis_hi = round(hi("Fis_"), 4),
                    pct_poly = round(gv("poly_"), 1), row.names = NULL)
  ## priv_total: the dataset-wide sum of privAr over whichever loci had a
  ## defined value for that population (see stat_mat() above). NA means no
  ## locus was defined for that population at all, not "zero private
  ## alleles" -- those are different claims.
  rich <- data.frame(population = names(pops),
                     Ar = round(gv("Ar_"), 4), Ar_se = round(sej("Ar_"), 4), Ar_n = ar_n,
                     Ar_lo = round(lo("Ar_"), 4), Ar_hi = round(hi("Ar_"), 4),
                     privAr = round(gv("Pr_"), 4), privAr_se = round(sej("Pr_"), 4), privAr_n = pr_n,
                     privAr_lo = round(lo("Pr_"), 4), privAr_hi = round(hi("Pr_"), 4),
                     priv_total = round(gv("PrTot_"), 1), priv_total_se = round(sej("PrTot_"), 1),
                     priv_total_lo = round(lo("PrTot_"), 1),
                     priv_total_hi = round(hi("PrTot_"), 1), row.names = NULL)
  ## Individual-jackknife SEs, each right after its locus-based counterpart.
  if (!is.null(se_ind)) {
    put_after <- function(d, after, name, v) {
      d[[name]] <- round(unname(v), 4)
      k <- match(after, names(d))
      d[c(names(d)[seq_len(k)], name, setdiff(names(d)[-seq_len(k)], name))]
    }
    sei <- function(pfx) se_ind[paste0(pfx, names(pops))]
    tab  <- put_after(tab, "Ho_se", "Ho_se_ind", sei("Ho_"))
    tab  <- put_after(tab, "He_se", "He_se_ind", sei("He_"))
    tab  <- put_after(tab, "Fis_se", "Fis_se_ind", sei("Fis_"))
    rich <- put_after(rich, "Ar_se", "Ar_se_ind", sei("Ar_"))
    rich <- put_after(rich, "privAr_se", "privAr_se_ind", sei("Pr_"))
  }
  ## The same parameter with the estimator behind Stacks' `Pi`, side by side.
  ## He can legitimately be 0 (a population monomorphic at every retained
  ## locus): print "n/a" then, rather than "NaN%"/"Inf%".
  he_ <- gv("He_"); he2n_ <- gv("He2n_")
  pct_diff <- ifelse(is.na(he_) | is.na(he2n_), NA_character_,
                     ifelse(he_ > 0, sprintf("%+.2f%%", 100 * (he2n_ / he_ - 1)),
                            "n/a (He=0)"))
  cmp <- data.frame(population = names(pops),
                    He_NeiChesser = round(gv("He_"), 4),
                    He_2n_corr    = round(gv("He2n_"), 4),
                    pct_diff = pct_diff,
                    Fis_NeiChesser = round(gv("Fis_"), 4),
                    Fis_2n_corr    = round(gv("Fis2n_"), 4),
                    row.names = NULL)
  ## Per-sequenced-site values. `ok_sites`: is each population's own sites
  ## value a plausible sequenced-site count (finite and >= n_rec, the variant
  ## records called across the whole dataset)? Only on a SNP VCF, where one
  ## record is one site; autosomal_het() itself gives NA wherever !ok_sites.
  ok_sites <- is.finite(sites_vec) & sites_vec >= n_rec
  aut <- NULL
  if (any(ok_sites) && !is_hapvcf)
    aut <- data.frame(population = names(pops),
      sites_used = ifelse(ok_sites, format(sites_vec, big.mark = ",", scientific = FALSE), NA),
      Ho_autosomal = signif(autosomal_het(gv("Ho_"), n_rec, sites_vec), 4),
      He_autosomal = signif(autosomal_het(gv("He_"), n_rec, sites_vec), 4),
      row.names = NULL)
  ## Two populations: the loci-bootstrap interval on their He difference, for
  ## the report (which also says why it is not the test to quote).
  he_diff <- NULL
  if (nboot > 0 && r == 2) {
    d <- bt[paste0("He_", names(pops)[1]), ] - bt[paste0("He_", names(pops)[2]), ]
    he_diff <- list(pops = names(pops), est = he_[[1]] - he_[[2]],
                    lo = stats::quantile(d, 0.025, names = FALSE),
                    hi = stats::quantile(d, 0.975, names = FALSE))
  }

  ## Files, only when asked for (`outdir`). The workflow runs this function
  ## TWICE, on the haplotype VCF and the SNP VCF, so the input's stem is in
  ## every file name (populations.haps.vcf -> diversity_per_population.haps.tsv)
  ## and the second run cannot overwrite the first; a non-default boot mode is
  ## in the name too, so a comparison run cannot overwrite the default output.
  written <- character(0)
  if (!is.null(outdir)) {
    dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
    stem <- .derive_stem(vcf_file, stem)
    boot_suffix <- if (boot_mode != "loci") paste0(".boot-", boot_mode) else ""
    out_tabs <- list(per_population = tab, richness = rich, autosomal = aut)
    for (nm in names(out_tabs)) {
      if (is.null(out_tabs[[nm]])) next
      f <- file.path(outdir, sprintf("diversity_%s.%s%s.tsv", nm, stem, boot_suffix))
      utils::write.table(out_tabs[[nm]], f, sep = "\t", quote = FALSE, row.names = FALSE)
      written <- c(written, f)
    }
  }

  structure(
    list(per_population = tab, richness = rich, autosomal = aut),
    class = "raddiv_diversity",
    report = list(n_pop = r, n_used = L, n_loci = nL, n_rec = n_rec, g = g,
                  complete_case = complete_case, min_n = min_n,
                  cells_used = cells_used, cells_total = cells_total,
                  nboot = nboot, boot = boot_mode, n_ind = nmax, pr_n = pr_n,
                  is_haplotype = is_hapvcf, sites = sites_vec, ok_sites = ok_sites,
                  estimators = cmp, he_diff = he_diff, se_individuals = !is.null(se_ind),
                  files = written))
}

#' @rdname diversity_stats
#' @param x A `raddiv_diversity` object, as returned by `diversity_stats()`.
#' @param ... Ignored.
#' @export
print.raddiv_diversity <- function(x, ...) {
  rp <- attr(x, "report")
  r <- rp$n_pop; L <- rp$n_used; n_rec <- rp$n_rec; g <- rp$g
  ## Printed tables must not line-wrap: this output is meant to be greppable
  ## and diffable, and 80 columns wraps the richness table. Restored on exit.
  old <- options(width = max(200L, getOption("width")))
  on.exit(options(old), add = TRUE)
  ## An indented paragraph, and a pointer to the README section that explains
  ## it: the report prints the numbers, the README holds the long reasoning.
  note <- function(...) cat(paste0("  ", c(...), "\n"), sep = "")
  see  <- function(sec) cat(paste0("  -> README.md, \"", sec, "\"\n"))

  cat("\n=====================================================================\n")
  cat("  DIVERSITY --", r, "populations,", format(L, big.mark = ","), "records on",
      format(rp$n_loci, big.mark = ","), "RAD loci\n")
  cat("  He = Nei & Chesser (1983) Hs, the hierfstat quantity\n")
  cat("  FIS = 1 - sum(Ho)/sum(Hs), a ratio of sums, monomorphic loci included\n")
  cat("  Ar / private Ar rarefied to", g, "GENE COPIES (Kalinowski 2004)\n")
  if (rp$complete_case) {
    cat(sprintf("  %s of %s variant records have complete data and were used (%.1f%%)\n",
                format(L, big.mark = ","), format(n_rec, big.mark = ","), 100 * L / n_rec))
    if (L / n_rec < 0.5)
      cat("  ^ most records were dropped for missing data. complete_case keeps a\n    record only if genotyped in EVERY individual (Schmidt et al. 2021), so a\n    few poor libraries can cost most of the dataset. Check the per-individual\n    table from het_between_pops() (loci_called column) before accepting\n    this -- or drop complete_case.\n")
  } else {
    cat(sprintf("  %s of %s locus-by-population cells meet min_n = %d (%.1f%%)\n",
                format(rp$cells_used, big.mark = ","), format(rp$cells_total, big.mark = ","),
                rp$min_n, 100 * rp$cells_used / rp$cells_total))
    cat("  Each population's Ho/He/Fis is a ratio of sums over whichever loci it\n    cleared min_n at -- see the per-population lines printed while it ran.\n")
  }
  if (rp$nboot > 0) cat("  95% CI from ", format(rp$nboot, big.mark = ","),
                        " replicates, boot=\"", rp$boot,
                        "\" (see README.md, \"Bootstrap mode\")\n", sep = "")
  cat("=====================================================================\n\n")
  print(x$per_population, row.names = FALSE)
  cat("  _se columns: delete-one-block jackknife over RAD loci (Weir 1996) --\n")
  cat("  always available (no nboot needed), on the same loci axis as _lo/_hi.\n")
  if (isTRUE(rp$se_individuals)) {
    cat("  _se_ind columns: delete-one-INDIVIDUAL jackknife -- the uncertainty from\n")
    cat("  which individuals were sampled, which the loci-based _se/_lo/_hi hold fixed.\n")
    cat("  Report _se_ind when individuals differ in inbreeding (g2 > 0; see\n")
    cat("  identity_disequilibrium() or the het_between_pops() report).\n")
  }
  cat("\n")
  print(x$richness, row.names = FALSE)
  cat("  _se columns: delete-one-block jackknife over RAD loci, same as above.\n")
  cat("  Ar_n / privAr_n: how many of the", format(L, big.mark = ","),
      "loci had a defined value for that population\n")
  cat("  (Ar needs only THAT population at >= g gene copies; privAr needs EVERY\n")
  cat("  population at >= g simultaneously, so privAr_n <= Ar_n, often far less\n")
  cat("  under available data -- see README.md if privAr_n is small).\n")
  if (min(rp$pr_n) < 0.5 * L)
    cat(sprintf("  WARNING: privAr_n is below 50%% of loci for at least one population --\n  its rarefied private richness rests on a minority of loci. Consider a\n  smaller g (rarefaction target), currently %d gene copies.\n", g))
  cat("  Ar and privAr are PER LOCUS (the HP-RARE / ADZE convention);\n")
  cat("  priv_total is the dataset-wide sum over privAr_n loci. Private = absent from all", r - 1, "others,\n")
  cat("  so it is not comparable to a dataset with a different number of populations.\n")
  arv <- x$richness$Ar
  if (all(is.na(arv))) {
    cat(sprintf("  NOTE: Ar is undefined for every population -- no locus has any\n  population at >= g = %d gene copies (see Ar_n above). That is too little\n  data at this rarefaction target, not a biallelic-SNP ceiling; try a smaller g.\n", g))
  } else if (max(arv, na.rm = TRUE) <= 2.001) {
    cat("  NOTE: Ar is capped at 2 because these are biallelic SNPs, so rarefaction\n  has almost nothing to correct. Run this on the HAPLOTYPE VCF for a\n  richness worth reporting; private allelic richness is still informative.\n")
  }

  ## Side by side with the estimator behind Stacks' `Pi` column
  cat("\nTwo estimators of the same parameter (gene diversity per variant record)\n")
  note("He_NeiChesser  (n/(n-1))(1 - sum p^2 - Ho/2n)   unbiased at any FIS",
       "He_2n_corr     (2n/(2n-1))(1 - sum p^2)         unbiased only if FIS = 0",
       "               -- this is the estimator behind Stacks' `Pi` column.")
  print(rp$estimators, row.names = FALSE)
  nmax <- rp$n_ind
  note("He barely moves; FIS moves ~(1-F)/F times as much, because FIS divides by",
       "the small quantity the two estimators disagree about. So Stacks' Pi is safe",
       "to quote and its Fis is not (and its Fis carries a second, larger error:",
       "it is a mean of per-locus ratios).",
       sprintf("Compare against `Pi`, never `Exp_Het`: Pi = Exp_Het * 2n/(2n-1) = %s.",
               paste(sprintf("x%.3f at n=%d", (2 * nmax) / (2 * nmax - 1), nmax),
                     collapse = ", ")))
  see("Formulas, and what maps onto what in Stacks")

  cat("\nWhat the Ho and He numbers are\n")
  sites_vec <- rp$sites; ok_sites <- rp$ok_sites
  mode_txt <- if (rp$complete_case) "complete-case"
              else sprintf("available-data, min_n = %d", rp$min_n)
  if (any(ok_sites) && rp$is_haplotype) {
    note("`sites` was supplied but this is a HAPLOTYPE VCF, so the autosomal",
         "conversion was SKIPPED: He here is per-tag gene diversity, not per-site,",
         "and dividing by sequenced sites would understate nucleotide diversity.",
         "Re-run on populations.snps.vcf with the same `sites` for He_autosomal.")
    see("Which file for which statistic")
  } else if (any(ok_sites)) {
    retain <- L / n_rec
    if (length(unique(sites_vec[ok_sites])) <= 1L) {
      note(sprintf("Scaled to %s sequenced sites (same for every population) via all %s variant records (%.3f%% variant);",
                   format(sites_vec[ok_sites][1], big.mark = ",", scientific = FALSE),
                   format(n_rec, big.mark = ","), 100 * n_rec / sites_vec[ok_sites][1]),
           sprintf("He estimated using %s of %s records (%.1f%% of them, %s).",
                   format(L, big.mark = ","), format(n_rec, big.mark = ","), 100 * retain,
                   mode_txt))
    } else {
      note(sprintf("Scaled to EACH POPULATION'S OWN sequenced-site count (%s-%s; see `sites_used`",
                   format(min(sites_vec[ok_sites]), big.mark = ",", scientific = FALSE),
                   format(max(sites_vec[ok_sites]), big.mark = ",", scientific = FALSE)),
           sprintf("above) via all %s variant records (%.3f%% variant, dataset-wide);",
                   format(n_rec, big.mark = ","), 100 * n_rec / mean(sites_vec[ok_sites])),
           sprintf("He estimated using %s of %s records (%.1f%% of them, %s).",
                   format(L, big.mark = ","), format(n_rec, big.mark = ","), 100 * retain,
                   mode_txt))
    }
    if (any(!ok_sites)) {
      reasons <- vapply(names(sites_vec)[!ok_sites], function(p)
        if (sites_vec[p] > 0)
          sprintf("%-22s sites = %s, less than the %s variant records called",
                  p, format(sites_vec[p], big.mark = ","), format(n_rec, big.mark = ","))
        else sprintf("%-22s no sites value supplied", p),
        character(1))
      note("No autosomal value for:", reasons)
    }
    print(x$autosomal, row.names = FALSE)
    note("He_autosomal is the SAME per-site quantity as He above, averaged over every",
         "sequenced site instead of over variant records. That denominator is what",
         "makes it nucleotide diversity (pi). There is no column named `pi`.",
         "Compare against Stacks' `Pi` from the All-positions block -- not the",
         "variant-positions block, and not `Exp_Het` (Schmidt et al. 2021).")
    see("Formulas, and what maps onto what in Stacks")
    if (rp$complete_case && retain < 0.5)
      note(sprintf("NOTE: only %.0f%% of variant records have complete data, so this rests on",
                   100 * retain),
           "them being representative. If missingness is concentrated in a few poor",
           "libraries it is not -- drop complete_case (the default already uses",
           "available data per population) or drop those individuals and re-run.")
  } else {
    if (any(sites_vec > 0))
      note(sprintf("`sites` was supplied but every value is LESS than the %s variant records",
                   format(n_rec, big.mark = ",")),
           "this run called, which is not a plausible sequenced-site count -- it was",
           "IGNORED rather than used. Pass the `Sites` column of the All-positions",
           "block of populations.sumstats_summary.tsv, which should be orders of",
           "magnitude larger than the number of variant records.")
    else
      note("Per ASCERTAINED record: the denominator is however many markers this run",
           "called. Valid for comparing these populations to each other and nothing",
           "else. Pass `sites` (the `Sites` column of the All-positions block of",
           "populations.sumstats_summary.tsv) to get the autosomal value.")
    see("Denominators, and where nucleotide diversity fits")
  }

  ## What to take from THIS run
  cat("\n---------------------------------------------------------------------\n")
  if (rp$is_haplotype) {
    note("TAKE FROM THIS RUN:  Fis, Ar, privAr.",
         "NOT Ho / He: haplotype gene diversity has no fixed ceiling, so it moves",
         "with read length, enzyme and filters and is not comparable across",
         "studies. Run populations.snps.vcf for those. Fis is a ratio, so the",
         "scale cancels -- which is why it belongs here.")
  } else {
    note("TAKE FROM THIS RUN:  Ho, He, pct_poly, Ho_autosomal, He_autosomal.",
         "Run populations.haps.vcf for Fis, Ar and privAr: on biallelic SNPs",
         "rarefied richness is bounded at 2, and Fis is more precise from",
         "multi-allelic loci.")
  }
  note("",
       "Fis is comparable between the two files; Ho and He are NOT (haplotype Ho",
       "asks 'is this tag heterozygous', per-site Ho asks 'is this site'). Label",
       "every table row with the file it came from.",
       "",
       "He and Ho above are averaged over ALL retained records including those",
       "monomorphic within a population -- the POOLED denominator, matching",
       "Stacks' Variant-positions block. Report pct_poly beside them:",
       "    He(pooled) = He(polymorphic within pop) x pct_poly/100   (exactly)",
       "Do not switch to the within-population denominator; it conditions on the",
       "outcome and can reverse a real difference.")
  see("Which file for which statistic")
  cat("---------------------------------------------------------------------\n")

  cat("\nFor Weir & Cockerham FST, Jost's D and Weir & Goudet's beta between these\n")
  cat("populations, run differentiation_stats() on the same file.\n")

  hd <- rp$he_diff
  if (!is.null(hd)) {
    cat(sprintf("\nHe difference (%s - %s): %+.4f  95%% CI [%+.4f, %+.4f]\n",
                hd$pops[1], hd$pops[2], hd$est, hd$lo, hd$hi))
    if (rp$boot == "loci") {
      cat("  CAUTION: this interval treats LOCI as the replicate (boot=\"loci\"). For a\n")
      cat("  comparison of population MEAN heterozygosity that is pseudoreplication --\n")
      cat("  the uncertainty is dominated by which INDIVIDUALS you sampled. Use\n")
      cat("  het_between_pops() for the p-value you report.\n")
    } else {
      cat(sprintf("  This interval also resamples INDIVIDUALS (boot=\"%s\"), unlike a plain\n", rp$boot))
      cat("  locus bootstrap -- but it is still not a paired individual-level test. Use\n")
      cat("  het_between_pops() for the p-value you report.\n")
    }
  }
  if (length(rp$files)) cat(sprintf("\nWrote %s\n", paste(rp$files, collapse = ", ")))
  cat("\n")
  invisible(x)
}

## Not exported. Resolves diversity_stats()'s `sites` argument into a named
## numeric vector aligned to `pop_names`, in that order -- mirrors
## .resolve_H()'s "accept several input shapes, fail loudly on mismatch"
## style (R/vcf_io.R). Four accepted shapes:
##   - a single non-negative number: broadcast to every population (the
##     original, still-default behavior; 0 means "off")
##   - a named numeric vector, names = population names
##   - a data frame or matrix with a population-name column
##     ("population"/"pop", case-insensitive) and a sites column
##     ("sites"/"site", case-insensitive)
##   - a path to populations.sumstats_summary.tsv: reads it via
##     read_sumstats_summary() and uses `all_positions$sites`, named by
##     `all_positions$population`
##
## Called twice from diversity_stats(): once with `pop_names = NULL`, right
## after argument parsing and before the (potentially slow) VCF read, purely
## to catch a malformed `sites` early (bad type, a negative number, a
## missing file) -- population names aren't known yet at that point, so
## name-matching is skipped and the return value is unused. Called again
## with the real `pop_names`, right after the popmap is read, to do the
## actual per-population matching and return the vector diversity_stats()
## uses.
.resolve_sites <- function(sites, pop_names = NULL) {
  if (is.character(sites)) {
    if (length(sites) != 1L) stop("sites: a file path must be a single string.")
    if (!file.exists(sites)) stop("sites: file not found: ", sites)
    ap <- read_sumstats_summary(sites)$all_positions
    v <- stats::setNames(ap$sites, ap$population)
  } else if (is.data.frame(sites) || is.matrix(sites)) {
    df <- as.data.frame(sites, stringsAsFactors = FALSE)
    nmL <- tolower(names(df))
    pop_col  <- which(nmL %in% c("population", "pop"))
    site_col <- which(nmL %in% c("sites", "site"))
    if (length(pop_col) != 1L || length(site_col) != 1L)
      stop("sites: a data frame/matrix must have exactly one population ",
           "column (\"population\" or \"pop\") and exactly one sites column ",
           "(\"sites\"). Found columns: ", paste(names(df), collapse = ", "))
    v <- stats::setNames(suppressWarnings(as.numeric(df[[site_col]])),
                         as.character(df[[pop_col]]))
  } else if (is.numeric(sites)) {
    v <- sites
  } else {
    stop("sites must be a non-negative number, a named numeric vector or ",
         "data frame keyed by population, or a path to ",
         "populations.sumstats_summary.tsv (got class: ",
         paste(class(sites), collapse = "/"), ").")
  }

  if (anyNA(v) || any(v < 0))
    stop("sites: every value must be a non-negative, non-missing number.")
  is_scalar <- length(v) == 1L && is.null(names(v))
  if (!is_scalar && (is.null(names(v)) || any(!nzchar(names(v)))))
    stop("sites: a vector/data frame/file of more than one value must be ",
         "named/keyed by population throughout.")

  if (is.null(pop_names)) return(invisible(NULL))
  if (is_scalar) v <- stats::setNames(rep(as.numeric(v), length(pop_names)), pop_names)
  .match_sites_to_pops(v, pop_names)
}

## Not exported. Matches a named `sites` vector against `pop_names` in BOTH
## directions: every population must have a value (error, naming which are
## missing) and a name in `sites` that matches no current population is
## reported (message, not silently dropped -- most likely a typo).
.match_sites_to_pops <- function(v, pop_names) {
  missing_pop <- setdiff(pop_names, names(v))
  if (length(missing_pop))
    stop("sites: no value given for population(s): ",
         paste(missing_pop, collapse = ", "),
         ". Every population in the popmap needs a sites value.")
  extra <- setdiff(names(v), pop_names)
  if (length(extra))
    message("sites: ", length(extra), " name(s) not among this run's ",
            "populations, ignored: ", paste(extra, collapse = ", "))
  v[pop_names]
}
