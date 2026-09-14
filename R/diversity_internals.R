###############################################################################
#
#  R/diversity_internals.R -- the building blocks of diversity_stats().
#
#  diversity_stats() (R/diversity_stats.R) reads as an outline; the work is
#  done here, in this order:
#
#    .locus_inclusion()        which records each population may use
#    .record_stats()           Ho, He, Ar, privAr for every record and
#                              population (calls .one_pop_record_stats())
#    .diversity_pieces()       those per-record values as sums-and-counts
#    .diversity_from_sums()    summed pieces -> the reported statistics
#    .jack_individuals()       delete-one-individual standard errors
#    .boot_individuals()       the boot = "individuals"/"both" bootstrap
#    .diversity_tables()       the result tables
#    .hierfstat_crosscheck()   optional comparison with hierfstat
#
#  Every statistic is computed by ONE function (.one_pop_record_stats() and
#  .diversity_from_sums()), whether for the point estimate, a jackknife
#  replicate or a bootstrap replicate. The formulas themselves are in
#  R/estimators.R.
#
#  SHAPES used throughout:
#    counts    a named list, one element per population, each a
#              records x alleles integer matrix of gene-copy counts
#              (see .pop_counts() in R/vcf_io.R)
#    n_typed   records x populations: individuals genotyped
#    n_het     records x populations: heterozygous individuals
#
###############################################################################

## ---------------------------------------------------------------------------
## Which records are used
## ---------------------------------------------------------------------------

## Not exported. Decides which records enter the analysis.
##   n_typed    records x populations matrix of genotyped individuals
##   pop_sizes  individuals per population
## Default rule (complete_case = FALSE): a record is kept when at least one
## population has >= min_n genotyped individuals there. Which populations
## actually use it is decided later, cell by cell (see `min_n` in
## .record_stats()).
## complete_case = TRUE: a record is kept only if EVERY individual of EVERY
## population is genotyped there (Schmidt et al. 2021). One population's
## missing individual then costs every population that record.
## Returns list(keep = logical per record, cells_used, cells_total).
.locus_inclusion <- function(n_typed, pop_sizes, min_n, complete_case, verbose) {
  n_rec <- nrow(n_typed)
  n_pops <- ncol(n_typed)
  cells_used <- NA_real_
  cells_total <- NA_real_

  if (complete_case) {
    fully_typed <- sweep(n_typed, 2L, pop_sizes, FUN = "==")   # records x populations
    keep <- rowSums(fully_typed) == n_pops
    .inform(verbose, sprintf("Complete-data record set (complete_case): %s of %s records (%.1f%%)",
                             .big(sum(keep)), .big(n_rec), 100 * mean(keep)))
    if (sum(keep) < 0.05 * n_rec)
      .inform(verbose, "  WARNING: under 5% of records survive. Report a looser rule alongside this.")
    if (!any(keep))
      stop("No record is genotyped in every individual. With many populations this is ",
           "common; drop complete_case or drop the worst individuals.", call. = FALSE)
  } else {
    meets_min_n <- n_typed >= min_n
    keep <- rowSums(meets_min_n) > 0
    cells_used <- sum(meets_min_n)
    cells_total <- n_rec * n_pops
    .inform(verbose, sprintf("Available-data mode (min_n = %d): %s of %s records usable by >= 1 population",
                             min_n, .big(sum(keep)), .big(n_rec)))
    .inform(verbose, sprintf("  %s of %s record-by-population cells meet min_n (%.1f%%)",
                             .big(cells_used), .big(cells_total), 100 * cells_used / cells_total))
    for (p in colnames(n_typed)) {
      usable <- meets_min_n[, p]
      ## The typed-n summary is printed only when the population has a usable
      ## record: min() of nothing is Inf, which "%d" cannot print.
      typed_summary <- if (any(usable))
        sprintf(" (typed-n there: mean %.1f, min %d, max %d)", mean(n_typed[usable, p]),
                min(n_typed[usable, p]), max(n_typed[usable, p]))
      else ""
      .inform(verbose, sprintf("    %-22s %s of %s records usable%s", p, .big(sum(usable)),
                               .big(n_rec), typed_summary))
    }
    if (!any(keep))
      stop("No record has >= min_n typed individuals in any population. Lower min_n ",
           "(currently ", min_n, ") or check that genotypes were parsed.", call. = FALSE)
  }
  list(keep = keep, cells_used = cells_used, cells_total = cells_total)
}

## Not exported. Heterozygous individuals per record and population, for the
## records `rows`: a records x populations matrix.
.het_by_pop <- function(H, pops, rows) {
  out <- vapply(pops, function(ids)
    rowSums(H$A1[rows, ids, drop = FALSE] != H$A2[rows, ids, drop = FALSE], na.rm = TRUE),
    numeric(length(rows)))
  matrix(out, nrow = length(rows), dimnames = list(NULL, names(pops)))
}

## ---------------------------------------------------------------------------
## Per-record statistics
## ---------------------------------------------------------------------------

## Not exported. Ho, He (two estimators), Ar and privAr of ONE population at
## every record. All inputs are for the same records:
##   counts            records x alleles gene-copy counts of this population
##   n_typed, n_het    genotyped and heterozygous individuals per record
##   absent_elsewhere  records x alleles: the probability that NONE of the
##                     other populations' g-copy samples contains the allele
##                     (see .absent_elsewhere()); needed for privAr
##   min_n             fewer genotyped individuals than this -> NA everywhere
##   g                 rarefaction size in gene copies
##   p_drawn           Pr(allele appears in g copies), if already computed
## Returns a list of per-record vectors:
##   Ho         observed heterozygosity
##   Hs         Nei & Chesser (1983) gene diversity (hs_nei_chesser())
##   Hp         the 2n/(2n-1)-corrected gene diversity behind Stacks' `Pi`
##              (gene_div_2n_counts()); reported only for comparison
##   Ar, Pr     rarefied allelic and private allelic richness
##              (rare_richness(), rare_private())
##   n_present  number of distinct alleles present
.one_pop_record_stats <- function(counts, n_typed, n_het, absent_elsewhere, min_n, g,
                                  p_drawn = .p_sampled_mat(counts, g)) {
  n_copies <- rowSums(counts)
  sum_p2 <- rowSums((counts / n_copies)^2)     # sum of squared allele frequencies

  Ho <- n_het / n_typed
  Hs <- hs_nei_chesser(sum_p2, Ho, n_typed)
  Hp <- (n_copies / (n_copies - 1)) * (1 - sum_p2)
  Hs[!(n_copies >= 2)] <- NA_real_
  Hp[!(n_copies >= 2)] <- NA_real_

  ## Rarefaction: expected number of alleles in g copies (Ar), and expected
  ## number drawn here AND absent from every other population (Pr). A record
  ## with fewer than g copies in any population has NA here (see
  ## .p_sampled_mat()).
  Ar <- rowSums(p_drawn)
  Pr <- rowSums(p_drawn * absent_elsewhere)

  too_few <- n_typed < min_n
  Ho[too_few] <- NA_real_
  Hs[too_few] <- NA_real_
  Hp[too_few] <- NA_real_
  Ar[too_few] <- NA_real_
  Pr[too_few] <- NA_real_
  list(Ho = Ho, Hs = Hs, Hp = Hp, Ar = Ar, Pr = Pr, n_present = rowSums(counts > 0))
}

## Not exported. For population `p`: records x alleles, the probability that
## the allele is absent from the g-copy sample of EVERY other population,
## i.e. the product of (1 - Pr(drawn)) over the other populations.
.absent_elsewhere <- function(p_drawn, p) {
  Reduce(`*`, lapply(p_drawn[-p], function(P) 1 - P))
}

## Not exported. .one_pop_record_stats() for every population, returned as
## records x populations matrices (Ho, Hs, Hp, Ar, Pr, n_present).
## `zero_low_n_counts = TRUE` sets a population's counts to zero wherever it
## has fewer than `min_n` genotyped individuals BEFORE rarefying. The
## boot = "individuals"/"both" bootstrap has always done this; the point
## estimate has not. The two differ only for privAr, and only at a record
## where one population is below min_n yet still has >= g gene copies (there,
## the point estimate keeps the other populations' privAr and the individual
## bootstrap does not).
.record_stats <- function(counts, n_typed, n_het, min_n, g, zero_low_n_counts = FALSE) {
  pop_names <- names(counts)
  if (zero_low_n_counts) {
    for (p in seq_along(counts)) counts[[p]][n_typed[, p] < min_n, ] <- 0L
  }
  p_drawn <- lapply(counts, .p_sampled_mat, g = g)
  per_pop <- lapply(seq_along(counts), function(p)
    .one_pop_record_stats(counts[[p]], n_typed[, p], n_het[, p],
                          .absent_elsewhere(p_drawn, p), min_n, g, p_drawn[[p]]))
  as_matrix <- function(name) {
    matrix(vapply(per_pop, `[[`, numeric(nrow(n_typed)), name),
           ncol = length(counts), dimnames = list(NULL, pop_names))
  }
  list(Ho = as_matrix("Ho"), Hs = as_matrix("Hs"), Hp = as_matrix("Hp"),
       Ar = as_matrix("Ar"), Pr = as_matrix("Pr"), n_present = as_matrix("n_present"))
}

## ---------------------------------------------------------------------------
## From per-record values to reported statistics
## ---------------------------------------------------------------------------

## Not exported. Every reported statistic is a ratio (or a plain sum) of
## COLUMN SUMS over records, so each record is turned into these 11 pieces
## per population:
##
##   piece        value at a record                  used for
##   n            1 if Ho and He are both defined    denominator of Ho, He, pct_poly
##   ho, hs       Ho, He (0 unless n = 1)            Ho, He, FIS
##   n_2n         1 if Ho, He and Hp are defined     denominator of He_2n
##   ho_2n, hp    Ho, Hp on that mask                He_2n, FIS_2n
##   n_pr, pr     1 / privAr where defined           privAr, priv_total
##   n_ar, ar     1 / Ar where defined               Ar
##   polymorphic  1 if >1 allele present (and n = 1) pct_poly
##
## An undefined value contributes 0 to its sum AND 0 to its count, so a mean
## always divides by the number of records that actually had a value. (Summing
## with na.rm = TRUE but dividing by the number of records would silently
## count an undefined record as a measured zero.)
##
## `rs` is the list from .record_stats() (records x populations matrices) or
## from .one_pop_record_stats() (vectors: one population). The result has one
## row per record and 11 * (number of populations) columns, grouped by piece:
## columns 1..r are piece 1 for populations 1..r, then piece 2, and so on.
.diversity_pieces <- function(rs) {
  usable <- is.finite(rs$Ho) & is.finite(rs$Hs)
  usable_2n <- usable & is.finite(rs$Hp)
  has_pr <- !is.na(rs$Pr)
  has_ar <- !is.na(rs$Ar)
  zero_unless <- function(x, ok) ifelse(ok, x, 0)
  pieces <- cbind(usable, zero_unless(rs$Ho, usable), zero_unless(rs$Hs, usable),
                  usable_2n, zero_unless(rs$Ho, usable_2n), zero_unless(rs$Hp, usable_2n),
                  has_pr, zero_unless(rs$Pr, has_pr),
                  has_ar, zero_unless(rs$Ar, has_ar),
                  zero_unless(rs$n_present > 1, usable))
  storage.mode(pieces) <- "double"
  pieces
}

## Not exported. Summed pieces -> statistics.
##   totals     one row per estimate (point estimate, jackknife or bootstrap
##              replicate), columns laid out as in .diversity_pieces()
##   pop_names  population names, in column order
## Returns one row per estimate and these columns for each population:
##   Ho_, He_, Fis_          Nei & Chesser He; FIS = 1 - sum(Ho)/sum(He)
##   He2n_, Fis2n_           the same with Stacks' Pi estimator
##   Ar_, Pr_                mean rarefied allelic / private allelic richness
##   poly_                   % of records polymorphic in the population
##   PrTot_                  privAr summed over records (priv_total)
## A mean with no contributing record is NA ("no information"), never the NaN
## that 0/0 would give. A population with no record where privAr is defined
## gets priv_total NA, which is a different claim from zero private alleles.
.diversity_from_sums <- function(totals, pop_names) {
  n_pops <- length(pop_names)
  totals <- matrix(totals, ncol = 11L * n_pops)
  piece <- function(k) totals[, (k - 1L) * n_pops + seq_len(n_pops), drop = FALSE]
  n <- piece(1)
  ho <- piece(2)
  hs <- piece(3)
  n_2n <- piece(4)
  ho_2n <- piece(5)
  hp <- piece(6)
  n_pr <- piece(7)
  pr <- piece(8)
  n_ar <- piece(9)
  ar <- piece(10)
  polymorphic <- piece(11)

  mean_of <- function(total, count) {
    out <- total / count
    out[count < 0.5] <- NA_real_
    out
  }
  fis_of <- function(ho_total, he_total) {
    ifelse(he_total > .zero_tol, 1 - ho_total / he_total, NA_real_)
  }
  private_total <- pr
  private_total[n_pr < 0.5] <- NA_real_

  out <- cbind(mean_of(ho, n), mean_of(hs, n), fis_of(ho, hs),
               mean_of(hp, n_2n), fis_of(ho_2n, hp),
               mean_of(ar, n_ar), mean_of(pr, n_pr),
               100 * mean_of(polymorphic, n), private_total)
  colnames(out) <- paste0(rep(c("Ho_", "He_", "Fis_", "He2n_", "Fis2n_", "Ar_", "Pr_",
                                "poly_", "PrTot_"), each = n_pops), pop_names)
  out
}

## ---------------------------------------------------------------------------
## Diagnostic: FIS by call rate (null alleles and allele dropout)
## ---------------------------------------------------------------------------

## Not exported. For each population, the records it uses are split by that
## population's own call rate there (100%, 90-99%, 75-89%, below 75%), and
## FIS (ratio of sums) and mean He are computed within each group, with a
## delete-one-locus jackknife SE.
##
## WHY. A null allele (e.g. a mutation in the restriction site) or allele
## dropout at low depth turns heterozygotes into apparent homozygotes, and
## when both copies drop out, into missing genotypes. Records affected this
## way have both more missing data and fewer heterozygotes, so FIS rising as
## call rate falls is their signature. Real inbreeding raises FIS at every
## record alike. Dropout is also more common in more diverse populations
## (Gautier et al. 2013), so it can create a FIS difference between
## populations as well as inflate FIS.
##   rs           per-record statistics from .record_stats() (records x pops)
##   n_typed      records x populations, genotyped individuals
##   pop_sizes    individuals per population
##   locus_index  RAD locus of each record, as 1..n_loci
## Returns a data frame: population, call_rate, n_records, He, Fis, Fis_se;
## only groups with at least one record are listed.
.fis_by_call_rate <- function(rs, n_typed, pop_sizes, locus_index) {
  pop_names <- colnames(n_typed)
  breaks <- c("100%" = 1, "90-99%" = 0.9, "75-89%" = 0.75, "<75%" = 0)
  one_pop <- function(p) {
    call_rate <- n_typed[, p] / pop_sizes[[p]]
    group <- ifelse(call_rate >= 1 - .threshold_tol, "100%",
             ifelse(call_rate >= 0.9 - .threshold_tol, "90-99%",
             ifelse(call_rate >= 0.75 - .threshold_tol, "75-89%", "<75%")))
    pieces <- .diversity_pieces(lapply(rs, function(m) m[, p]))
    rows <- lapply(names(breaks), function(g) {
      in_group <- group == g
      n_records <- sum(in_group & pieces[, 1] > 0)
      if (!n_records) return(NULL)
      S <- rowsum(pieces * in_group, locus_index, reorder = TRUE)
      f <- function(totals) .diversity_from_sums(totals, "x")
      point <- f(colSums(S))[1L, ]
      se <- .jack_block_sums(S, f)
      data.frame(population = p, call_rate = g, n_records = n_records,
                 He = unname(point["He_x"]), Fis = unname(point["Fis_x"]),
                 Fis_se = unname(se["Fis_x"]))
    })
    do.call(rbind, rows)
  }
  out <- do.call(rbind, lapply(pop_names, one_pop))
  rownames(out) <- NULL
  out
}

## ---------------------------------------------------------------------------
## Uncertainty over individuals
## ---------------------------------------------------------------------------

## Not exported. Delete-one-INDIVIDUAL jackknife standard errors of each
## population's Ho, He, FIS, Ar and privAr (diversity_stats(se_individuals =
## TRUE)): the uncertainty from which individuals were sampled, which the
## locus jackknife and bootstrap hold fixed. For each individual of
## population p, its two gene copies and its heterozygosity are removed from
## p's per-record counts, and p's statistics are recomputed with the same
## functions as the point estimate. The other populations do not change.
## Nobody is duplicated, so this avoids the small-sample bias that rules out
## a bootstrap over individuals (vignette("rationale"), "Why not bootstrap individuals?").
##   SE = sqrt( (n - 1)/n * sum_i (theta_(-i) - mean_i theta_(-i))^2 )
## Inputs are for the records `rows` of H (see the SHAPES note at the top).
## Returns list(se = named SEs, full = the same statistics recomputed from
## everyone, which must equal the point estimates).
.jack_individuals <- function(H, pops, rows, counts, n_typed, n_het, min_n, g) {
  p_drawn <- lapply(counts, .p_sampled_mat, g = g)
  reported <- c("Ho_", "He_", "Fis_", "Ar_", "Pr_")
  se <- full <- numeric(0)

  for (p in seq_along(pops)) {
    absent <- .absent_elsewhere(p_drawn, p)
    ## Population p's five statistics from its (possibly reduced) counts.
    stats_of <- function(C, n, h) {
      rs <- .one_pop_record_stats(C, n, h, absent, min_n, g)
      stats <- .diversity_from_sums(colSums(.diversity_pieces(rs)), "x")
      stats[1L, paste0(reported, "x")]
    }
    C <- counts[[p]]
    replicates <- vapply(pops[[p]], function(id) {
      a <- H$A1[rows, id]
      b <- H$A2[rows, id]
      typed <- !is.na(a)
      typed_rows <- which(typed)
      C_without <- C
      ## Two separate steps, so a homozygote loses both of its copies (a
      ## single assignment with a repeated index would subtract only once).
      C_without[cbind(typed_rows, a[typed_rows])] <- C_without[cbind(typed_rows, a[typed_rows])] - 1L
      C_without[cbind(typed_rows, b[typed_rows])] <- C_without[cbind(typed_rows, b[typed_rows])] - 1L
      stats_of(C_without, n_typed[, p] - typed, n_het[, p] - (typed & a != b))
    }, numeric(length(reported)))
    n_ind <- ncol(replicates)
    names_p <- paste0(reported, names(pops)[p])
    deviations <- replicates - rowMeans(replicates)
    se <- c(se, stats::setNames(sqrt((n_ind - 1) / n_ind * rowSums(deviations^2)), names_p))
    full <- c(full, stats::setNames(stats_of(C, n_typed[, p], n_het[, p]), names_p))
  }
  list(se = se, full = full)
}

## Not exported. The boot = "individuals" / "both" bootstrap (comparison
## modes only; see ?diversity_stats). Each replicate:
##   1. resamples individuals with replacement within each population (a
##      drawn-twice individual is simply included twice);
##   2. recomputes every record's statistics from those individuals, with the
##      same functions as the point estimate;
##   3. for "both", also resamples RAD loci (as boot = "loci" does).
## Why this is not the default: a duplicated individual biases He, FIS, Ar and
## privAr low, because their small-sample corrections assume n DISTINCT
## individuals (vignette("rationale"), "Why not bootstrap individuals?").
## Returns an nboot x (statistics) matrix.
.boot_individuals <- function(H, pops, rows, locus_index, n_loci, min_n, g, nboot,
                              resample_loci) {
  pop_names <- names(pops)
  k <- max(1L, H$n_alleles[rows])
  one_replicate <- function() {
    resampled <- lapply(pops, function(ids) {
      times_drawn <- tabulate(sample.int(length(ids), length(ids), replace = TRUE),
                              nbins = length(ids))
      rep(ids, times_drawn)
    })
    n_typed <- .typed_by_pop(H, resampled, rows)
    rs <- .record_stats(.pop_counts(H, resampled, rows, k), n_typed,
                        .het_by_pop(H, resampled, rows), min_n, g, zero_low_n_counts = TRUE)
    S <- rowsum(.diversity_pieces(rs), locus_index, reorder = TRUE)
    totals <- if (resample_loci) {
      times_drawn <- tabulate(sample.int(n_loci, n_loci, replace = TRUE), nbins = n_loci)
      crossprod(times_drawn, S)
    } else {
      colSums(S)
    }
    .diversity_from_sums(totals, pop_names)[1L, ]
  }
  t(vapply(seq_len(nboot), function(i) one_replicate(), numeric(9L * length(pops))))
}

## ---------------------------------------------------------------------------
## Result tables
## ---------------------------------------------------------------------------

## Not exported. How diversity_stats() tables are rounded when printed or
## written to a file (the returned object keeps full precision).
.diversity_rounding <- list(
  per_population = list(digits = 4, round_cols = c(pct_poly = 1)),
  richness = list(digits = 4, round_cols = c(priv_total = 1, priv_total_se = 1,
                                             priv_total_lo = 1, priv_total_hi = 1)),
  autosomal = list(digits = 4, signif_cols = c(Ho_autosomal = 4, He_autosomal = 4)),
  estimator_comparison = list(digits = 4, round_cols = c(pct_diff = 2)),
  he_difference = list(digits = 4),
  fis_by_call_rate = list(digits = 4)
)

## Not exported. Applies .diversity_rounding to one table, for printing or
## writing. `sites_used` and `variant_records` become plain digits, because R
## would otherwise print and write 100000 as "1e+05".
.round_diversity_table <- function(df, table_name) {
  if (is.null(df)) return(NULL)
  rules <- .diversity_rounding[[table_name]]
  df <- .round_table(df, digits = rules$digits, round_cols = rules$round_cols,
                     signif_cols = rules$signif_cols)
  for (count_col in intersect(c("sites_used", "variant_records"), names(df)))
    df[[count_col]] <- ifelse(is.na(df[[count_col]]), NA_character_,
                              format(df[[count_col]], scientific = FALSE, trim = TRUE))
  df
}

## Not exported. Builds diversity_stats()' result tables, full precision.
##   point, se_loci  named vectors from .diversity_from_sums() / jackknife
##   ci              statistics x 2 matrix (bootstrap interval; NA if none)
##   se_ind          individual-jackknife SEs, or NULL
##   ar_n, pr_n      records with a defined Ar / privAr, per population
##   sites, variant_records, ok_sites
##                   per population: sequenced sites, variant records with data
##                   there, and whether the sites value is usable
.diversity_tables <- function(pop_names, pop_sizes, point, se_loci, ci, se_ind,
                              ar_n, pr_n, sites, ok_sites, variant_records, is_haplotype,
                              boot_matrix) {
  value <- function(prefix) unname(point[paste0(prefix, pop_names)])
  se    <- function(prefix) unname(se_loci[paste0(prefix, pop_names)])
  lo    <- function(prefix) unname(ci[paste0(prefix, pop_names), 1])
  hi    <- function(prefix) unname(ci[paste0(prefix, pop_names), 2])

  per_population <- data.frame(
    population = pop_names, n = as.integer(pop_sizes),
    Ho = value("Ho_"), Ho_se = se("Ho_"), Ho_lo = lo("Ho_"), Ho_hi = hi("Ho_"),
    He = value("He_"), He_se = se("He_"), He_lo = lo("He_"), He_hi = hi("He_"),
    Fis = value("Fis_"), Fis_se = se("Fis_"), Fis_lo = lo("Fis_"), Fis_hi = hi("Fis_"),
    pct_poly = value("poly_"), row.names = NULL)

  richness <- data.frame(
    population = pop_names,
    Ar = value("Ar_"), Ar_se = se("Ar_"), Ar_n = as.integer(ar_n),
    Ar_lo = lo("Ar_"), Ar_hi = hi("Ar_"),
    privAr = value("Pr_"), privAr_se = se("Pr_"), privAr_n = as.integer(pr_n),
    privAr_lo = lo("Pr_"), privAr_hi = hi("Pr_"),
    priv_total = value("PrTot_"), priv_total_se = se("PrTot_"),
    priv_total_lo = lo("PrTot_"), priv_total_hi = hi("PrTot_"), row.names = NULL)

  ## Individual-jackknife SEs go right after their locus-based counterparts.
  if (!is.null(se_ind)) {
    insert_after <- function(df, after, name, values) {
      df[[name]] <- unname(values)
      position <- match(after, names(df))
      df[c(names(df)[seq_len(position)], name,
           setdiff(names(df)[-seq_len(position)], name))]
    }
    se_i <- function(prefix) se_ind[paste0(prefix, pop_names)]
    per_population <- insert_after(per_population, "Ho_se", "Ho_se_ind", se_i("Ho_"))
    per_population <- insert_after(per_population, "He_se", "He_se_ind", se_i("He_"))
    per_population <- insert_after(per_population, "Fis_se", "Fis_se_ind", se_i("Fis_"))
    richness <- insert_after(richness, "Ar_se", "Ar_se_ind", se_i("Ar_"))
    richness <- insert_after(richness, "privAr_se", "privAr_se_ind", se_i("Pr_"))
  }

  ## The same parameter with the estimator behind Stacks' `Pi`, side by side.
  he <- value("He_")
  he_2n <- value("He2n_")
  estimator_comparison <- data.frame(
    population = pop_names,
    He_NeiChesser = he, He_2n_corr = he_2n,
    pct_diff = ifelse(is.finite(he) & he > 0, 100 * (he_2n / he - 1), NA_real_),
    Fis_NeiChesser = value("Fis_"), Fis_2n_corr = value("Fis2n_"),
    row.names = NULL)

  ## Per sequenced site: only on a SNP VCF (one record = one site), and only
  ## for populations with a plausible sites value (see ok_sites). Each
  ## population uses its OWN sequenced sites and its OWN count of variant
  ## records (see .autosomal_counts() in R/diversity_stats.R).
  autosomal <- NULL
  if (any(ok_sites) && !is_haplotype) {
    records <- unname(variant_records)
    sites_ok <- ifelse(ok_sites, unname(as.vector(sites)), NA_real_)
    autosomal <- data.frame(
      population = pop_names,
      sites_used = sites_ok,
      variant_records = ifelse(ok_sites, records, NA_real_),
      Ho_autosomal = autosomal_het(value("Ho_"), records, sites_ok),
      He_autosomal = autosomal_het(value("He_"), records, sites_ok),
      row.names = NULL)
  }

  ## Two populations: the locus-bootstrap interval of their He difference.
  he_difference <- NULL
  if (!is.null(boot_matrix) && length(pop_names) == 2L) {
    d <- boot_matrix[, paste0("He_", pop_names[1])] - boot_matrix[, paste0("He_", pop_names[2])]
    he_difference <- data.frame(pop1 = pop_names[1], pop2 = pop_names[2],
                                He_diff = he[1] - he[2],
                                lo = stats::quantile(d, 0.025, names = FALSE, na.rm = TRUE),
                                hi = stats::quantile(d, 0.975, names = FALSE, na.rm = TRUE))
  }

  list(per_population = per_population, richness = richness, autosomal = autosomal,
       estimator_comparison = estimator_comparison, he_difference = he_difference)
}

## ---------------------------------------------------------------------------
## Optional cross-check against hierfstat
## ---------------------------------------------------------------------------

## Not exported. diversity_stats(hierfstat_check = TRUE): computes per-record
## Ho/Hs with hierfstat::basic.stats() and allelic richness with
## hierfstat::allelic.richness(), and reports the largest difference from
## this package's own values. Every reported number still comes from this
## package. A failed or incomparable hierfstat call is reported the same way
## for both functions: a warning saying the check was skipped.
.hierfstat_crosscheck <- function(H, pops, rows, rs, n_typed, g, min_n, verbose) {
  dat <- .to_hierfstat_df(H, pops, rows)
  compared <- n_typed >= min_n          # the cells this package reports

  bs <- tryCatch(hierfstat::basic.stats(dat, diploid = TRUE, digits = 12),
                 error = function(e) e)
  if (inherits(bs, "error")) {
    warning("Skipped the hierfstat::basic.stats() cross-check: it failed with \"",
            conditionMessage(bs), "\".", call. = FALSE)
  } else {
    d_ho <- max(abs(as.matrix(bs$Ho)[compared] - rs$Ho[compared]), na.rm = TRUE)
    d_hs <- max(abs(as.matrix(bs$Hs)[compared] - rs$Hs[compared]), na.rm = TRUE)
    .inform(verbose, sprintf("  cross-check vs hierfstat::basic.stats(): max diff Ho %.2e, Hs %.2e",
                             d_ho, d_hs))
    if (max(d_ho, d_hs) > 1e-8) warning("Ho/Hs disagree with hierfstat::basic.stats().", call. = FALSE)
  }

  ## hierfstat adds a "dummy.loc" row when given exactly one locus, so the
  ## shapes then differ and the values cannot be compared.
  ar_hf <- tryCatch(hierfstat::allelic.richness(dat, min.n = g)$Ar, error = function(e) e)
  if (inherits(ar_hf, "error") || !identical(dim(as.matrix(ar_hf)), dim(rs$Ar))) {
    warning("Skipped the hierfstat::allelic.richness() cross-check: ",
            if (inherits(ar_hf, "error")) paste0("it failed with \"", conditionMessage(ar_hf), "\".")
            else paste0("its result has a different shape than this run's (", length(rows),
                        " records). This usually means too few records survived to trust ",
                        "these numbers."), call. = FALSE)
  } else {
    d_ar <- max(abs(as.matrix(ar_hf) - rs$Ar), na.rm = TRUE)
    .inform(verbose, sprintf("  cross-check vs hierfstat::allelic.richness(): max diff %.2e", d_ar))
    if (d_ar > 1e-8) warning("Allelic richness disagrees with hierfstat.", call. = FALSE)
  }
  invisible(NULL)
}
