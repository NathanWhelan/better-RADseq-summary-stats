###############################################################################
#
#  R/hwe_test.R -- per-locus test of Hardy-Weinberg equilibrium (HWE). A
#  REPORT, never a filter: no function in this package removes loci based
#  on this one's output (see @details in hwe_test() below for why -- the
#  short version is that the common practice of pooling populations, testing
#  HWE, and dropping loci that fail measurably UNDERESTIMATES real population
#  structure, by removing exactly the loci most informative about it).
#
#  THE TEST ITSELF: Levene (1949)/Haldane (1954)'s exact conditional
#  probability test -- "given how many copies of each allele this locus has,
#  how likely is the exact genotype arrangement observed, compared to every
#  OTHER arrangement with the same allele counts?" `method = "exact"`
#  (the default) always uses this; `method = "chisq"` uses the older,
#  faster, approximate chi-square goodness-of-fit test instead, kept
#  available as a quick first pass over a large dataset.
#
#  TWO WAYS TO GET THE EXACT P-VALUE, chosen automatically per locus
#  depending on how many alleles were actually observed there (this matters
#  because a RAD-haplotype locus can have many alleles, while a SNP locus
#  never has more than 2):
#
#    BIALLELIC   (submethod "exact-enum") -- full enumeration. With only 2
#                alleles, the entire genotype table is pinned down by ONE
#                number (the heterozygote count), which can only take a
#                handful of feasible values -- cheap to list out completely,
#                so the exact p-value is computed with NO Monte Carlo error
#                at all.
#    MULTIALLELIC (submethod "exact-mc")  -- Guo & Thompson's (1992) DIRECT
#                Monte Carlo method. With 3+ alleles, listing every feasible
#                genotype table is not tractable in general (unlike the
#                biallelic case, the number of feasible tables grows with
#                SAMPLE SIZE as well as allele count), so instead: draw a
#                uniformly random re-pairing of this locus's own allele
#                copies, over and over, and see how often the resulting
#                table is at least as improbable as the one actually
#                observed. Guo & Thompson (1992) proposed TWO Monte Carlo
#                methods for this -- this package implements their DIRECT
#                one (repeatedly draw a fresh random permutation and count),
#                not the Markov-chain ("switching") one GENEPOP uses. That
#                choice was made deliberately, not by default: each direct
#                draw is an EXACT, independent sample from the true null
#                distribution (no burn-in, no "has the chain mixed yet?"
#                question, unlike a Markov chain, whose own original paper
#                describes its mixing time as "guessed at but never fully
#                determined" -- Huber et al. 2006), and its per-draw cost
#                (linear in sample size) is irrelevant at the sample sizes
#                RADseq studies actually have (tens to low hundreds of
#                individuals) -- GENEPOP's own switching method exists to
#                help at the much larger sample sizes typical of human
#                genetics, not the setting this package targets.
#
#  VALIDATION. The formula both branches share (.levene_log_weight() below)
#  was checked two ways before being used anywhere in this file, directly in
#  R, not assumed from the papers' notation: (1) for a small worked
#  biallelic example, the normalized probabilities of every feasible table
#  summed to exactly 1; (2) for a small worked 3-allele example, the
#  formula's probabilities matched brute-force enumeration of every gene-
#  copy arrangement to machine precision (max difference 1.1e-16, and they
#  also summed to 1). See tests/testthat/test-hwe-test.R for the equivalent
#  checks kept as permanent regression tests, plus a cross-check against
#  pegas::hw.test() (which implements Guo & Thompson's OTHER, Markov-chain
#  method) where installed.
#
#  REFERENCES
#    Levene, H. (1949) On a matching problem arising in genetics. Annals of
#      Mathematical Statistics 20: 91-94.
#    Haldane, J.B.S. (1954) An exact test for randomness of mating. Journal
#      of Genetics 52: 631-635.
#    Guo, S.-W. & Thompson, E.A. (1992) Performing the exact test of
#      Hardy-Weinberg proportion for multiple alleles. Biometrics 48:
#      361-372. https://doi.org/10.2307/2532296 (see the paper itself, or
#      Huber et al. 2006 below, for both of the two Monte Carlo methods it
#      proposes -- this file implements the DIRECT one, not the chain).
#    Huber, M., Chen, Y., Dinwoodie, I., Dobra, A. & Nicholas, M. (2006)
#      Monte Carlo algorithms for Hardy-Weinberg proportions. Biometrics
#      62: 49-53. https://doi.org/10.1111/j.1541-0420.2005.00418.x (the
#      source used here for how Guo & Thompson's direct method works, and
#      for Levene's formula, eqs. 1-2 there).
#    Besag, J. & Clifford, P. (1991) Sequential Monte Carlo p-values.
#      Biometrika 78: 301-304. (Early stopping in the multiallelic branch.)
#    Phipson, B. & Smyth, G.K. (2010) Permutation p-values should never be
#      zero. Statistical Applications in Genetics and Molecular Biology 9:
#      Article 39. (Why a Monte Carlo p-value is (ge + 1) / (n + 1).)
#    Raymond, M. & Rousset, F. (1995) GENEPOP (version 1.2): population
#      genetics software for exact tests and ecumenicism. Journal of
#      Heredity 86: 248-249. (What GENEPOP does differently -- a modified
#      Markov-chain/switching version of Guo & Thompson's method -- and
#      why this package chose the direct method instead; see above.)
#    Pearman, W.S., Urban, L. & Alexander, A. (2022) Commonly used
#      Hardy-Weinberg equilibrium filtering schemes impact population
#      structure inferences using RADseq data. Molecular Ecology Resources
#      22: 2599-2613. (Why this function only ever reports -- see
#      hwe_test()'s @details.)
#
###############################################################################

## ---------------------------------------------------------------------------
## Shared building blocks
## ---------------------------------------------------------------------------

## Not exported. Builds the (k x k) genotype-count table for one locus, one
## group of samples: cell [g, h] for g <= h holds how many TYPED individuals
## have the unordered genotype {allele g, allele h} -- g == h is a
## homozygote count, g < h a heterozygote count; cells with g > h are always
## 0 and unused. Untyped (missing) individuals are simply excluded first.
.hwe_geno_table <- function(a, b, k) {
  ok <- !is.na(a) & !is.na(b)
  a <- a[ok]; b <- b[ok]
  gi <- pmin(a, b); gj <- pmax(a, b)   # gi <= gj always, by construction
  idx <- (gj - 1L) * k + gi            # linear index into a k x k matrix
  tab <- tabulate(idx, nbins = k * k)
  dim(tab) <- c(k, k)
  tab
}

## Not exported. Levene's (1949) exact conditional probability of a genotype
## table -- or rather, the part of it that actually varies from one table to
## another when every table being compared shares the same typed-individual
## count and the same allele-count vector (always true below: sampling only
## ever reshuffles which allele pairs with which, never how many copies of
## each allele exist). In full:
##
##   P(table) = N! * 2^h * prod_i(n_i!)  /  [ (2N)! * prod_{i<=j}(n_ij!) ]
##
## N = typed individuals, n_i = allele i's copy count, n_ij = this table's
## own genotype counts, h = total heterozygote count = sum of the off-
## diagonal (i<j) counts. N!, (2N)! and prod_i(n_i!) are IDENTICAL for every
## table considered here (fixed by the data, not by which table this is),
## so they are left out below -- what remains is proportional to the true
## probability, in a way that a RATIO of sums over several tables (all a
## p-value needs) is completely unaffected by.
.levene_log_weight <- function(tab) {
  hom <- diag(tab)
  het <- tab[upper.tri(tab)]
  -sum(lfactorial(hom)) - sum(lfactorial(het)) + log(2) * sum(het)
}

## Not exported. Tolerance for the "at least as extreme" comparison of log
## weights. Two tables with the same true probability can get log weights
## that differ in the last few bits when the sums are taken in a different
## order (e.g. the vectorised Monte Carlo branch below vs
## .levene_log_weight()), and without a tolerance such a tie would be
## counted as LESS extreme. 1e-7 on the log scale is a relative probability
## difference of 1e-7 -- far below anything that could change a p-value.
.hwe_tol <- 1e-7

## Not exported. Given a1/a2 (this group's two allele-index vectors at one
## locus, some entries possibly NA for a missing genotype) and `k` (the
## locus's declared allele count), returns a compact, RENUMBERED version
## using only the alleles actually seen among typed individuals here --
## e.g. if only alleles 1 and 3 of a declared 3-allele locus are ever
## observed in this particular group, they are relabelled 1 and 2 so the
## biallelic-vs-multiallelic decision below (and the exact-enum branch,
## which assumes alleles 1 and 2) is based on what is REALLY present in
## this group, not on how many alleles the VCF happened to declare -- the
## same principle write_plink()/kinship.R already apply via
## locus_allele_stats()$n_observed_alleles.
.hwe_compact <- function(a, b, k) {
  ok <- !is.na(a) & !is.na(b)
  counts <- tabulate(c(a[ok], b[ok]), nbins = k)
  observed <- which(counts > 0L)
  list(a = match(a, observed), b = match(b, observed), k_obs = length(observed))
}

## ---------------------------------------------------------------------------
## Exact test: biallelic branch (full enumeration, no Monte Carlo error)
## ---------------------------------------------------------------------------

## Not exported. `tab` is a 2x2 genotype table (tab[1,1]/tab[2,2] the two
## homozygote counts, tab[1,2] the heterozygote count). Every FEASIBLE
## heterozygote count is enumerated directly -- there are at most
## min(n1, n2) + 1 of them (n1/n2 the two allele copy counts), never more
## than n_called + 1, so this is always cheap, unlike full enumeration of a
## multi-allele table (see the multiallelic branch below for why that one
## needs Monte Carlo instead). The p-value is the total probability of
## every feasible table AT LEAST AS EXTREME as (no more probable than) the
## one actually observed -- the standard "exact test" definition (Levene
## 1949; Haldane 1954).
.hwe_exact_biallelic <- function(tab) {
  n1 <- 2L * tab[1, 1] + tab[1, 2]   # allele-1 copy count
  n2 <- 2L * tab[2, 2] + tab[1, 2]   # allele-2 copy count
  feasible <- seq.int(n1 %% 2L, min(n1, n2), by = 2L)
  n11 <- (n1 - feasible) / 2; n22 <- (n2 - feasible) / 2
  w <- -lfactorial(n11) - lfactorial(feasible) - lfactorial(n22) + log(2) * feasible
  obs_w <- w[feasible == tab[1, 2]]
  ## Shift by the largest weight before exponentiating -- w can be a very
  ## large negative number at realistic sample sizes (lfactorial(200) is
  ## already in the thousands), and exp() of that alone would silently
  ## underflow to 0 for every table at once, turning a well-defined ratio
  ## into 0/0. Subtracting the max first leaves the RATIO (all this
  ## function reports) exactly unchanged while keeping every exponentiated
  ## term safely inside floating-point range.
  mx <- max(w)
  p <- sum(exp(w[w <= obs_w + .hwe_tol] - mx)) / sum(exp(w - mx))
  list(statistic = NA_real_, df = NA_real_, p_value = p, p_value_se = NA_real_,
       n_draws_used = NA_integer_, pct_low_expected = NA_real_,
       submethod = "exact-enum")
}

## ---------------------------------------------------------------------------
## Exact test: multiallelic branch (Guo & Thompson 1992's direct Monte Carlo)
## ---------------------------------------------------------------------------

## Not exported. `a`/`b` are this group's two (already NA-free, already
## compactly renumbered) allele-index vectors at one locus with k_obs >= 3
## observed alleles. Draws independent uniformly-random re-pairings of this
## locus's own 2*n allele copies (n = typed individuals) -- each one an
## EXACT sample from the true null distribution by construction, not an
## approximation that needs time to "mix" -- and counts the draws at least
## as extreme as (no more probable than) the table actually observed.
##
## SEQUENTIAL STOPPING (Besag & Clifford 1991). Drawing stops as soon as
## `stop_after` extreme draws have been seen, or after `n_draws` draws,
## whichever comes first:
##   * stopped early, after `used` draws     -> p = stop_after / used
##   * reached n_draws with fewer exceedances -> p = (ge + 1) / (n_draws + 1)
## Both are valid p-values. The second form also means p is never exactly
## 0: a finite simulation cannot show that (Phipson & Smyth 2010). A locus
## nowhere near significance stops after a few dozen draws instead of all
## n_draws, which is where almost all the time used to go.
## `stop_after = Inf` switches stopping off (always n_draws draws). The
## standard error is the binomial one over the draws actually used.
##
## VECTORISED: draws are made a batch at a time. Ordering a row of uniform
## random numbers gives a uniformly random permutation of that row, and one
## order() call handles every row of a batch at once; each draw's genotype
## table then comes from a single tabulate() call, and its Levene log weight
## is computed row-wise. Same test and same null distribution as one
## sample() per draw, without an R-level loop per draw. The batch starts
## small (most loci stop early) and doubles each round.
.hwe_exact_multiallelic <- function(a, b, k_obs, n_draws, stop_after = 20) {
  obs_w <- .levene_log_weight(.hwe_geno_table(a, b, k_obs))
  n <- length(a); n2 <- 2L * n
  copies <- c(a, b)   # the locus's own 2n allele copies, as a flat multiset
  kk <- k_obs * k_obs
  hom_cells <- which(diag(k_obs) == 1)           # (i, i) cells of a k x k table
  het_cells <- which(upper.tri(diag(k_obs)))     # (i, j) cells, i < j
  odd <- seq.int(1L, n2, 2L); even <- odd + 1L
  ge <- 0L; used <- 0L; batch <- 64L
  while (used < n_draws && ge < stop_after) {
    m <- as.integer(min(batch, n_draws - used))
    ## Row r of `pos` is a random permutation of 1..2n: sorting on
    ## (row number + uniform) keeps each row's block together and shuffles
    ## within it.
    o <- order(rep(seq_len(m), each = n2) + stats::runif(m * n2))
    pos <- matrix(o, m, n2, byrow = TRUE) - (seq_len(m) - 1L) * n2
    perm <- matrix(copies[pos], m, n2)
    g1 <- perm[, odd, drop = FALSE]; g2 <- perm[, even, drop = FALSE]
    cell <- (pmax(g1, g2) - 1L) * k_obs + pmin(g1, g2)   # as in .hwe_geno_table()
    tab <- matrix(tabulate(as.vector(cell) + rep.int((seq_len(m) - 1L) * kk, n),
                           nbins = m * kk), m, kk, byrow = TRUE)
    het <- tab[, het_cells, drop = FALSE]
    w <- -rowSums(lfactorial(tab[, hom_cells, drop = FALSE])) -
         rowSums(lfactorial(het)) + log(2) * rowSums(het)
    hits <- cumsum(w <= obs_w + .hwe_tol)
    if (ge + hits[m] >= stop_after) {
      used <- used + which(ge + hits >= stop_after)[1L]   # stop exactly at the hit
      ge <- as.integer(stop_after)
    } else {
      used <- used + m; ge <- ge + hits[m]
    }
    batch <- min(2L * batch, 4096L)
  }
  p <- if (ge >= stop_after) ge / used else (ge + 1) / (used + 1)
  se <- sqrt(p * (1 - p) / used)
  list(statistic = NA_real_, df = NA_real_, p_value = p, p_value_se = se,
       n_draws_used = as.integer(used), pct_low_expected = NA_real_,
       submethod = "exact-mc")
}

## ---------------------------------------------------------------------------
## Chi-square goodness-of-fit test (method = "chisq")
## ---------------------------------------------------------------------------

## Not exported. `tab` is a k_obs x k_obs genotype table (compactly
## renumbered, i.e. every allele 1..k_obs is actually observed). Standard
## asymptotic HWE chi-square test: expected genotype counts from the
## observed allele frequencies, df = k(k-1)/2 (the number of genotype
## classes, k(k+1)/2, minus the k-1 allele-frequency parameters estimated
## from the data, minus 1 more for the counts summing to n -- a standard,
## textbook result). `pct_low_expected` flags the classic caveat: the
## chi-square approximation itself is unreliable when expected counts are
## small (conventionally, < 5), so a user can see when that is happening
## rather than trusting a p-value blindly -- the same role this column
## plays nowhere else in this package except here, since it is specific to
## the asymptotic (not exact) test.
.hwe_chisq <- function(tab, k_obs) {
  n <- sum(tab)
  ## Allele i's total copy count: its own diagonal (homozygote) cell counts
  ## TWICE (2 copies per homozygous individual), every off-diagonal cell in
  ## its row or column ONCE (1 copy per heterozygous individual) -- which is
  ## exactly what rowSums(tab)[i] + colSums(tab)[i] adds up to, since `tab`
  ## is upper-triangular (tab[i,i] is counted once in each sum, i.e. twice
  ## total, and every off-diagonal cell touching allele i is counted once).
  allele_n <- rowSums(tab) + colSums(tab)
  p <- allele_n / (2 * n)
  exp_tab <- matrix(0, k_obs, k_obs)
  for (i in seq_len(k_obs)) for (j in i:k_obs)
    exp_tab[i, j] <- if (i == j) n * p[i]^2 else 2 * n * p[i] * p[j]
  obs_cells <- c(diag(tab), tab[upper.tri(tab)])
  exp_cells <- c(diag(exp_tab), exp_tab[upper.tri(exp_tab)])
  stat <- sum((obs_cells - exp_cells)^2 / exp_cells)
  df <- k_obs * (k_obs - 1) / 2
  list(statistic = stat, df = df, p_value = stats::pchisq(stat, df, lower.tail = FALSE),
       p_value_se = NA_real_, n_draws_used = NA_integer_,
       pct_low_expected = 100 * mean(exp_cells < 5), submethod = "chisq")
}

## ---------------------------------------------------------------------------
## The exported function
## ---------------------------------------------------------------------------

#' Per-locus test of Hardy-Weinberg equilibrium
#'
#' Tests every locus in `H` for departure from Hardy-Weinberg equilibrium
#' (HWE), pooled across all samples (`pops = NULL`, the default) or
#' separately per population (`pops` given, same shape as
#' [filter_call_rate()]'s `pops` argument). This is a **report, never a
#' filter**: no function in this package removes loci based on this one's
#' output -- see `@details` for why.
#'
#' `method = "exact"` (the default) is Levene (1949)/Haldane (1954)'s exact
#' conditional-probability test, computed exactly (no Monte Carlo error) for
#' biallelic loci via full enumeration, and via Guo & Thompson's (1992)
#' direct Monte Carlo method for loci with 3+ observed alleles (RAD
#' haplotypes). `method = "chisq"` is the older, approximate chi-square
#' goodness-of-fit test, kept available as a fast first pass. See
#' `R/hwe_test.R`'s file header for the full algorithm description, and
#' `tests/testthat/test-hwe-test.R` for how both are validated.
#'
#' @details
#' **Why this only ever reports, never filters.** Pearman, Urban & Alexander
#' (2022) show that the most common HWE-filtering scheme -- pool every
#' population together, test, and drop any locus that fails -- measurably
#' **underestimates real population structure**, because it tends to remove
#' exactly the loci most informative about that structure (a locus that
#' looks like it deviates from HWE when pooled is often one where allele
#' frequencies genuinely differ between populations, the Wahlund effect,
#' not a genotyping artifact). Testing separately per population (`pops`
#' given) avoids conflating the two, which is why that option exists here --
#' but even then, this function stops at reporting.
#'
#' **Monte Carlo p-values for multiallelic loci.** The p-value for a locus
#' with 3+ observed alleles is estimated by simulation and carries its own
#' standard error (`p_value_se`, over `n_draws_used` draws). Simulation
#' stops early once `stop_after` draws at least as extreme as the observed
#' table have been seen (Besag & Clifford 1991): a locus nowhere near
#' significance needs only a few dozen draws, while a locus with a small
#' p-value runs to `n_draws`, so precision goes where it matters. A p-value
#' is never reported as exactly 0; the smallest possible value is
#' `1 / (n_draws + 1)` (Phipson & Smyth 2010). Biallelic loci
#' (`submethod = "exact-enum"`) have no Monte Carlo error -- they are exact.
#'
#' **References.** Levene, H. (1949) On a matching problem arising in
#' genetics. *Annals of Mathematical Statistics* 20:91-94. -- Haldane, J.B.S.
#' (1954) An exact test for randomness of mating. *Journal of Genetics*
#' 52:631-635. -- Guo, S.-W. & Thompson, E.A. (1992) Performing the exact
#' test of Hardy-Weinberg proportion for multiple alleles. *Biometrics*
#' 48:361-372. -- Huber, M., Chen, Y., Dinwoodie, I., Dobra, A. & Nicholas,
#' M. (2006) Monte Carlo algorithms for Hardy-Weinberg proportions.
#' *Biometrics* 62:49-53. <https://doi.org/10.1111/j.1541-0420.2005.00418.x>
#' -- Pearman, W.S., Urban, L. & Alexander, A. (2022) Commonly used
#' Hardy-Weinberg equilibrium filtering schemes impact population structure
#' inferences using RADseq data. *Molecular Ecology Resources* 22:2599-2613.
#' -- Besag, J. & Clifford, P. (1991) Sequential Monte Carlo p-values.
#' *Biometrika* 78:301-304. -- Phipson, B. & Smyth, G.K. (2010) Permutation
#' p-values should never be zero. *Statistical Applications in Genetics and
#' Molecular Biology* 9:39.
#'
#' @param H A list as returned by [read_stacks_vcf()] (or by another filter in
#'   this package, since they all return the same shape).
#' @param pops Optional named list of sample-ID vectors (from
#'   [read_popmap()]). `NULL` (the default) pools every sample into one
#'   group named `"pooled"`; given, tests each population separately.
#' @param method `"exact"` (default) or `"chisq"`.
#' @param n_draws Maximum Monte Carlo draws per multiallelic locus under
#'   `method = "exact"`. Default `10000L`. Ignored for biallelic loci
#'   (exact by enumeration) and under `method = "chisq"`.
#' @param stop_after Stop drawing for a locus once this many draws at least
#'   as extreme as the observed table have been seen (Besag & Clifford
#'   1991). Default `20`. `Inf` always uses all `n_draws` draws.
#' @param seed Random seed for the Monte Carlo draws, for reproducibility.
#'   Default `2024`. The caller's own RNG state is restored when this
#'   function returns (see [diversity_stats()] for the same convention).
#' @param verbose Print a short summary. Default `TRUE`.
#' @return A data frame, one row per (locus, population) combination:
#'   `locus`, `population`, `method`, `submethod` (`"exact-enum"`,
#'   `"exact-mc"`, or `"chisq"`), `n_alleles` (declared), `n_observed_alleles`,
#'   `n_called` (typed individuals in that group), `statistic`, `df`,
#'   `p_value`, `p_value_se` and `n_draws_used` (only defined for
#'   `"exact-mc"` rows),
#'   `pct_low_expected` (only defined for `"chisq"` rows). A locus/group
#'   combination with fewer than 2 observed alleles (monomorphic there) or
#'   fewer than 2 typed individuals gets `NA` throughout except `n_called`/
#'   `n_observed_alleles` themselves.
#' @examples
#' # 2 loci: locus_1 biallelic (12 individuals, exact-enum path), locus_2
#' # triallelic (exact-mc path).
#' samp <- paste0("s", 1:12)
#' H <- list(
#'   A1 = rbind(locus_1 = rep(c(1L,1L,2L), 4), locus_2 = rep(c(1L,2L,3L), 4)),
#'   A2 = rbind(locus_1 = rep(c(1L,2L,2L), 4), locus_2 = rep(c(2L,3L,1L), 4)),
#'   locus = c("locus_1","locus_2"), locus_raw = c("locus_1","locus_2"),
#'   n_alleles = c(2L, 3L), alleles = list(c("A","C"), c("A","C","T")),
#'   samples = samp
#' )
#' colnames(H$A1) <- colnames(H$A2) <- samp
#' hwe_test(H, n_draws = 500, verbose = FALSE)
#' @export
hwe_test <- function(H, pops = NULL, method = "exact", n_draws = 10000L,
                      stop_after = 20, seed = 2024, verbose = TRUE) {
  if (!(length(method) == 1L && method %in% c("exact", "chisq")))
    stop("method must be exactly \"exact\" or \"chisq\" (got: ", paste(method, collapse = ", "), ").")
  n_draws <- as.integer(n_draws)
  if (is.na(n_draws) || n_draws < 1L)
    stop("n_draws must be a positive integer.")
  if (length(stop_after) != 1L || is.na(stop_after) || stop_after < 1)
    stop("stop_after must be a single number >= 1 (Inf for a fixed n_draws).")

  ## Same RNG-preservation convention as every other stochastic function in
  ## this package: restore the caller's own random-number state on exit.
  restore_rng <- .save_rng_state()
  on.exit(restore_rng(), add = TRUE)
  set.seed(seed)

  groups <- if (is.null(pops)) stats::setNames(list(H$samples), "pooled") else pops
  n_rec <- nrow(H$A1)
  n_out <- n_rec * length(groups)

  ## One slot per (population, locus) row, filled in place and turned into a
  ## data frame once at the end -- building a one-row data.frame() per locus
  ## and rbind()-ing tens of thousands of them was a large share of the run
  ## time on a real dataset.
  submethod <- rep(NA_character_, n_out)
  n_obs_all <- n_called <- integer(n_out)
  statistic <- df <- p_value <- p_value_se <- pct_low <- rep(NA_real_, n_out)
  n_used <- rep(NA_integer_, n_out)
  ri <- 0L

  for (g in names(groups)) {
    ids <- groups[[g]]
    for (j in seq_len(n_rec)) {
      ri <- ri + 1L
      cp <- .hwe_compact(H$A1[j, ids], H$A2[j, ids], H$n_alleles[j])
      ok <- !is.na(cp$a) & !is.na(cp$b)
      n_obs_all[ri] <- cp$k_obs
      n_called[ri]  <- sum(ok)
      if (cp$k_obs < 2L || n_called[ri] < 2L) next  # monomorphic here, or < 2 typed
      a <- cp$a[ok]; b <- cp$b[ok]
      res <- if (method == "chisq") {
        .hwe_chisq(.hwe_geno_table(a, b, cp$k_obs), cp$k_obs)
      } else if (cp$k_obs == 2L) {
        .hwe_exact_biallelic(.hwe_geno_table(a, b, cp$k_obs))
      } else {
        .hwe_exact_multiallelic(a, b, cp$k_obs, n_draws, stop_after = stop_after)
      }
      submethod[ri] <- res$submethod
      statistic[ri] <- res$statistic; df[ri] <- res$df
      p_value[ri] <- res$p_value; p_value_se[ri] <- res$p_value_se
      n_used[ri] <- res$n_draws_used; pct_low[ri] <- res$pct_low_expected
    }
  }
  out <- data.frame(
    locus = rep(H$locus, length(groups)),
    population = rep(names(groups), each = n_rec),
    method = rep(method, n_out), submethod = submethod,
    n_alleles = rep(as.integer(H$n_alleles), length(groups)),
    n_observed_alleles = n_obs_all, n_called = n_called,
    statistic = statistic, df = df, p_value = p_value, p_value_se = p_value_se,
    n_draws_used = n_used, pct_low_expected = pct_low,
    stringsAsFactors = FALSE)

  if (verbose) {
    n_bad <- sum(is.na(out$p_value))
    message(sprintf(
      "hwe_test() [method = \"%s\"]: %s (locus, population) combinations, %s skipped (monomorphic or < 2 typed individuals there)",
      method, format(nrow(out), big.mark = ","), format(n_bad, big.mark = ",")))
    if (method == "exact") {
      n_enum <- sum(out$submethod == "exact-enum", na.rm = TRUE)
      n_mc   <- sum(out$submethod == "exact-mc", na.rm = TRUE)
      message(sprintf("  %s biallelic (exact-enum, no Monte Carlo error), %s multiallelic (exact-mc, up to %s draws, stopping after %s extreme draws%s)",
                      format(n_enum, big.mark = ","), format(n_mc, big.mark = ","),
                      format(n_draws, big.mark = ","), format(stop_after),
                      if (n_mc) sprintf("; median %s draws used",
                                        format(stats::median(out$n_draws_used, na.rm = TRUE), big.mark = ","))
                      else ""))
    }
  }
  out
}
