## Same make_H() helper convention as the other test files in this package.
make_H <- function(A1, A2, locus = NULL, locus_raw = NULL, n_alleles = NULL,
                    alleles = NULL, samples = NULL) {
  n_rec <- nrow(A1)
  if (is.null(samples)) samples <- colnames(A1)
  if (is.null(samples)) samples <- paste0("s", seq_len(ncol(A1)))
  if (is.null(locus)) locus <- rownames(A1)
  dimnames(A1) <- dimnames(A2) <- list(NULL, samples)
  if (is.null(locus)) locus <- paste0("locus_", seq_len(n_rec))
  if (is.null(locus_raw)) locus_raw <- locus
  if (is.null(n_alleles)) n_alleles <- rep(2L, n_rec)
  if (is.null(alleles)) alleles <- replicate(n_rec, c("A", "C"), simplify = FALSE)
  list(A1 = A1, A2 = A2, locus = locus, locus_raw = locus_raw,
       alleles = alleles, n_alleles = n_alleles, samples = samples)
}

## ---------------------------------------------------------------------------
## .levene_log_weight(): the formula both exact branches share
## ---------------------------------------------------------------------------

test_that(".levene_log_weight() normalized probabilities sum to 1 (biallelic)", {
  ## 4 individuals, allele counts n1=5, n2=3 -- feasible heterozygote counts
  ## are {1, 3} (same parity as n1, bounded by min(n1,n2)=3).
  feasible <- c(1, 3)
  w <- vapply(feasible, function(n12) {
    n11 <- (5 - n12) / 2; n22 <- (3 - n12) / 2
    tab <- matrix(c(n11, 0, n12, n22), 2, 2)
    RADdiversity:::.levene_log_weight(tab)
  }, numeric(1))
  p <- exp(w - max(w)); p <- p / sum(p)
  expect_equal(sum(p), 1, tolerance = 1e-10)
  ## Cross-checked (in this project's own development session, not in this
  ## test) against brute-force enumeration of all C(8,3)=56 gene-copy
  ## arrangements: p(n12=1) = 0.4285714, p(n12=3) = 0.5714286.
  expect_equal(sort(p), c(0.4285714, 0.5714286), tolerance = 1e-6)
})

test_that(".levene_log_weight() matches brute-force enumeration (3 alleles)", {
  ## N=5 individuals (10 gene copies), allele counts n1=4, n2=4, n3=2.
  ## Ground truth: every C(10,4)*C(6,4) = 210*15 = 3150 way of assigning
  ## the 3 allele labels to the 10 gene-copy "slots" is equally likely;
  ## brute-force count how often each genotype table occurs, and compare
  ## to the Levene-formula-implied probability.
  slots <- 1:10
  combosA <- combn(10, 4, simplify = FALSE)
  geno_of <- function(a_slots) {
    rest <- setdiff(slots, a_slots)
    combosB <- combn(rest, 4, simplify = FALSE)
    t(vapply(combosB, function(b_slots) {
      c_slots <- setdiff(rest, b_slots)
      lab <- integer(10); lab[a_slots] <- 1L; lab[b_slots] <- 2L; lab[c_slots] <- 3L
      g1 <- lab[seq(1, 10, 2)]; g2 <- lab[seq(2, 10, 2)]
      tab <- RADdiversity:::.hwe_geno_table(g1, g2, 3)
      c(tab[1,1], tab[2,2], tab[3,3], tab[1,2], tab[1,3], tab[2,3])
    }, numeric(6)))
  }
  all_tabs <- do.call(rbind, lapply(combosA, geno_of))
  key <- apply(all_tabs, 1, paste, collapse = ",")
  emp <- table(key) / nrow(all_tabs)

  uniq <- unique(all_tabs)
  w <- apply(uniq, 1, function(row) {
    tab <- matrix(0, 3, 3)
    tab[1,1] <- row[1]; tab[2,2] <- row[2]; tab[3,3] <- row[3]
    tab[1,2] <- row[4]; tab[1,3] <- row[5]; tab[2,3] <- row[6]
    RADdiversity:::.levene_log_weight(tab)
  })
  p_formula <- exp(w) / sum(exp(w))
  names(p_formula) <- apply(uniq, 1, paste, collapse = ",")

  cmp <- p_formula[names(emp)]
  expect_equal(max(abs(as.numeric(emp) - as.numeric(cmp))), 0, tolerance = 1e-10)
  expect_equal(sum(p_formula), 1, tolerance = 1e-10)
})

## ---------------------------------------------------------------------------
## Biallelic exact-enum branch
## ---------------------------------------------------------------------------

test_that("the biallelic exact-enum branch matches a hand-computed example", {
  ## Same 4-individual, n1=5/n2=3 example as above. Observed table: n12=3
  ## (n11=1, n22=0) -- both feasible tables (n12=1, n12=3) are AT LEAST AS
  ## extreme as n12=3 itself (it's the higher-probability one), so the
  ## exact p-value should be exactly 1.
  tab <- matrix(c(1, 0, 3, 0), 2, 2)
  res <- RADdiversity:::.hwe_exact_biallelic(tab)
  expect_equal(res$p_value, 1, tolerance = 1e-10)
  expect_identical(res$submethod, "exact-enum")
  expect_true(is.na(res$p_value_se))

  ## Observed table: n12=1 (n11=2, n22=1) -- the LESS probable of the two
  ## feasible tables, so its own p-value should equal its own probability,
  ## 0.4285714 (from the sum-to-1 test above).
  tab2 <- matrix(c(2, 0, 1, 1), 2, 2)
  res2 <- RADdiversity:::.hwe_exact_biallelic(tab2)
  expect_equal(res2$p_value, 0.4285714, tolerance = 1e-6)
})

test_that("hwe_test() routes biallelic loci to exact-enum and reports p_value_se = NA there", {
  H <- make_H(A1 = rbind(c(1,1,1,1,2,2,2,1)), A2 = rbind(c(1,1,2,1,2,2,1,2)))
  res <- suppressMessages(hwe_test(H, verbose = FALSE))
  expect_equal(res$submethod, "exact-enum")
  expect_true(is.na(res$p_value_se))
  expect_equal(res$n_observed_alleles, 2L)
})

## ---------------------------------------------------------------------------
## Multiallelic exact-mc branch
## ---------------------------------------------------------------------------

test_that("hwe_test() routes multiallelic loci to exact-mc with a binomial SE", {
  set.seed(1)
  n <- 15
  a1 <- sample(1:3, n, replace = TRUE); a2 <- sample(1:3, n, replace = TRUE)
  H <- make_H(A1 = rbind(a1), A2 = rbind(a2), n_alleles = 3L,
              alleles = list(c("A","C","T")))
  res <- suppressMessages(hwe_test(H, n_draws = 2000, verbose = FALSE))
  expect_equal(res$submethod, "exact-mc")
  expect_true(res$p_value_se > 0 && res$p_value_se < 0.02)
  expect_equal(res$p_value_se, sqrt(res$p_value * (1 - res$p_value) / 2000), tolerance = 1e-10)
})

test_that("the exact-mc branch agrees with pegas::hw.test() (Guo & Thompson's own direct method)", {
  testthat::skip_if_not_installed("pegas")
  set.seed(42)
  n <- 20
  a1 <- sample(1:3, n, replace = TRUE, prob = c(.5, .3, .2))
  a2 <- sample(1:3, n, replace = TRUE, prob = c(.5, .3, .2))
  H <- make_H(A1 = rbind(a1), A2 = rbind(a2), n_alleles = 3L,
              alleles = list(c("A","C","T")))

  ## pegas conversion, written only here (test-only, not a runtime
  ## dependency of the package) -- pegas::alleles2loci() wants a matrix
  ## with one column per allele copy.
  m <- cbind(as.character(a1), as.character(a2))
  colnames(m) <- c("L_a1", "L_a2")
  lo <- pegas::alleles2loci(m)
  peg <- pegas::hw.test(lo, B = 20000)

  mine <- suppressMessages(hwe_test(H, n_draws = 20000, seed = 1, verbose = FALSE))

  ## Both implementations estimate the SAME exact quantity via independent
  ## direct Monte Carlo draws (pegas's own hw.test.loci() source, read
  ## this session, uses the identical "sample() the flat allele-copy
  ## multiset, reshape into pairs, rank by -sum(lfactorial(counts)) +
  ## log(2)*sum(heterozygote counts)" procedure this package implements
  ## independently) -- so they should agree within a few combined standard
  ## errors, not to Monte Carlo-defying precision.
  combined_se <- sqrt(mine$p_value_se^2 + (0.5 / sqrt(20000))^2)
  expect_lt(abs(mine$p_value - peg[1, "Pr.exact"]), 6 * combined_se)
})

test_that("small-table agreement: exact-mc lands within its own SE of brute-force truth", {
  ## N=5 individuals, 3 alleles, allele counts n1=4/n2=4/n3=2 (same table
  ## family validated above) -- compute the TRUE p-value for one specific
  ## observed table by brute-force enumeration (independent of both
  ## .levene_log_weight() and .hwe_exact_multiallelic()), then check the
  ## Monte Carlo estimate lands within a generous multiple of its own
  ## reported SE.
  slots <- 1:10
  combosA <- combn(10, 4, simplify = FALSE)
  all_tabs <- do.call(rbind, lapply(combosA, function(a_slots) {
    rest <- setdiff(slots, a_slots)
    combosB <- combn(rest, 4, simplify = FALSE)
    t(vapply(combosB, function(b_slots) {
      c_slots <- setdiff(rest, b_slots)
      lab <- integer(10); lab[a_slots] <- 1L; lab[b_slots] <- 2L; lab[c_slots] <- 3L
      g1 <- lab[seq(1,10,2)]; g2 <- lab[seq(2,10,2)]
      tab <- RADdiversity:::.hwe_geno_table(g1, g2, 3)
      c(tab)
    }, numeric(9)))
  }))
  ## IMPORTANT: normalize from the EMPIRICAL arrangement counts (how many
  ## of the 3150 raw arrangements produced each distinct table), not by
  ## re-applying exp(weight) row-by-row over the (heavily duplicated) raw
  ## arrangement list -- that would double-count each table's own
  ## multiplicity (once via its repeated rows, once again via
  ## exp(weight)), squaring its effective weight. This mirrors exactly how
  ## the "matches brute-force enumeration" test above validates the
  ## formula -- deduplicate first, then compare.
  key <- apply(all_tabs, 1, paste, collapse = ",")
  emp <- table(key) / nrow(all_tabs)
  uniq_tabs <- all_tabs[!duplicated(key), , drop = FALSE]
  uniq_key  <- key[!duplicated(key)]
  w_uniq <- apply(uniq_tabs, 1, function(v) RADdiversity:::.levene_log_weight(matrix(v, 3, 3)))
  names(w_uniq) <- uniq_key

  ## Pick one specific table as "observed": a1=(1,1,1,1,2,2,2,2,3,3) paired
  ## sequentially -> individuals (1,1)(1,1)(2,2)(2,2)(3,3): a genuinely
  ## non-extreme table (2 hom-1, 2 hom-2, 1 hom-3, no heterozygotes).
  a1 <- c(1,1,2,2,3); a2 <- c(1,1,2,2,3)
  obs_w <- RADdiversity:::.levene_log_weight(RADdiversity:::.hwe_geno_table(a1, a2, 3))
  true_p <- sum(emp[names(w_uniq)[w_uniq <= obs_w]])

  set.seed(7)
  mc <- RADdiversity:::.hwe_exact_multiallelic(a1, a2, 3, n_draws = 20000)
  expect_lt(abs(mc$p_value - true_p), 5 * mc$p_value_se + 1e-6)
})

## ---------------------------------------------------------------------------
## method = "chisq"
## ---------------------------------------------------------------------------

test_that("method = \"chisq\" statistic and df match a hand-computed example", {
  ## Biallelic, n=10 individuals: 2 hom-1, 6 het, 2 hom-2 (allele freq 0.5
  ## each) -- expected counts under HWE: E11=2.5, E12=5, E22=2.5 (n*p^2,
  ## 2*n*p*q, n*q^2). chi^2 = (2-2.5)^2/2.5 + (6-5)^2/5 + (2-2.5)^2/2.5
  ##        = 0.1 + 0.2 + 0.1 = 0.4, df = 2*1/2 = 1.
  a1 <- c(rep(1,2), rep(1,6), rep(2,2)); a2 <- c(rep(1,2), rep(2,6), rep(2,2))
  H <- make_H(A1 = rbind(a1), A2 = rbind(a2))
  res <- suppressMessages(hwe_test(H, method = "chisq", verbose = FALSE))
  expect_equal(res$statistic, 0.4, tolerance = 1e-8)
  expect_equal(res$df, 1)
  expect_equal(res$p_value, stats::pchisq(0.4, 1, lower.tail = FALSE), tolerance = 1e-10)
  expect_equal(res$submethod, "chisq")
})

test_that("pct_low_expected flags small expected counts under method = \"chisq\"", {
  ## Tiny n and a rare allele -> some expected genotype counts < 5.
  a1 <- c(1,1,1,1,1,1,1,1,1,2); a2 <- c(1,1,1,1,1,1,1,1,1,1)
  H <- make_H(A1 = rbind(a1), A2 = rbind(a2))
  res <- suppressMessages(hwe_test(H, method = "chisq", verbose = FALSE))
  expect_gt(res$pct_low_expected, 0)
})

## ---------------------------------------------------------------------------
## Uniformity under the null (both exact branches) -- the check most likely
## to catch a wrong weight formula or a wrong "at least as extreme" ranking,
## per the advisor consult that shaped this validation plan.
## ---------------------------------------------------------------------------

test_that("exact p-values do not show excess small values under the null (biallelic)", {
  ## Exact discrete tests are conservative by construction (their null
  ## p-value distribution stochastically dominates Uniform[0,1], a well
  ## known property -- NOT a bug), so this checks the type-I-error-control
  ## direction that actually matters (P(p <= 0.05) is not inflated ABOVE
  ## 0.05), rather than a strict uniformity claim that a correct exact test
  ## would fail.
  set.seed(11)
  n <- 20
  ps <- replicate(300, {
    a1 <- sample(1:2, n, replace = TRUE, prob = c(.6, .4))
    a2 <- sample(1:2, n, replace = TRUE, prob = c(.6, .4))
    tab <- RADdiversity:::.hwe_geno_table(a1, a2, 2)
    RADdiversity:::.hwe_exact_biallelic(tab)$p_value
  })
  expect_lt(mean(ps <= 0.05), 0.10)  # generous margin around the 0.05 target
})

test_that("exact-mc p-values do not show excess small values under the null (multiallelic)", {
  set.seed(12)
  n <- 20
  ps <- replicate(100, {
    a1 <- sample(1:3, n, replace = TRUE, prob = c(.5, .3, .2))
    a2 <- sample(1:3, n, replace = TRUE, prob = c(.5, .3, .2))
    RADdiversity:::.hwe_exact_multiallelic(a1, a2, 3, n_draws = 2000)$p_value
  })
  expect_lt(mean(ps <= 0.05), 0.12)
})

## ---------------------------------------------------------------------------
## API/plumbing: pops, RNG discipline, edge cases
## ---------------------------------------------------------------------------

test_that("pops = NULL pools every sample; pops given tests per population", {
  H <- make_H(A1 = rbind(c(1,1,1,1,2,2,2,1)), A2 = rbind(c(1,1,2,1,2,2,1,2)),
              samples = paste0("s", 1:8))
  pooled <- suppressMessages(hwe_test(H, verbose = FALSE))
  expect_equal(pooled$population, "pooled")
  expect_equal(pooled$n_called, 8L)

  pops <- list(popA = paste0("s", 1:4), popB = paste0("s", 5:8))
  bypop <- suppressMessages(hwe_test(H, pops = pops, verbose = FALSE))
  expect_setequal(bypop$population, c("popA", "popB"))
  expect_equal(nrow(bypop), 2L)
})

test_that("a monomorphic locus/group is reported as NA, not silently computed", {
  H <- make_H(A1 = rbind(c(1,1,1,1)), A2 = rbind(c(1,1,1,1)))
  res <- suppressMessages(hwe_test(H, verbose = FALSE))
  expect_true(is.na(res$p_value))
  expect_equal(res$n_observed_alleles, 1L)
  expect_true(is.na(res$submethod))
})

test_that("fewer than 2 typed individuals in a group is reported as NA", {
  H <- make_H(A1 = rbind(c(1, 2, NA, NA)), A2 = rbind(c(1, 2, NA, NA)))
  res <- suppressMessages(hwe_test(H, verbose = FALSE))
  expect_equal(res$n_called, 2L)  # only 2 typed, fine on its own
  H2 <- make_H(A1 = rbind(c(1, NA, NA, NA)), A2 = rbind(c(1, NA, NA, NA)))
  res2 <- suppressMessages(hwe_test(H2, verbose = FALSE))
  expect_true(is.na(res2$p_value))
  expect_equal(res2$n_called, 1L)
})

test_that("method must be exactly \"exact\" or \"chisq\"", {
  H <- make_H(A1 = rbind(c(1,1,1,1)), A2 = rbind(c(1,2,1,2)))
  expect_error(hwe_test(H, method = "exac", verbose = FALSE), "method must be")
})

test_that("hwe_test() restores the caller's RNG state (uses Monte Carlo internally)", {
  H <- make_H(A1 = rbind(c(1,2,3,1,2,3,1,2,3,1)), A2 = rbind(c(2,3,1,3,1,2,1,2,3,2)),
              n_alleles = 3L, alleles = list(c("A","C","T")))
  set.seed(123)
  before <- .Random.seed
  suppressMessages(hwe_test(H, n_draws = 500, verbose = FALSE))
  after <- .Random.seed
  expect_identical(before, after)
})
