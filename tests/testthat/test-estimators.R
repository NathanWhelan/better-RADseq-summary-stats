## The exported low-level estimators (R/estimators.R). diversity_core_selftest()
## checks them against brute-force Monte Carlo, but that test is skipped on CRAN
## (it takes about 20 seconds), so these closed-form checks run everywhere: each
## one is a value that can be worked out by hand, so a regression in any
## estimator fails here even where the Monte Carlo checks do not run.

test_that("hs_nei_chesser() and hs_biallelic() match the formula, and are vectorised", {
  ## (n/(n-1)) * (1 - sum p^2 - Ho/(2n)) with p = 0.3, Ho = 0.35, n = 10:
  ## sum p^2 = 0.58, so (10/9) * (0.42 - 0.0175).
  expect_equal(hs_biallelic(0.3, ho = 0.35, n = 10), (10 / 9) * (0.42 - 0.35 / 20))
  expect_equal(hs_nei_chesser(0.58, ho = 0.35, n = 10), hs_biallelic(0.3, ho = 0.35, n = 10))
  ## One value per input element, and the same answer element by element.
  many <- hs_biallelic(c(0.3, 0.5), ho = c(0.35, 0.60), n = c(10, 20))
  expect_equal(length(many), 2L)
  expect_equal(many[2], hs_biallelic(0.5, ho = 0.6, n = 20))
  ## Undefined below 2 individuals.
  expect_true(is.na(hs_biallelic(0.3, ho = 0.35, n = 1)))
})

test_that("FIS from Nei & Chesser reaches its bounds, where Stacks' Pi cannot", {
  n <- 10
  ## Every individual heterozygous at a biallelic locus: the true FIS is -1.
  expect_equal(1 - 1 / hs_biallelic(0.5, ho = 1, n = n), -1)
  ## The 2n-corrected estimator behind Stacks' `Pi` stops at -(n-1)/n.
  expect_equal(1 - 1 / hs_stacks_pi(0.5, n = n), -(n - 1) / n)
  ## No heterozygote at p = 0.5: FIS = 1.
  expect_equal(1 - 0 / hs_biallelic(0.5, ho = 0, n = n), 1)
  ## A minor allele seen once must sit in a heterozygote, which forces Hs = Ho
  ## and so FIS = 0, whatever n.
  expect_equal(hs_biallelic(1 / (2 * n), ho = 1 / n, n = n), 1 / n)
})

test_that("the two forms of the 2n-corrected gene diversity agree, and match published Stacks output", {
  for (n in c(2L, 10L, 15L)) {
    for (p in c(0.1, 0.35, 0.5)) {
      counts <- c(round(2 * n * p), 2 * n - round(2 * n * p))
      expect_equal(gene_div_2n_counts(counts), hs_stacks_pi(counts[1] / sum(counts), n = n))
    }
  }
  ## A real Stacks summary row: Exp_Het 0.40000 at Num_Indv 2 is reported as
  ## Pi 0.53333, i.e. Exp_Het * 2n/(2n - 1).
  expect_equal(0.4 * (2 * 2) / (2 * 2 - 1), 0.53333, tolerance = 1e-5)
  ## Fewer than 2 gene copies: undefined, not 0.
  expect_true(is.na(gene_div_2n_counts(c(1, 0))))
})

test_that("hs_from_counts() is hs_nei_chesser() from allele counts", {
  counts <- c(6, 14)
  expect_equal(hs_from_counts(counts, ho = 0.35, n = 10),
               hs_nei_chesser(sum((counts / 20)^2), ho = 0.35, n = 10))
  expect_true(is.na(hs_from_counts(c(1, 0), ho = 0, n = 10)))   # < 2 gene copies
  expect_true(is.na(hs_from_counts(counts, ho = 0.35, n = 1)))  # < 2 individuals
})

test_that("fis_ratio_of_sums() ignores monomorphic loci and undefined values", {
  ho <- c(0.30, 0.10, 0)
  he <- c(0.42, 0.15, 0)
  expect_equal(fis_ratio_of_sums(ho, he), 1 - sum(ho) / sum(he))
  ## A monomorphic locus adds 0 to both sums, so appending any number of them
  ## changes nothing. This is the whole reason for a ratio of sums.
  expect_equal(fis_ratio_of_sums(c(ho, rep(0, 800)), c(he, rep(0, 800))),
               fis_ratio_of_sums(ho, he))
  ## A locus with no defined value is left out, not counted as a zero.
  expect_equal(fis_ratio_of_sums(c(ho, NA), c(he, NA)), fis_ratio_of_sums(ho, he))
  ## No expected heterozygosity anywhere: undefined.
  expect_true(is.na(fis_ratio_of_sums(c(0, 0), c(0, 0))))
})

test_that("autosomal_het() rescales per-record heterozygosity to per sequenced site", {
  expect_equal(autosomal_het(0.29, n_snps_used = 10348, n_sites_sequenced = 1e6),
               0.29 * 10348 / 1e6)
  ## One value per population, each with its own sites count.
  expect_equal(autosomal_het(c(0.3, 0.2), c(1000, 500), c(1e5, 2e5)),
               c(0.3 * 1000 / 1e5, 0.2 * 500 / 2e5))
  ## Fewer sequenced sites than variant records is impossible, so it is NA
  ## rather than a number larger than the per-record value.
  expect_true(is.na(autosomal_het(0.29, n_snps_used = 1000, n_sites_sequenced = 999)))
  expect_true(is.na(autosomal_het(0.29, n_snps_used = 1000, n_sites_sequenced = NA)))
})

test_that("p_sampled() and rare_richness() match the rarefaction formula by hand", {
  ## 2 copies of each of 2 alleles, drawing g = 2 of the 4 copies: an allele is
  ## missed only if both drawn copies are the other one, C(2,2)/C(4,2) = 1/6.
  expect_equal(p_sampled(c(2, 2), g = 2), c(5 / 6, 5 / 6))
  expect_equal(rare_richness(c(2, 2), g = 2), 5 / 3)
  ## An allele nobody carries is drawn with probability 0, so zero-padding the
  ## counts (as diversity_stats() does, to give every record the same columns)
  ## cannot change the answer.
  expect_equal(rare_richness(c(12, 5, 2, 1), g = 10),
               rare_richness(c(12, 5, 2, 1, 0, 0), g = 10))
  ## Drawing every copy must recover the alleles actually present.
  expect_equal(rare_richness(c(7, 3, 2, 1), g = 13), 4)
  ## Fewer copies than g: the draw is impossible, so NA, not 0.
  expect_true(is.na(rare_richness(c(3, 2), g = 10)))
  expect_true(all(is.na(p_sampled(c(3, 2), g = 10))))
  expect_true(all(is.na(p_sampled(c(3, 2), g = 0))))
})

test_that("rare_private() and rare_private_all() match a hand calculation", {
  ## Population 1 has both alleles (2 copies each), population 2 only allele 1
  ## (4 copies); g = 2. Pr(drawn) is 5/6 for each of population 1's alleles,
  ## and 1 and 0 for population 2's. Private to 1: 5/6 * (1 - 0) = 5/6 (only
  ## allele 2 can be private). Private to 2: 1 * (1 - 5/6) = 1/6.
  counts <- rbind(c(2, 2), c(4, 0))
  expect_equal(rare_private(counts, j = 1, g = 2), 5 / 6)
  expect_equal(rare_private(counts, j = 2, g = 2), 1 / 6)
  expect_equal(rare_private_all(counts, g = 2), c(5 / 6, 1 / 6))
  ## Three populations, each sampled completely: allele 1 is private to
  ## population 1, allele 2 is shared, allele 3 is private to population 3.
  three <- rbind(c(4, 4, 0), c(0, 8, 0), c(0, 4, 4))
  expect_equal(rare_private_all(three, g = 8), c(1, 0, 1))
  ## "Private" needs every population sampled at the same depth, so one
  ## population with too few copies makes it undefined for all of them.
  expect_true(is.na(rare_private(rbind(c(2, 2), c(1, 0)), j = 1, g = 2)))
})
