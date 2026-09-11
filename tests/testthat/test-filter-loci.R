fx <- function(name) test_path("fixtures", name)

## A small hand-built H, for tests that don't need a real VCF file. Mirrors
## exactly the shape read_stacks_vcf() returns, including sample-named
## columns on A1/A2 (several filters look samples up by name).
make_H <- function(A1, A2, locus = NULL, locus_raw = NULL, n_alleles = NULL,
                    alleles = NULL, samples = NULL, fields = NULL) {
  n_rec <- nrow(A1)
  if (is.null(samples)) samples <- colnames(A1)
  ## rbind(name1 = ..., name2 = ...) is a convenient way to label loci in a
  ## test's A1 -- pick those row names up as `locus` by default, if present.
  if (is.null(locus)) locus <- rownames(A1)
  dimnames(A1) <- dimnames(A2) <- list(NULL, samples)
  if (is.null(locus)) locus <- paste0("locus_", seq_len(n_rec))
  if (is.null(locus_raw)) locus_raw <- locus
  if (is.null(n_alleles)) n_alleles <- rep(2L, n_rec)
  if (is.null(alleles)) alleles <- replicate(n_rec, c("A", "C"), simplify = FALSE)
  list(A1 = A1, A2 = A2, locus = locus, locus_raw = locus_raw,
       alleles = alleles, n_alleles = n_alleles, samples = samples, fields = fields)
}

## ---------------------------------------------------------------------------
## .subset_H()
## ---------------------------------------------------------------------------

test_that(".subset_H() subsets every per-locus element consistently", {
  H <- make_H(A1 = rbind(c(1L, 1L), c(1L, 2L), c(2L, 2L)),
              A2 = rbind(c(1L, 1L), c(1L, 2L), c(2L, 2L)),
              samples = c("s1", "s2"))
  H2 <- RADdiversity:::.subset_H(H, c(TRUE, FALSE, TRUE))
  expect_equal(nrow(H2$A1), 2L)
  expect_equal(H2$locus, H$locus[c(1, 3)])
  expect_equal(H2$locus_raw, H$locus_raw[c(1, 3)])
  expect_equal(H2$n_alleles, H$n_alleles[c(1, 3)])
  expect_equal(H2$alleles, H$alleles[c(1, 3)])
  expect_equal(H2$samples, H$samples)  # samples untouched by locus subsetting
})

test_that(".subset_H() does not re-uniquify locus names", {
  H <- make_H(A1 = rbind(c(1L, 1L), c(1L, 1L)), A2 = rbind(c(1L, 1L), c(1L, 1L)),
              locus = c("dup", "dup.1"), samples = c("s1", "s2"))
  H2 <- RADdiversity:::.subset_H(H, c(TRUE, TRUE))
  expect_equal(H2$locus, c("dup", "dup.1"))
})

test_that(".subset_H() rejects a malformed `keep`", {
  H <- make_H(A1 = rbind(c(1L, 1L)), A2 = rbind(c(1L, 1L)), samples = c("s1", "s2"))
  expect_error(RADdiversity:::.subset_H(H, c(TRUE, TRUE)), "exactly")
  expect_error(RADdiversity:::.subset_H(H, c(TRUE, NA)), "missing")
  expect_error(RADdiversity:::.subset_H(H, 5L), "between 1 and")
})

## ---------------------------------------------------------------------------
## locus_allele_stats()
## ---------------------------------------------------------------------------

test_that("locus_allele_stats() distinguishes 'no data' (NA) from 'monomorphic' (0)", {
  H <- make_H(
    A1 = rbind(monomorphic = c(1L, 1L, 1L), no_data = c(NA, NA, NA)),
    A2 = rbind(monomorphic = c(1L, 1L, 1L), no_data = c(NA, NA, NA)),
    samples = c("s1", "s2", "s3")
  )
  st <- locus_allele_stats(H)
  expect_equal(st$maf, c(0, NA_real_))
  expect_equal(st$mac, c(0L, NA_integer_))
  expect_equal(st$n_called, c(6, 0))
})

test_that("locus_allele_stats() tabulates a multiallelic locus against its own n_alleles", {
  H <- make_H(
    A1 = rbind(c(1L, 2L, 3L)), A2 = rbind(c(1L, 2L, 3L)),
    n_alleles = 3L, samples = c("s1", "s2", "s3")
  )
  st <- locus_allele_stats(H)
  expect_equal(st$n_observed_alleles, 3L)
  expect_equal(st$n_called, 6)
  expect_equal(st$major_af, 2 / 6)  # each allele seen exactly twice -- a 3-way tie
})

## ---------------------------------------------------------------------------
## filter_maf() / filter_mac()
## ---------------------------------------------------------------------------

test_that("filter_maf()/filter_mac() keep a locus exactly AT the threshold", {
  H <- make_H(
    A1 = rbind(locus_1 = c(1L, 1L, 1L), locus_2 = c(1L, 2L, 2L)),
    A2 = rbind(locus_1 = c(1L, 1L, 1L), locus_2 = c(1L, 2L, 2L)),
    samples = c("s1", "s2", "s3")
  )
  ## locus_2: maf = 2/6 = 1/3, mac = 2 -- both thresholds set exactly there.
  expect_equal(filter_maf(H, min_maf = 1 / 3, verbose = FALSE)$locus, "locus_2")
  expect_equal(filter_mac(H, min_mac = 2, verbose = FALSE)$locus, "locus_2")
})

test_that("filter_maf()/filter_mac() drop a no-data locus and accept a precomputed `stats`", {
  H <- make_H(
    A1 = rbind(has_data = c(1L, 2L), no_data = c(NA, NA)),
    A2 = rbind(has_data = c(1L, 2L), no_data = c(NA, NA)),
    samples = c("s1", "s2")
  )
  st <- locus_allele_stats(H)
  expect_equal(filter_maf(H, min_maf = 0, stats = st, verbose = FALSE)$locus, "has_data")
  expect_equal(filter_mac(H, min_mac = 0, stats = st, verbose = FALSE)$locus, "has_data")
})

test_that("filter_maf()/filter_mac() reject a `stats` that doesn't match H, and error when nothing survives", {
  H <- make_H(A1 = rbind(c(1L, 1L)), A2 = rbind(c(1L, 1L)), samples = c("s1", "s2"))
  wrong_stats <- locus_allele_stats(H)[0, ]
  expect_error(filter_maf(H, min_maf = 0, stats = wrong_stats), "line up")
  expect_error(filter_mac(H, min_mac = 0, stats = wrong_stats), "line up")
  expect_error(filter_maf(H, min_maf = 0.9, verbose = FALSE), "No locus")
  expect_error(filter_mac(H, min_mac = 99, verbose = FALSE), "No locus")
})

## ---------------------------------------------------------------------------
## filter_call_rate()
## ---------------------------------------------------------------------------

test_that("filter_call_rate() with pops = NULL matches a plain pooled call rate", {
  H <- make_H(
    A1 = rbind(well_typed = c(1L, 1L, 1L, NA), poorly_typed = c(1L, NA, NA, NA)),
    A2 = rbind(well_typed = c(1L, 1L, 1L, NA), poorly_typed = c(1L, NA, NA, NA)),
    samples = c("s1", "s2", "s3", "s4")
  )
  expect_equal(filter_call_rate(H, min_call = 0.75, verbose = FALSE)$locus, "well_typed")
})

test_that("filter_call_rate() rule='all' vs 'any' genuinely differ under asymmetric coverage", {
  ## locus_1: population A fully typed, population B not at all.
  H <- make_H(
    A1 = rbind(c(1L, 1L, NA, NA)), A2 = rbind(c(1L, 1L, NA, NA)),
    samples = c("a1", "a2", "b1", "b2")
  )
  pops <- list(A = c("a1", "a2"), B = c("b1", "b2"))
  expect_length(filter_call_rate(H, min_call = 0.5, pops = pops, rule = "any", verbose = FALSE)$locus, 1L)
  expect_error(filter_call_rate(H, min_call = 0.5, pops = pops, rule = "all", verbose = FALSE), "No locus")
})

test_that("filter_call_rate() rejects an invalid `rule`", {
  H <- make_H(A1 = rbind(c(1L, 1L)), A2 = rbind(c(1L, 1L)), samples = c("s1", "s2"))
  expect_error(filter_call_rate(H, min_call = 0.5, rule = "some"), "rule must be")
})

## ---------------------------------------------------------------------------
## filter_max_het()
## ---------------------------------------------------------------------------

test_that("filter_max_het() drops an excess-het locus but keeps an all-missing one", {
  H <- make_H(
    A1 = rbind(low_het = c(1L, 1L, 1L, 1L), high_het = c(1L, 2L, 1L, 2L), no_data = c(NA, NA, NA, NA)),
    A2 = rbind(low_het = c(1L, 1L, 1L, 1L), high_het = c(2L, 1L, 2L, 1L), no_data = c(NA, NA, NA, NA)),
    samples = c("s1", "s2", "s3", "s4")
  )
  kept <- filter_max_het(H, max_ho = 0.5, verbose = FALSE)$locus
  expect_setequal(kept, c("low_het", "no_data"))
})

test_that("filter_max_het() keeps a locus exactly AT the threshold", {
  H <- make_H(A1 = rbind(c(1L, 1L, 1L, 1L)), A2 = rbind(c(1L, 2L, 1L, 1L)),
              samples = c("s1", "s2", "s3", "s4"))  # Ho = 0.25 exactly
  expect_length(filter_max_het(H, max_ho = 0.25, verbose = FALSE)$locus, 1L)
})

## ---------------------------------------------------------------------------
## filter_thin_one_snp()
## ---------------------------------------------------------------------------

test_that("filter_thin_one_snp() keeps one record per locus_raw, method='first' is deterministic", {
  H <- make_H(
    A1 = rbind(c(1L, 1L), c(1L, 1L), c(1L, 1L)),
    A2 = rbind(c(1L, 1L), c(1L, 1L), c(1L, 1L)),
    locus = c("tag1_a", "tag1_b", "tag2"), locus_raw = c("tag1", "tag1", "tag2"),
    samples = c("s1", "s2")
  )
  out1 <- filter_thin_one_snp(H, method = "first", verbose = FALSE)
  out2 <- filter_thin_one_snp(H, method = "first", verbose = FALSE)
  expect_equal(out1$locus, c("tag1_a", "tag2"))  # first record of tag1, plus tag2, in original order
  expect_equal(out1$locus, out2$locus)
})

test_that("filter_thin_one_snp() method='random' never mis-samples a singleton group", {
  ## A locus_raw group of size 1 must always return THAT record -- sample(i, 1)
  ## would instead draw a random integer in 1:i for a numeric i, which this
  ## regression test guards against.
  H <- make_H(A1 = rbind(c(1L, 1L)), A2 = rbind(c(1L, 1L)),
              locus = "solo", locus_raw = "solo", samples = c("s1", "s2"))
  out <- filter_thin_one_snp(H, method = "random", seed = 1, verbose = FALSE)
  expect_equal(out$locus, "solo")
})

test_that("filter_thin_one_snp() is a documented no-op when locus_raw is already unique", {
  H <- make_H(A1 = rbind(c(1L, 1L), c(1L, 1L)), A2 = rbind(c(1L, 1L), c(1L, 1L)),
              samples = c("s1", "s2"))
  expect_message(out <- filter_thin_one_snp(H, verbose = TRUE), "nothing to thin")
  expect_equal(out$locus, H$locus)
})

test_that("filter_thin_one_snp() restores the caller's RNG state", {
  H <- make_H(
    A1 = rbind(c(1L, 1L), c(1L, 1L), c(1L, 1L)),
    A2 = rbind(c(1L, 1L), c(1L, 1L), c(1L, 1L)),
    locus = c("tag1_a", "tag1_b", "tag2"), locus_raw = c("tag1", "tag1", "tag2"),
    samples = c("s1", "s2")
  )
  set.seed(123); expected <- runif(1)
  set.seed(123); invisible(filter_thin_one_snp(H, method = "random", seed = 999, verbose = FALSE))
  expect_equal(runif(1), expected)
})

test_that("filter_thin_one_snp() rejects an inexact `method`", {
  H <- make_H(A1 = rbind(c(1L, 1L)), A2 = rbind(c(1L, 1L)), samples = c("s1", "s2"))
  expect_error(filter_thin_one_snp(H, method = "firs"), "method must be")
})

## ---------------------------------------------------------------------------
## filter_low_conf_alt() / low_conf_alt_sensitivity() -- uses the small_ad
## fixture, which was built specifically to exercise every documented edge
## case (see tests/testthat/fixtures/small_ad.haps.vcf).
## ---------------------------------------------------------------------------

test_that("filter_low_conf_alt() flags each fixture row exactly as documented", {
  H <- suppressMessages(read_stacks_vcf(fx("small_ad.haps.vcf")))
  res <- filter_low_conf_alt(H, min_alt_reads = 2, mode = "mask", verbose = FALSE)
  ls <- res$locus_summary
  row <- function(loc) ls[ls$locus == loc, ]

  expect_equal(row("locus_ad1")$n_flagged, 1L)   # 1/1, AD_alt = 2 (not double-counted) -- flagged
  expect_equal(row("locus_ad2")$n_flagged, 1L)   # 1/2 multiallelic, AD_alt = 1+1 = 2 -- flagged
  expect_equal(row("locus_ad3")$n_alt_calls, 1L) # AD dropped via GT:DP:AD ordering
  expect_equal(row("locus_ad3")$n_flagged, 0L)   # unusable calls are never flagged
  expect_equal(row("locus_ad4")$n_flagged, 0L)   # AD = "."
  expect_equal(row("locus_ad5")$n_alt_calls, 0L) # homozygous REF is never an ALT call, however low its AD
  expect_equal(row("locus_ad6")$n_flagged, 0L)   # no AD anywhere in FORMAT
  expect_equal(row("locus_ad7")$n_flagged, 0L)   # AD_alt = 4 > 2, via DP-absent fallback to sum(AD)
  expect_equal(row("locus_ad8")$n_alt_calls, 4L)
  expect_equal(row("locus_ad8")$n_flagged, 1L)
  expect_equal(row("locus_ad8")$flagged_fraction, 0.25)
  expect_equal(row("locus_ad9")$n_alt_calls, 1L) # the ./. sample there is never counted at all
})

test_that("filter_low_conf_alt() mode='mask' changes only the flagged cells", {
  H <- suppressMessages(read_stacks_vcf(fx("small_ad.haps.vcf")))
  res <- filter_low_conf_alt(H, min_alt_reads = 2, mode = "mask", verbose = FALSE)
  j <- match("locus_ad1", H$locus)
  expect_true(is.na(res$H$A1[j, "popA_1"]))   # the flagged cell is now missing
  expect_true(is.na(res$H$A2[j, "popA_1"]))
  expect_equal(nrow(res$H$A1), nrow(H$A1))    # mask never removes a record
  expect_true(all(res$locus_summary$kept))
})

test_that("filter_low_conf_alt() mode='drop' respects drop_frac at its exact boundary", {
  H <- suppressMessages(read_stacks_vcf(fx("small_ad.haps.vcf")))
  in_result <- function(drop_frac) {
    r <- filter_low_conf_alt(H, min_alt_reads = 2, mode = "drop", drop_frac = drop_frac, verbose = FALSE)
    "locus_ad8" %in% r$H$locus
  }
  expect_false(in_result(0))
  expect_false(in_result(0.2))
  expect_true(in_result(0.25))  # flagged_fraction is exactly 0.25 -- boundary is inclusive
  expect_true(in_result(0.3))
})

test_that("locus_summary$kept agrees with the drop rule for every locus", {
  H <- suppressMessages(read_stacks_vcf(fx("small_ad.haps.vcf")))
  res <- filter_low_conf_alt(H, min_alt_reads = 2, mode = "drop", drop_frac = 0.2, verbose = FALSE)
  expect_equal(res$locus_summary$kept, res$locus_summary$flagged_fraction <= 0.2 + 1e-9)
})

test_that("masking then dropping does not reproduce the two independent Python-script views", {
  H <- suppressMessages(read_stacks_vcf(fx("small_ad.haps.vcf")))
  masked <- filter_low_conf_alt(H, min_alt_reads = 2, mode = "mask", verbose = FALSE)$H
  ## After masking, the flagged calls are gone -- a further drop pass, even
  ## at drop_frac = 0, has nothing left to drop.
  redrop <- filter_low_conf_alt(masked, min_alt_reads = 2, mode = "drop", drop_frac = 0, verbose = FALSE)
  expect_equal(nrow(redrop$H$A1), nrow(masked$A1))
})

test_that("low_conf_alt_sensitivity() is monotonic and agrees with filter_low_conf_alt() at its own threshold", {
  H <- suppressMessages(read_stacks_vcf(fx("small_ad.haps.vcf")))
  sens <- low_conf_alt_sensitivity(H, thresholds = c(1, 2, 3, 4, 5, 10))
  expect_true(all(diff(sens$n_flagged) >= 0))
  res <- filter_low_conf_alt(H, min_alt_reads = 2, mode = "mask", verbose = FALSE)
  expect_equal(sens$n_flagged[sens$threshold == 2], sum(res$locus_summary$n_flagged))
})
