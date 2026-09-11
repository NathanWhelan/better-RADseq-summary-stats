## Same make_H() helper convention as the other test files in this package.
make_H <- function(A1, A2, locus = NULL, locus_raw = NULL, n_alleles = NULL,
                    alleles = NULL, samples = NULL) {
  n_rec <- nrow(A1)
  if (is.null(samples)) samples <- colnames(A1)
  if (is.null(locus)) locus <- rownames(A1)
  dimnames(A1) <- dimnames(A2) <- list(NULL, samples)
  if (is.null(locus)) locus <- paste0("locus_", seq_len(n_rec))
  if (is.null(locus_raw)) locus_raw <- locus
  if (is.null(n_alleles)) n_alleles <- rep(2L, n_rec)
  if (is.null(alleles)) alleles <- replicate(n_rec, c("A", "C"), simplify = FALSE)
  list(A1 = A1, A2 = A2, locus = locus, locus_raw = locus_raw,
       alleles = alleles, n_alleles = n_alleles, samples = samples)
}

## A dataset big enough (40 loci) to clear the "< 30 shared loci -> NA"
## floor, with one PAIR of samples constructed as genotype-identical (s1/s2)
## and the rest built from independent random draws so they read as
## effectively unrelated.
make_kinship_H <- function(n_loc = 40, seed = 1) {
  set.seed(seed)
  samples <- c("s1", "s2", "s3", "s4")
  A1 <- matrix(sample(1:2, n_loc * 4, replace = TRUE), n_loc, 4)
  A2 <- matrix(sample(1:2, n_loc * 4, replace = TRUE), n_loc, 4)
  ## s2 becomes an exact genotype copy of s1.
  A1[, 2] <- A1[, 1]; A2[, 2] <- A2[, 1]
  make_H(A1 = A1, A2 = A2, samples = samples)
}

test_that("method = \"king\" needs no dependency and runs by default", {
  H <- make_kinship_H()
  res <- suppressMessages(kinship_check(H, verbose = FALSE))
  expect_named(res, c("pairwise", "flagged_pairs"))
  expect_true(all(c("sample1", "sample2", "kinship", "n_loci_used", "low_confidence") %in% names(res$pairwise)))
})

test_that("an identical pair reads near 0.5 under both methods, on a distinct pair", {
  H <- make_kinship_H()
  king <- suppressMessages(kinship_check(H, method = "king", verbose = FALSE))$pairwise
  row <- king[king$sample1 == "s1" & king$sample2 == "s2", ]
  expect_equal(row$kinship, 0.5, tolerance = 1e-8)

  testthat::skip_if_not_installed("hierfstat")
  beta <- suppressMessages(kinship_check(H, method = "beta", verbose = FALSE))$pairwise
  rowb <- beta[beta$sample1 == "s1" & beta$sample2 == "s2", ]
  expect_gt(rowb$kinship, 0.4)  # beta.dosage()'s own scale isn't guaranteed exactly 0.5
})

test_that("method = \"beta\" stops with an actionable message when hierfstat is unavailable", {
  H <- make_kinship_H()
  testthat::with_mocked_bindings(
    .hierfstat_available = function() FALSE,
    expect_error(kinship_check(H, method = "beta", verbose = FALSE), "hierfstat"),
    .package = "RADdiversity"
  )
})

test_that("multiallelic loci are excluded (with a message) before either method runs", {
  H <- make_kinship_H()
  ## Make locus 1 triallelic -- both methods should silently drop it and say so.
  H$A1[1, 1] <- 3L; H$n_alleles[1] <- 3L; H$alleles[[1]] <- c("A", "C", "T")
  expect_message(kinship_check(H, verbose = TRUE), "excluding 1 of")
})

test_that("n_loci_used correctly restricts to the shared both-called locus set", {
  H <- make_kinship_H(n_loc = 50)
  ## Wipe out 10 DIFFERENT loci for s1 and s3 (no overlap), so their shared
  ## both-called count is 50 - 10 - 10 = 30, distinct from every other pair's
  ## full 50.
  H$A1[1:10, "s1"] <- NA; H$A2[1:10, "s1"] <- NA
  H$A1[11:20, "s3"] <- NA; H$A2[11:20, "s3"] <- NA
  res <- suppressMessages(kinship_check(H, verbose = FALSE))$pairwise
  row13 <- res[res$sample1 == "s1" & res$sample2 == "s3", ]
  row24 <- res[res$sample1 == "s2" & res$sample2 == "s4", ]  # untouched pair
  expect_equal(row13$n_loci_used, 30L)
  expect_equal(row24$n_loci_used, 50L)
})

test_that("a pair below the low-confidence floor is flagged/NA rather than silently reported", {
  H <- make_kinship_H(n_loc = 50)
  ## Leave only 20 shared loci between s1 and s4 -- below the 30-locus floor,
  ## so kinship should come back NA (not a number) for that pair specifically.
  H$A1[1:30, "s1"] <- NA; H$A2[1:30, "s1"] <- NA
  res <- suppressMessages(kinship_check(H, verbose = FALSE))$pairwise
  row14 <- res[res$sample1 == "s1" & res$sample2 == "s4", ]
  expect_equal(row14$n_loci_used, 20L)
  expect_true(is.na(row14$kinship))
  expect_true(row14$low_confidence)
  ## An NA-kinship pair must never appear in flagged_pairs, even if it would
  ## have exceeded threshold (it can't be compared to a threshold at all).
  flagged <- suppressMessages(kinship_check(H, verbose = FALSE))$flagged_pairs
  expect_false(any(is.na(flagged$kinship)))
})

test_that("threshold controls which pairs land in flagged_pairs", {
  H <- make_kinship_H()
  res_strict <- suppressMessages(kinship_check(H, threshold = 0.9, verbose = FALSE))
  res_mid    <- suppressMessages(kinship_check(H, threshold = 0.4, verbose = FALSE))
  res_loose  <- suppressMessages(kinship_check(H, threshold = -1, verbose = FALSE))
  expect_true(nrow(res_strict$flagged_pairs) < nrow(res_mid$flagged_pairs))
  expect_true(nrow(res_mid$flagged_pairs) < nrow(res_loose$flagged_pairs))
  ## s1/s2 (constructed identical, kinship 0.5) must be flagged once
  ## threshold drops below 0.5, but not above it.
  expect_false(any(res_strict$flagged_pairs$sample1 == "s1" & res_strict$flagged_pairs$sample2 == "s2"))
  expect_true(any(res_mid$flagged_pairs$sample1 == "s1" & res_mid$flagged_pairs$sample2 == "s2"))
})

test_that("outdir writes kinship_pairwise.tsv", {
  H <- make_kinship_H()
  outdir <- tempfile("kinship-")
  suppressMessages(kinship_check(H, outdir = outdir, verbose = FALSE))
  expect_true(file.exists(file.path(outdir, "kinship_pairwise.tsv")))
})

test_that("a pair with no heterozygous locus between them gets NA, not NaN/-Inf", {
  ## Regression: KING's denominator is the pair's heterozygous-locus count;
  ## when it is 0 the ratio is 0/0 or -x/0, which is no information at all.
  H <- make_kinship_H(n_loc = 40)
  H$A1[, "s1"] <- 1L; H$A2[, "s1"] <- 1L    # s1 homozygous REF everywhere
  H$A1[, "s2"] <- 2L; H$A2[, "s2"] <- 2L    # s2 homozygous ALT everywhere
  res <- suppressMessages(kinship_check(H, verbose = FALSE))$pairwise
  k12 <- res$kinship[res$sample1 == "s1" & res$sample2 == "s2"]
  expect_true(is.na(k12))
  expect_false(is.nan(k12))
})

test_that("method must be exactly \"king\" or \"beta\"", {
  H <- make_kinship_H()
  expect_error(kinship_check(H, method = "kinship", verbose = FALSE), "method must be")
})
