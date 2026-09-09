fx <- function(name) test_path("fixtures", name)

test_that("diversity_stats() returns the expected structure on a small fixture", {
  outdir <- tempfile("raddiversity-test-")
  dir.create(outdir)
  on.exit(unlink(outdir, recursive = TRUE), add = TRUE)
  res <- suppressMessages(capture.output(
    result <- diversity_stats(fx("small.haps.vcf"), fx("small_popmap.tsv"),
                               g = 4, nboot = 50, outdir = outdir)
  ))
  expect_named(result, c("per_population", "richness"))
  expect_setequal(result$per_population$population, c("popA", "popB"))
  expect_true(all(c("Ho", "He", "Fis", "pct_poly") %in% names(result$per_population)))
  expect_true(all(c("Ar", "privAr", "priv_total") %in% names(result$richness)))
  expect_true(file.exists(file.path(outdir, "diversity_per_population.haps.tsv")))
  expect_true(file.exists(file.path(outdir, "diversity_richness.haps.tsv")))
})

test_that("diversity_stats() requires g and validates its arguments", {
  expect_error(diversity_stats(fx("small.haps.vcf"), fx("small_popmap.tsv")),
               "missing")
  expect_error(
    suppressMessages(diversity_stats(fx("small.haps.vcf"), fx("small_popmap.tsv"),
                                      g = 4, min_n = 1)),
    "min_n"
  )
})

test_that("diversity_stats() rejects an inexact boot value instead of partial-matching it", {
  expect_error(
    diversity_stats(fx("small.haps.vcf"), fx("small_popmap.tsv"),
                     g = 4, nboot = 5, boot = "individual"),
    "boot must be one of"
  )
})

test_that("diversity_stats() and het_between_pops() restore the caller's RNG state", {
  outdir <- tempfile("raddiversity-test-")
  dir.create(outdir)
  on.exit(unlink(outdir, recursive = TRUE), add = TRUE)

  set.seed(999)
  expected <- runif(1)

  set.seed(999)
  invisible(capture.output(suppressMessages(
    diversity_stats(fx("small.haps.vcf"), fx("small_popmap.tsv"),
                     g = 4, nboot = 20, outdir = outdir)
  )))
  expect_equal(runif(1), expected)

  set.seed(999)
  invisible(capture.output(suppressMessages(
    het_between_pops(fx("small.haps.vcf"), fx("small_popmap.tsv"),
                      min_call = 0.5, outdir = outdir)
  )))
  expect_equal(runif(1), expected)
})

test_that("het_between_pops() returns the expected structure on a small fixture", {
  outdir <- tempfile("raddiversity-test-")
  dir.create(outdir)
  on.exit(unlink(outdir, recursive = TRUE), add = TRUE)
  invisible(capture.output(
    result <- suppressMessages(
      het_between_pops(fx("small.haps.vcf"), fx("small_popmap.tsv"),
                        min_call = 0.5, outdir = outdir)
    )
  ))
  expect_named(result, c("individual_heterozygosity", "pairwise_tests"))
  expect_equal(nrow(result$individual_heterozygosity), 7L)
  expect_true(all(c("p_welch", "p_wilcox", "hedges_g") %in% names(result$pairwise_tests)))
  expect_true(file.exists(file.path(outdir, "individual_heterozygosity.tsv")))
  expect_true(file.exists(file.path(outdir, "het_between_pops_tests.tsv")))
})
