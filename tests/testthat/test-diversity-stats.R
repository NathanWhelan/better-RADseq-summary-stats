fx <- function(name) test_path("fixtures", name)

test_that("diversity_stats() returns the expected structure on a small fixture", {
  outdir <- tempfile("raddiversity-test-")
  dir.create(outdir)
  on.exit(unlink(outdir, recursive = TRUE), add = TRUE)
  res <- suppressMessages(capture.output(
    result <- diversity_stats(fx("small.haps.vcf"), fx("small_popmap.tsv"),
                               g = 4, nboot = 50, outdir = outdir)
  ))
  expect_named(result, c("per_population", "richness", "autosomal"))
  expect_setequal(result$per_population$population, c("popA", "popB"))
  expect_true(all(c("Ho", "He", "Fis", "pct_poly") %in% names(result$per_population)))
  expect_true(all(c("Ar", "privAr", "priv_total") %in% names(result$richness)))
  expect_true(file.exists(file.path(outdir, "diversity_per_population.haps.tsv")))
  expect_true(file.exists(file.path(outdir, "diversity_richness.haps.tsv")))
})

test_that("diversity_stats(hierfstat_check = TRUE) agrees with hierfstat", {
  skip_if_not_installed("hierfstat")
  outdir <- tempfile("raddiversity-test-")
  dir.create(outdir)
  on.exit(unlink(outdir, recursive = TRUE), add = TRUE)
  msgs <- capture_messages(expect_no_warning(invisible(capture.output(
    diversity_stats(fx("small.haps.vcf"), fx("small_popmap.tsv"), g = 4, nboot = 0,
                    outdir = outdir, hierfstat_check = TRUE)))))
  expect_true(any(grepl("cross-check vs hierfstat::basic.stats", msgs)))
  expect_true(any(grepl("cross-check vs hierfstat::allelic.richness", msgs)))
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

test_that("diversity_stats() sites= accepts a scalar, a named vector, a data frame, or a sumstats_summary path -- all equivalent when the value is the same for every population", {
  outdir <- tempfile("raddiversity-test-")
  dir.create(outdir)
  on.exit(unlink(outdir, recursive = TRUE), add = TRUE)
  go <- function(sites) {
    out <- NULL
    invisible(capture.output(suppressMessages(
      out <- diversity_stats(fx("small.snps.vcf"), fx("small_popmap.tsv"),
                              g = 4, nboot = 0, outdir = outdir, sites = sites)
    )))
    out
  }

  r_scalar <- go(100000)
  r_vector <- go(c(popA = 100000, popB = 100000))
  r_df     <- go(data.frame(population = c("popA", "popB"), sites = c(100000, 100000)))
  r_path   <- go(fx("sumstats_summary_small.tsv"))

  expect_false(is.null(r_scalar$autosomal))
  expect_identical(r_scalar$autosomal, r_vector$autosomal)
  expect_identical(r_scalar$autosomal, r_df$autosomal)
  expect_identical(r_scalar$autosomal, r_path$autosomal)
  expect_true(all(c("sites_used", "Ho_autosomal", "He_autosomal") %in% names(r_scalar$autosomal)))
  expect_true(file.exists(file.path(outdir, "diversity_autosomal.snps.tsv")))
})

test_that("diversity_stats() sites= as a per-population vector gives each population its own denominator", {
  outdir <- tempfile("raddiversity-test-")
  dir.create(outdir)
  on.exit(unlink(outdir, recursive = TRUE), add = TRUE)
  out <- NULL
  invisible(capture.output(suppressMessages(
    out <- diversity_stats(fx("small.snps.vcf"), fx("small_popmap.tsv"), g = 4, nboot = 0,
                            outdir = outdir, sites = c(popA = 100000, popB = 50000))
  )))
  aut <- out$autosomal
  n_rec <- 15L  # tests/testthat/fixtures/small.snps.vcf has 15 records
  ## Each population's own He (per marker, from per_population) scaled by ITS
  ## OWN sites value -- not a shared scalar, which would give both
  ## populations the same denominator regardless of their true sites.
  for (p in c("popA", "popB")) {
    he_marker <- out$per_population$He[out$per_population$population == p]
    sites_p   <- c(popA = 100000, popB = 50000)[[p]]
    ## he_marker (from per_population) is itself already rounded to 4 decimal
    ## places, so a loose relative tolerance absorbs that rounding rather than
    ## re-deriving diversity_stats()'s full-precision internal value here.
    expect_equal(aut$He_autosomal[aut$population == p],
                 signif(he_marker * n_rec / sites_p, 4), tolerance = 1e-3)
  }
  ## And the two populations' own values are NOT equal to each other (unlike
  ## the equal-sites case above) -- confirms the per-population denominator
  ## actually took effect rather than silently falling back to one shared value.
  expect_false(isTRUE(all.equal(aut$He_autosomal[aut$population == "popA"],
                                 aut$He_autosomal[aut$population == "popB"])))
})

test_that("diversity_stats() sites= flags an implausible per-population value as NA rather than erroring the whole run", {
  outdir <- tempfile("raddiversity-test-")
  dir.create(outdir)
  on.exit(unlink(outdir, recursive = TRUE), add = TRUE)
  out <- NULL
  invisible(capture.output(suppressMessages(
    out <- diversity_stats(fx("small.snps.vcf"), fx("small_popmap.tsv"), g = 4, nboot = 0,
                            outdir = outdir, sites = c(popA = 100000, popB = 5))
  )))
  aut <- out$autosomal
  expect_false(is.na(aut$He_autosomal[aut$population == "popA"]))
  expect_true(is.na(aut$He_autosomal[aut$population == "popB"]))
})

test_that("diversity_stats() sites= errors when a population has no value, naming it", {
  expect_error(
    suppressMessages(diversity_stats(fx("small.snps.vcf"), fx("small_popmap.tsv"),
                                      g = 4, nboot = 0, sites = c(popA = 100000))),
    "popB"
  )
})

test_that("diversity_stats() sites= rejects a bad shape before the VCF is even read", {
  expect_error(diversity_stats("nonexistent.vcf", fx("small_popmap.tsv"), g = 4, sites = -5),
               "non-negative")
  expect_error(diversity_stats("nonexistent.vcf", fx("small_popmap.tsv"), g = 4,
                                sites = "no/such/file.tsv"),
               "not found")
})

test_that("diversity_stats() ignores sites= entirely (autosomal NULL) on a haplotype VCF", {
  outdir <- tempfile("raddiversity-test-")
  dir.create(outdir)
  on.exit(unlink(outdir, recursive = TRUE), add = TRUE)
  out <- NULL
  invisible(capture.output(suppressMessages(
    out <- diversity_stats(fx("small.haps.vcf"), fx("small_popmap.tsv"), g = 4, nboot = 0,
                            outdir = outdir, sites = 100000)
  )))
  expect_null(out$autosomal)
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
  expect_true(file.exists(file.path(outdir, "individual_heterozygosity.haps.tsv")))
  expect_true(file.exists(file.path(outdir, "het_between_pops_tests.haps.tsv")))
})

test_that("het_between_pops() runs on the SNP and haplotype VCFs into one outdir without overwriting", {
  ## Regression: the output names used to be fixed, so the second run
  ## silently replaced the first run's files.
  outdir <- tempfile("raddiversity-test-")
  dir.create(outdir)
  on.exit(unlink(outdir, recursive = TRUE), add = TRUE)
  ## The same fixture under the two Stacks file names -- only the name matters.
  vcfs <- file.path(outdir, c("populations.snps.vcf", "populations.haps.vcf"))
  file.copy(fx("small.haps.vcf"), vcfs)
  for (v in vcfs)
    invisible(capture.output(suppressMessages(
      het_between_pops(v, fx("small_popmap.tsv"), min_call = 0.5, outdir = outdir))))
  expect_true(all(file.exists(file.path(outdir, c(
    "individual_heterozygosity.snps.tsv", "individual_heterozygosity.haps.tsv",
    "het_between_pops_tests.snps.tsv", "het_between_pops_tests.haps.tsv")))))
})
