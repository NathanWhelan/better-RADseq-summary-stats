fx <- function(name) test_path("fixtures", name)

test_that("read_stacks_vcf() parses the small fixture", {
  H <- suppressMessages(read_stacks_vcf(fx("small.haps.vcf")))
  expect_setequal(H$samples,
                   c(paste0("popA_", 1:4), paste0("popB_", 1:3)))
  expect_equal(nrow(H$A1), 80L)
  expect_equal(ncol(H$A1), 7L)
  expect_true(all(H$n_alleles %in% c(2L, 3L)))
})

test_that("read_popmap() splits samples by population", {
  H <- suppressMessages(read_stacks_vcf(fx("small.haps.vcf"), verbose = FALSE))
  pops <- suppressMessages(read_popmap(fx("small_popmap.tsv"), H$samples))
  expect_equal(names(pops), c("popA", "popB"))
  expect_equal(lengths(pops), c(popA = 4L, popB = 3L))
})

test_that("read_popmap() keeps IDs as text (leading zeros and T/F survive)", {
  ## Regression: read.delim()'s type guessing turned "001" into 1 (no longer
  ## matching the VCF) and a population named "T" into TRUE.
  pm <- tempfile(fileext = ".tsv")
  on.exit(unlink(pm), add = TRUE)
  writeLines(c("001\tT", "002\tT", "010\tF"), pm)
  pops <- read_popmap(pm, c("001", "002", "010"), verbose = FALSE)
  expect_identical(pops, list(T = c("001", "002"), F = "010"))
})

test_that(".resolve_H() passes a path through to read_stacks_vcf() and an H list straight through", {
  H_from_path <- suppressMessages(RADdiversity:::.resolve_H(fx("small.haps.vcf"), verbose = FALSE))
  expect_equal(H_from_path$samples, suppressMessages(read_stacks_vcf(fx("small.haps.vcf"), verbose = FALSE))$samples)

  H <- suppressMessages(read_stacks_vcf(fx("small.haps.vcf"), verbose = FALSE))
  expect_identical(RADdiversity:::.resolve_H(H, verbose = FALSE), H)
})

test_that(".resolve_H() rejects something that isn't a path or a proper H list", {
  expect_error(RADdiversity:::.resolve_H(42), "path to a VCF file")
  expect_error(RADdiversity:::.resolve_H(list(A1 = matrix(1))), "missing element")
  bad <- list(A1 = matrix(1), A2 = matrix(1:2), locus = "x", locus_raw = "x",
              alleles = list("A"), n_alleles = 1L, samples = "s1")
  expect_error(RADdiversity:::.resolve_H(bad), "same size")
})

test_that("diversity_stats()/het_between_pops() accept a pre-parsed H list in place of a path", {
  outdir <- tempfile("raddiversity-test-")
  dir.create(outdir)
  on.exit(unlink(outdir, recursive = TRUE), add = TRUE)
  H <- suppressMessages(read_stacks_vcf(fx("small.haps.vcf")))

  expect_error(
    suppressMessages(diversity_stats(H, fx("small_popmap.tsv"), g = 4, nboot = 10, outdir = outdir)),
    "stem"
  )
  res <- suppressMessages(capture.output(
    result <- diversity_stats(H, fx("small_popmap.tsv"), g = 4, nboot = 10, outdir = outdir, stem = "test")
  ))
  expect_named(result, c("per_population", "richness", "autosomal"))
  expect_true(file.exists(file.path(outdir, "diversity_per_population.test.tsv")))

  expect_error(
    suppressMessages(het_between_pops(H, fx("small_popmap.tsv"), min_call = 0.5, outdir = outdir)),
    "stem"
  )
  het_res <- suppressMessages(capture.output(
    het_result <- het_between_pops(H, fx("small_popmap.tsv"), min_call = 0.5,
                                   outdir = outdir, stem = "test")
  ))
  expect_named(het_result, c("individual_heterozygosity", "pairwise_tests"))
  expect_true(file.exists(file.path(outdir, "individual_heterozygosity.test.tsv")))
})
