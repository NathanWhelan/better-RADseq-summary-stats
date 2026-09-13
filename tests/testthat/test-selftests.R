## The self-tests are ~25 s of Monte Carlo between them -- too slow for CRAN's
## test budget -- so they are skipped there. devtools::test() and the GitHub
## Actions check (NOT_CRAN = "true") still run them.
test_that("diversity_core_selftest() reports all checks passing", {
  skip_on_cran()
  result <- diversity_core_selftest(verbose = FALSE)
  expect_s3_class(result, "data.frame")
  expect_true(all(result$pass))
})

test_that("het_between_pops_selftest() holds the nominal rate for the individual-level tests", {
  skip_on_cran()
  result <- het_between_pops_selftest(verbose = FALSE)
  individual_level <- result[result$method %in% c("welch", "wilcoxon"), ]
  expect_true(all(individual_level$rejection_rate < 0.10))
  expect_silent(het_between_pops_selftest(verbose = FALSE))
})
