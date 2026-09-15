## The self-tests are ~25 s of Monte Carlo between them -- too slow for CRAN's
## test budget -- so they are skipped there. devtools::test() and the GitHub
## Actions check (NOT_CRAN = "true") still run them.
test_that("diversity_core_selftest() reports all checks passing", {
  skip_on_cran()
  result <- diversity_core_selftest(verbose = FALSE)
  expect_s3_class(result, "data.frame")
  expect_true(all(result$pass))
})

test_that("het_between_pops_selftest(): the combined test holds its rate where the others fail", {
  skip_on_cran()
  expect_silent(result <- het_between_pops_selftest(verbose = FALSE))
  rate <- function(method, fst, sd_F)
    result$rejection_rate[result$method == method & result$fst == fst & result$sd_F == sd_F]
  expect_true(all(result$rejection_rate[result$method == "combined"] < 0.10))
  ## Loci alone fail when individuals differ in inbreeding ...
  expect_gt(rate("locus_bootstrap", 0, 0.10), 0.25)
  ## ... and individuals alone when populations are differentiated.
  expect_gt(rate("welch", 0.10, 0), 0.10)
})
