## The self-tests are ~25 s of Monte Carlo between them -- too slow for CRAN's
## test budget -- so they are skipped there. devtools::test() and the GitHub
## Actions check (NOT_CRAN = "true") still run them.
test_that("diversity_core_selftest() reports all checks passing", {
  skip_on_cran()
  invisible(capture.output(ok <- diversity_core_selftest()))
  expect_true(ok)
})

test_that("het_between_pops_selftest() runs without error and returns TRUE", {
  skip_on_cran()
  invisible(capture.output(ok <- het_between_pops_selftest()))
  expect_true(ok)
})
