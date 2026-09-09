test_that("diversity_core_selftest() reports all checks passing", {
  invisible(capture.output(ok <- diversity_core_selftest()))
  expect_true(ok)
})

test_that("het_between_pops_selftest() runs without error and returns TRUE", {
  invisible(capture.output(ok <- het_between_pops_selftest()))
  expect_true(ok)
})
