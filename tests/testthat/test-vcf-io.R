fx <- function(name) test_path("fixtures", name)

test_that("read_haps_vcf() parses the small fixture", {
  H <- suppressMessages(read_haps_vcf(fx("small.haps.vcf")))
  expect_setequal(H$samples,
                   c(paste0("popA_", 1:4), paste0("popB_", 1:3)))
  expect_equal(nrow(H$A1), 80L)
  expect_equal(ncol(H$A1), 7L)
  expect_true(all(H$n_alleles %in% c(2L, 3L)))
})

test_that("read_popmap() splits samples by population", {
  H <- suppressMessages(read_haps_vcf(fx("small.haps.vcf"), verbose = FALSE))
  pops <- suppressMessages(read_popmap(fx("small_popmap.tsv"), H$samples))
  expect_equal(names(pops), c("popA", "popB"))
  expect_equal(lengths(pops), c(popA = 4L, popB = 3L))
})
