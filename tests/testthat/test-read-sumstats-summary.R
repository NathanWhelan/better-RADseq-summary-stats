fx <- function(name) test_path("fixtures", name)

test_that("read_sumstats_summary() parses both blocks of the small fixture", {
  x <- read_sumstats_summary(fx("sumstats_summary_small.tsv"))
  expect_named(x, c("variant_positions", "all_positions"))

  vp <- x$variant_positions
  expect_equal(vp$population, c("popA", "popB"))
  expect_equal(vp$private, c(2, 1))
  expect_equal(vp$obs_het, c(0.20000, 0.18000))
  expect_equal(vp$pi, c(0.19500, 0.17500))
  expect_true(all(c("obs_het_var", "obs_het_se", "fis", "fis_var", "fis_se") %in% names(vp)))
  ## "All positions" is a different block with the SAME statistic names
  ## meaning different things (per-sequenced-site, not per-variant-record) --
  ## the two must never be merged into one data frame.
  expect_false("sites" %in% names(vp))

  ap <- x$all_positions
  expect_equal(ap$population, c("popA", "popB"))
  expect_equal(ap$sites, c(100000, 100000))
  expect_equal(ap$variant_sites, c(15, 15))
  expect_equal(ap$pct_polymorphic_loci, c(40.0, 33.3))
  expect_true(all(ap$obs_het < vp$obs_het))  # per-sequenced-site << per-variant-record
})

test_that("read_sumstats_summary() reports found-vs-expected on a header it can't read positionally", {
  bad <- c(
    "# Variant positions",
    "# Pop ID\tPrivate\tNum_Indv\tVar\tStdErr",  # truncated header, missing statistics
    "popA\t2\t4\t0\t0",
    "# All positions (variant and fixed)",
    "# Pop ID\tPrivate\tSites\tVariant_Sites\tPolymorphic_Sites\t%Polymorphic_Loci\tNum_Indv\tVar\tStdErr\tP\tVar\tStdErr\tObs_Het\tVar\tStdErr\tObs_Hom\tVar\tStdErr\tExp_Het\tVar\tStdErr\tExp_Hom\tVar\tStdErr\tPi\tVar\tStdErr\tFis\tVar\tStdErr",
    "popA\t2\t100000\t15\t6\t40\t4\t0\t0\t0.9\t0\t0\t0.2\t0\t0\t0.8\t0\t0\t0.19\t0\t0\t0.81\t0\t0\t0.195\t0\t0\t-0.05\t0\t0"
  )
  f <- tempfile(fileext = ".tsv")
  writeLines(bad, f)
  expect_error(read_sumstats_summary(f), "does not.*match the expected")
})

test_that("read_sumstats_summary() requires both blocks, in order", {
  one_block <- c(
    "# All positions (variant and fixed)",
    "# Pop ID\tPrivate\tSites\tVariant_Sites\tPolymorphic_Sites\t%Polymorphic_Loci\tNum_Indv\tVar\tStdErr\tP\tVar\tStdErr\tObs_Het\tVar\tStdErr\tObs_Hom\tVar\tStdErr\tExp_Het\tVar\tStdErr\tExp_Hom\tVar\tStdErr\tPi\tVar\tStdErr\tFis\tVar\tStdErr",
    "popA\t2\t100000\t15\t6\t40\t4\t0\t0\t0.9\t0\t0\t0.2\t0\t0\t0.8\t0\t0\t0.19\t0\t0\t0.81\t0\t0\t0.195\t0\t0\t-0.05\t0\t0"
  )
  f <- tempfile(fileext = ".tsv")
  writeLines(one_block, f)
  expect_error(read_sumstats_summary(f), "expected exactly 2 block title lines")
})

test_that("read_sumstats_summary() errors on a nonexistent file", {
  expect_error(read_sumstats_summary("no/such/file.tsv"), "File not found")
})
