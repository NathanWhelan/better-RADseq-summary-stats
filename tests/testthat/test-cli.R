fx <- function(name) test_path("fixtures", name)

## Run one command-line command in-process; returns its exit status plus
## whatever it printed and messaged.
run_cli <- function(cmd, ...) {
  status <- NULL
  msgs <- NULL
  out <- capture.output(msgs <- capture_messages(
    status <- RADdiversity:::.cli_main(cmd, c(...))))
  list(status = status, out = out, msgs = msgs)
}

test_that(".cli_parse() separates flags, switches and positionals, and names the fix for mistakes", {
  p <- RADdiversity:::.cli_parse(c("a.vcf", "--g=20", "--complete-case", "pm.tsv"),
                                 c("g", "nboot"), "complete-case")
  expect_equal(p$positional, c("a.vcf", "pm.tsv"))
  expect_equal(p$flags, list(g = "20"))
  expect_equal(p$switches, "complete-case")
  expect_error(RADdiversity:::.cli_parse("--bogus=1", "g"), "Unrecognized flag: --bogus")
  expect_error(RADdiversity:::.cli_parse("--g", "g"), "--g=VALUE")
  expect_error(RADdiversity:::.cli_parse("--nope", "g"), "Unrecognized flag: --nope")
})

test_that("diversity_stats CLI: --sites=N is a number, not a file path", {
  ## Regression: command-line values arrive as text and diversity_stats()
  ## reads text as a path, so --sites=100000 failed with "file not found".
  outdir <- tempfile("cli-"); dir.create(outdir)
  on.exit(unlink(outdir, recursive = TRUE), add = TRUE)
  res <- run_cli("diversity_stats", fx("small.snps.vcf"), fx("small_popmap.tsv"),
                 "--g=4", "--nboot=0", "--sites=100000", paste0("--outdir=", outdir))
  expect_equal(res$status, 0L)
  expect_true(file.exists(file.path(outdir, "diversity_autosomal.snps.tsv")))
})

test_that("diversity_stats CLI: --sites= also accepts a sumstats_summary path", {
  outdir <- tempfile("cli-"); dir.create(outdir)
  on.exit(unlink(outdir, recursive = TRUE), add = TRUE)
  res <- run_cli("diversity_stats", fx("small.snps.vcf"), fx("small_popmap.tsv"),
                 "--g=4", "--nboot=0", paste0("--sites=", fx("sumstats_summary_small.tsv")),
                 paste0("--outdir=", outdir))
  expect_equal(res$status, 0L)
  expect_true(file.exists(file.path(outdir, "diversity_autosomal.snps.tsv")))
})

test_that("diversity_stats CLI: bad invocations exit 1 with a usable message", {
  expect_equal(run_cli("diversity_stats")$status, 1L)                     # no args
  no_g <- run_cli("diversity_stats", fx("small.snps.vcf"), fx("small_popmap.tsv"))
  expect_equal(no_g$status, 1L)
  expect_true(any(grepl("--g=N", no_g$msgs)))
  expect_equal(run_cli("diversity_stats", fx("small.snps.vcf"), fx("small_popmap.tsv"),
                       "--g=4", "extra")$status, 1L)                      # leftover positional
  expect_equal(run_cli("diversity_stats", fx("small.snps.vcf"), fx("small_popmap.tsv"),
                       "--g", "4")$status, 1L)                            # bare --g
})

test_that("het_between_pops CLI: flag form and the older positional form agree", {
  d1 <- tempfile("cli-"); d2 <- tempfile("cli-"); dir.create(d1); dir.create(d2)
  on.exit(unlink(c(d1, d2), recursive = TRUE), add = TRUE)
  a <- run_cli("het_between_pops", fx("small.haps.vcf"), fx("small_popmap.tsv"),
               "--min-call=0.5", paste0("--outdir=", d1))
  b <- run_cli("het_between_pops", fx("small.haps.vcf"), fx("small_popmap.tsv"), "0.5", d2)
  expect_equal(c(a$status, b$status), c(0L, 0L))
  expect_identical(readLines(file.path(d1, "het_between_pops_tests.haps.tsv")),
                   readLines(file.path(d2, "het_between_pops_tests.haps.tsv")))
  expect_equal(run_cli("het_between_pops", fx("small.haps.vcf"), fx("small_popmap.tsv"),
                       "--min-call=5")$status, 1L)                       # out of range
})

test_that("diversity_core CLI needs --selftest", {
  expect_equal(run_cli("diversity_core")$status, 1L)
})
