## Ported from the pre-package test/run_tests.sh harness: golden numeric
## values for diversity_stats() and het_between_pops(), every edge case the
## functions must reject or accept, and the guard that the printed report
## never names a column it does not print. (The command-line cases -- no
## arguments, a missing --g, unknown or bare flags -- are in test-cli.R.)
##
## The fixtures in fixtures/legacy/ were written by
## fixtures/legacy/make_fixtures.R (set.seed(4242): 2 populations of 15 and
## 10 individuals with FIS = 0.10, 400 RAD tags of 1-4 SNPs each, plus
## deliberately broken files) and then gzipped. The golden TSVs are the
## output of the version that the original harness validated.

lg <- function(name) test_path("fixtures", "legacy", name)

run_div <- function(vcf, popmap = "popmap.tsv", g = 20, ...) {
  od <- tempfile("legacy-"); dir.create(od)
  res <- NULL
  suppressMessages(
    res <- diversity_stats(lg(vcf), lg(popmap), g = g, nboot = 0, outdir = od, ...))
  ## The report is what print() on the result shows.
  list(res = res, out = capture.output(print(res)), dir = od)
}

## Same rule as the old check_golden.R: identical columns and rows, numbers
## within 1e-6 (the files are already rounded to 4-5 decimals), NA where NA --
## except p_wilcox/p_wilcox_BH, which get a looser tolerance (below). Those
## come straight out of stats::wilcox.test() with its default exact/normal-
## approximation selection; when the pooled sample has ties, R's own tie
## handling there has shifted slightly between R versions (observed: R 4.5.3
## vs 4.6.1, same input `a`/`b` -- mean1/mean2/diff/p_welch/hedges_g all
## matched exactly, only p_wilcox moved by ~0.004-0.006). That is an R
## implementation detail this package does not control, not a computation
## bug, so asserting it to 1e-6 makes the suite fail on R updates rather than
## on package regressions.
wide_tol_cols <- c("p_wilcox", "p_wilcox_BH")
expect_golden <- function(actual_file, golden_name) {
  a <- utils::read.delim(actual_file, check.names = FALSE, stringsAsFactors = FALSE)
  g <- utils::read.delim(lg(file.path("golden", golden_name)), check.names = FALSE,
                         stringsAsFactors = FALSE)
  ## Columns added since the golden files were made are allowed (after the
  ## golden ones); every golden column must still be there, in order.
  added_since <- c("F")
  expect_identical(intersect(names(a), names(g)), names(g), info = golden_name)
  expect_true(all(setdiff(names(a), names(g)) %in% added_since), info = golden_name)
  expect_identical(nrow(a), nrow(g), info = golden_name)
  for (col in names(g)) {
    tol <- if (col %in% wide_tol_cols) 0.02 else 1e-6
    expect_equal(a[[col]], g[[col]], tolerance = tol, info = paste0(golden_name, ": ", col))
  }
}

test_that("diversity_stats() reproduces the golden values", {
  for (stem in c("allsnps", "haps", "miss10")) {
    vcf <- if (stem == "miss10") "miss10.vcf.gz" else paste0("sim.", stem, ".vcf.gz")
    r <- run_div(vcf)
    for (kind in c("per_population", "richness")) {
      f <- sprintf("diversity_%s.%s.tsv", kind, stem)
      expect_golden(file.path(r$dir, f), f)
    }
    unlink(r$dir, recursive = TRUE)
  }
})

test_that("het_between_pops() reproduces the golden values", {
  cases <- list(c(vcf = "sim.allsnps.vcf.gz",  stem = "allsnps",      golden = "hbp_allsnps"),
                c(vcf = "one_dead_ind.vcf.gz", stem = "one_dead_ind", golden = "hbp_onedead"))
  for (cs in cases) {
    od <- tempfile("legacy-"); dir.create(od)
    invisible(capture.output(suppressMessages(
      het_between_pops(lg(cs[["vcf"]]), lg("popmap.tsv"), min_call = 0.9, outdir = od))))
    for (kind in c("individual_heterozygosity", "het_between_pops_tests"))
      expect_golden(file.path(od, sprintf("%s.%s.tsv", kind, cs[["stem"]])),
                    sprintf("%s.%s.tsv", cs[["golden"]], kind))
    unlink(od, recursive = TRUE)
  }
})

test_that("diversity_stats() rejects bad input with a message naming the problem", {
  bad <- list(
    list("sim.allsnps.vcf.gz",       "popmap.tsv",     list(g = 999),  "exceeds the smallest population"),
    list("sim.allsnps.vcf.gz",       "popmap.tsv",     list(g = 1),    "at least 2 gene copies"),
    list("nope.vcf",                 "popmap.tsv",     list(),         "VCF file not found"),
    list("sim.allsnps.vcf.gz",       "nope.tsv",       list(),         "Popmap file not found"),
    list("sim.allsnps.vcf.gz",       "pm_dup.tsv",     list(),         "more than once"),
    list("sim.allsnps.vcf.gz",       "pm_nomatch.tsv", list(g = 2),    "No popmap sample names match"),
    list("sim.allsnps.vcf.gz",       "pm_onepop.tsv",  list(g = 2),    "at least 2 populations"),
    list("sim.allsnps.vcf.gz",       "pm_n1.tsv",      list(g = 2),    "fewer than 2 individuals"),
    list("bad_nohdr.vcf.gz",         "popmap.tsv",     list(),         "#CHROM"),
    list("bad_norec.vcf.gz",         "popmap.tsv",     list(),         "No variant records"),
    list("bad_allmiss.vcf.gz",       "popmap.tsv",     list(),         "No record has >= min_n"),
    list("bad_nocomplete.vcf.gz",    "popmap.tsv",     list(complete_case = TRUE),
         "No record is genotyped in every individual"),
    list("bad_mono.vcf.gz",          "popmap.tsv",     list(),         "No locus is polymorphic"),
    list("locus_driven_miss.vcf.gz", "popmap.tsv",     list(min_n = 1), "min_n must be"))
  for (b in bad) {
    args <- utils::modifyList(list(vcf_file = lg(b[[1]]), popmap_f = lg(b[[2]]), g = 20,
                                   nboot = 0, outdir = tempfile("legacy-")), b[[3]])
    expect_error(suppressMessages(capture.output(do.call(diversity_stats, args))),
                 b[[4]], fixed = TRUE, info = paste(b[[1]], b[[2]]))
  }
})

test_that("diversity_stats() accepts the edge cases it should", {
  good <- list(
    list("sim.allsnps.vcf.gz",       "popmap.tsv",        list(sites = 10)),  # implausible sites: ignored
    list("sim.allsnps.vcf.gz",       "pm_subset.tsv",     list(g = 16)),
    list("sim.allsnps.vcf.gz",       "pm_numeric.tsv",    list()),
    list("sim.allsnps.vcf.gz",       "pm_underscore.tsv", list()),
    list("sim.allsnps.vcf.gz",       "pm_three.tsv",      list(g = 16)),
    list("bad_nocomplete.vcf.gz",    "popmap.tsv",        list()),            # available data survives
    list("miss10.vcf.gz",            "popmap.tsv",        list()),
    list("miss10.vcf.gz",            "popmap.tsv",        list(complete_case = TRUE)),
    list("locus_driven_miss.vcf.gz", "popmap.tsv",        list()),
    list("locus_driven_miss.vcf.gz", "popmap.tsv",        list(complete_case = TRUE)),
    list("locus_driven_miss.vcf.gz", "popmap.tsv",        list(min_n = 3)))
  for (o in good) {
    args <- utils::modifyList(list(vcf_file = lg(o[[1]]), popmap_f = lg(o[[2]]), g = 20,
                                   nboot = 0, outdir = tempfile("legacy-")), o[[3]])
    expect_no_error(suppressMessages(capture.output(do.call(diversity_stats, args))),
                    message = paste(o[[1]], o[[2]]))
  }
})

test_that("per-site values are refused on a haplotype VCF and given on a SNP VCF", {
  expect_null(run_div("sim.haps.vcf.gz", sites = 1e6)$res$autosomal)
  expect_false(is.null(run_div("sim.allsnps.vcf.gz", sites = 1e6)$res$autosomal))
  ## One SNP per tag is still a SNP file (single-base alleles, biallelic).
  msgs <- capture_messages(invisible(capture.output(
    diversity_stats(lg("sim.onesnp.vcf.gz"), lg("popmap.tsv"), g = 20, nboot = 0,
                    sites = 1e6, outdir = tempfile("legacy-")))))
  expect_true(any(grepl("record type: SNP", msgs)))
})

test_that("het_between_pops() edge cases", {
  od <- tempfile("legacy-")
  expect_error(het_between_pops(lg("sim.allsnps.vcf.gz"), lg("popmap.tsv"), min_call = 5,
                                outdir = od), "min_call must be in 0-1")
  for (cs in list(c("sim.haps.vcf.gz", "popmap.tsv"), c("sim.allsnps.vcf.gz", "pm_three.tsv")))
    expect_no_error(suppressMessages(capture.output(
      het_between_pops(lg(cs[1]), lg(cs[2]), min_call = 0.9, outdir = od))))
  ## One individual with no genotype anywhere is excluded, loudly.
  msgs <- capture_messages(invisible(capture.output(
    het_between_pops(lg("one_dead_ind.vcf.gz"), lg("popmap.tsv"), min_call = 0.9, outdir = od))))
  expect_true(any(grepl("EXCLUDED", msgs)))
})

test_that("the printed report never names a column it does not print", {
  ## Ported from test/check_labels.R: every Ho_x/He_x/Fis_x/Ar_x/priv_x/pct_x
  ## identifier mentioned in the report's prose must be a printed column
  ## header -- except mentions that explicitly point to the OTHER run.
  id_re <- "\\b(?:Ho|He|Fis|Ar|priv|pct)[A-Za-z0-9_]*\\b"
  for (f in c("sim.allsnps.vcf.gz", "sim.haps.vcf.gz")) {
    out <- run_div(f, sites = 1e6)$out
    txt <- paste(out, collapse = "\n")
    cols <- unique(unlist(strsplit(trimws(grep("^ *population ", out, value = TRUE)), " +")))
    expect_gt(length(cols), 5)                       # a vacuous pass would find none
    ment <- unique(unlist(regmatches(txt, gregexpr(id_re, txt, perl = TRUE))))
    ment <- ment[grepl("_", ment)]
    cross_txt <- paste(grep("SKIPPED|Re-run|Run populations|TAKE FROM|DO NOT take",
                            out, value = TRUE), collapse = "\n")
    cross <- unique(unlist(regmatches(cross_txt, gregexpr(id_re, cross_txt, perl = TRUE))))
    expect_equal(setdiff(setdiff(ment, cross), cols), character(0), info = f)
    expect_false("pi" %in% cols)                     # nothing is ever called `pi`
    expect_false(any(grepl("TAKE FROM THIS RUN.*\\bpi\\b|Ho, He, pi", out)))
  }
})
