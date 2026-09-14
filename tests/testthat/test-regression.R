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
  ## The full report, printed as the command-line script prints it (wide, so
  ## no table is split across lines).
  old <- options(width = 200)
  on.exit(options(old), add = TRUE)
  list(res = res, out = capture.output(print(summary(res))), dir = od)
}

## Same rule as the old check_golden.R: identical columns and rows, numbers
## within 1e-6 (the files are already rounded to 4-5 decimals), NA where NA.
## p_wilcox once needed a looser tolerance: R 4.6.0 changed wilcox.test()'s
## default for tied values from the normal approximation to exact
## conditional inference, which moved these p-values by ~0.005 on
## win-builder. .two_sample() now passes `exact` explicitly (the pre-4.6
## rule), so the golden values should hold on every R version.
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
    expect_equal(a[[col]], g[[col]], tolerance = 1e-6, info = paste0(golden_name, ": ", col))
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
    run <- function() invisible(capture.output(suppressMessages(
      het_between_pops(lg(cs[["vcf"]]), lg("popmap.tsv"), min_call = 0.9, outdir = od))))
    ## one_dead_ind has an individual with no genotype: it is kept (NA in the
    ## individual table, left out of every test) and named in a warning.
    if (cs[["stem"]] == "one_dead_ind") expect_warning(run(), "A04") else run()
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
    args <- utils::modifyList(list(vcf = lg(b[[1]]), popmap = lg(b[[2]]), g = 20,
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
    args <- utils::modifyList(list(vcf = lg(o[[1]]), popmap = lg(o[[2]]), g = 20,
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
  ## One individual with no genotype anywhere is kept, and named in a warning.
  expect_warning(suppressMessages(invisible(capture.output(
    res <- het_between_pops(lg("one_dead_ind.vcf.gz"), lg("popmap.tsv"), min_call = 0.9,
                            outdir = od)))), "look like failed libraries")
  expect_identical(res$settings$flagged_individuals, "A04")
  expect_true(is.na(res$individual_heterozygosity$heterozygosity[
    res$individual_heterozygosity$sample == "A04"]))
})

test_that("the Wilcoxon p-value uses the same exact/approximate rule on every R version", {
  ## No ties, small groups: the exact p-value.
  a <- c(0.31, 0.29, 0.35, 0.33, 0.30)
  b <- c(0.25, 0.27, 0.24, 0.28)
  row <- RADdiversity:::.two_sample(a, b, "A", "B", 100, "heterozygosity", verbose = FALSE)
  expect_equal(row$p_wilcox, stats::wilcox.test(a, b, exact = TRUE)$p.value)
  ## Ties: the normal approximation with continuity correction (R >= 4.6.0
  ## would otherwise switch to exact conditional inference).
  b_tied <- c(0.25, 0.29, 0.24, 0.28)
  row <- RADdiversity:::.two_sample(a, b_tied, "A", "B", 100, "heterozygosity", verbose = FALSE)
  expect_equal(row$p_wilcox,
               suppressWarnings(stats::wilcox.test(a, b_tied, exact = FALSE, correct = TRUE))$p.value)
})

test_that("het_between_pops() gives an overall test for 3 or more populations", {
  res <- suppressMessages(het_between_pops(lg("sim.allsnps.vcf.gz"), lg("pm_three.tsv"),
                                           min_call = 0.9, nboot_g2 = 0, verbose = FALSE))
  om <- res$omnibus
  expect_equal(om$statistic, c("heterozygosity", "F"))
  expect_true(all(is.finite(om$p_welch)) && all(is.finite(om$p_kruskal)))

  ## The same numbers as oneway.test() and kruskal.test() on per-individual
  ## heterozygosity over the loci every population clears.
  H <- read_stacks_vcf(lg("sim.allsnps.vcf.gz"), verbose = FALSE)
  pops <- read_popmap(lg("pm_three.tsv"), H$samples, verbose = FALSE)
  sets <- RADdiversity:::.population_locus_sets(H, pops, 0.9)
  loci <- rowSums(sets) == ncol(sets)
  het <- unlist(lapply(pops, function(ids)
    colMeans(H$A1[loci, ids, drop = FALSE] != H$A2[loci, ids, drop = FALSE], na.rm = TRUE)))
  group <- factor(rep(names(pops), lengths(pops)), levels = names(pops))
  expect_equal(om$n_loci[1], sum(loci))
  expect_equal(om$p_welch[1], stats::oneway.test(het ~ group, var.equal = FALSE)$p.value)
  expect_equal(om$p_kruskal[1], stats::kruskal.test(het, group)$p.value)

  ## Too few loci shared by every population: NA rows, not an error.
  none <- RADdiversity:::.omnibus_tests(H, pops, sets, 0.9, sum(loci) + 1L, verbose = FALSE)
  expect_true(all(is.na(none$p_welch)))
  ## Written to its own file with outdir.
  od <- tempfile("legacy-")
  suppressMessages(het_between_pops(lg("sim.allsnps.vcf.gz"), lg("pm_three.tsv"), min_call = 0.9,
                                    nboot_g2 = 0, outdir = od, verbose = FALSE))
  expect_true(file.exists(file.path(od, "het_between_pops_omnibus.allsnps.tsv")))
})

test_that("the full report never names a column it does not print", {
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
