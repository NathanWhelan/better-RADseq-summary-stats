## summary() prints a SHORT view by default and the full report with
## `details = TRUE`. These tests check what the short views promise, for all
## four functions (diversity_stats(), differentiation_stats(),
## het_between_pops(), pi_allsites()) and for the helpers they share:
##   - the short view is short and no line is wider than 80 characters;
##   - every column of "estimate (SE)" cells has "(SE)" in its header, and a
##     cell is exactly the estimate and the ONE recommended SE, formatted with
##     the package's number rule; columns with no SE say so;
##   - none of the other uncertainty columns is shown by default;
##   - populations with long names never lose their row label.

lg <- function(name) test_path("fixtures", "legacy", name)

## Print at the default console width, and return the lines.
lines_of <- function(x, ...) {
  old <- options(width = 80)
  on.exit(options(old), add = TRUE)
  capture.output(print(summary(x, ...)))
}

## The 15- and 10-individual simulated data, with popmap names of the length
## real studies use.
long_popmap <- function() {
  pm <- readLines(lg("popmap.tsv"))
  pm <- sub("\tpopA$", "\tSipsey_Fork_2019", pm)
  pm <- sub("\tpopB$", "\tCahaba_River_2020", pm)
  path <- tempfile(fileext = ".tsv")
  writeLines(pm, path)
  path
}

div_run <- function(vcf, popmap = lg("popmap.tsv"), ...)
  suppressMessages(diversity_stats(lg(vcf), popmap, g = 20, nboot = 100, seed = 1, ...))

## ---- the shared helpers -------------------------------------------------------

test_that(".est_se() gives the SE 2 significant digits and the estimate the same decimals", {
  est_se <- RADdiversity:::.est_se
  expect_identical(est_se(c(0.313421, 0.246091), c(0.00988, 0.0162)),
                   c("0.3134 (0.0099)", "0.2461 (0.0162)"))
  ## Decimals come from the smallest SE *after* rounding: 0.0000996 counts as 0.00010.
  expect_identical(est_se(c(0.005858, 0.004996), c(0.0000996, 0.000121)),
                   c("0.00586 (0.00010)", "0.00500 (0.00012)"))
  expect_identical(est_se(c(0.5, NA), c(NA, 0.1)), c("0.50 (NA)", "NA (0.10)"))
  expect_identical(est_se(0.3, NA), "0.300 (NA)")            # no usable SE: 3 decimals
})

test_that("an independent copy of the estimate (SE) rule agrees with .est_se()", {
  ## The rule in words: SE to 2 significant digits, the estimate to the same
  ## number of decimals, fixed within a column (see .est_se()).
  readme_est_se <- function(x, se) {
    d <- max(0, 1 - floor(log10(signif(min(se[se > 0]), 2))))
    sprintf("%.*f (%.*f)", d, x, d, se)
  }
  res <- div_run("sim.allsnps.vcf.gz", sites = 1e6)
  for (cols in list(c("Ho", "Ho_se_combined"), c("He", "He_se_combined"), c("Fis", "Fis_se_combined")))
    expect_identical(RADdiversity:::.est_se(res$per_population[[cols[1]]], res$per_population[[cols[2]]]),
                     readme_est_se(res$per_population[[cols[1]]], res$per_population[[cols[2]]]))
  a <- res$autosomal
  expect_identical(RADdiversity:::.est_se(a$He_autosomal, a$He_autosomal_se_combined),
                   readme_est_se(a$He_autosomal, a$He_autosomal_se_combined))
})

test_that(".stat_blocks() keeps Ar apart from privAr and puts the label first in every block", {
  df <- data.frame(population = c("a", "b"), n = 1:2, Ar = 1:2, Ar_se = 3:4, Ar_n = 5:6,
                   privAr = 7:8, privAr_se = 9:10, priv_total = 11:12, priv_total_se = 13:14,
                   Ho = 1:2, Ho_autosomal = 3:4)
  b <- RADdiversity:::.stat_blocks(df, c("Ar", "privAr", "priv_total", "Ho", "Ho_autosomal"),
                                   always = "n")
  expect_identical(names(b$Ar), c("population", "n", "Ar", "Ar_se", "Ar_n"))
  expect_identical(names(b$privAr), c("population", "n", "privAr", "privAr_se"))
  expect_identical(names(b$priv_total), c("population", "n", "priv_total", "priv_total_se"))
  expect_identical(names(b$Ho), c("population", "n", "Ho"))         # not Ho_autosomal
  expect_identical(names(b$Ho_autosomal), c("population", "n", "Ho_autosomal"))
  ## with only "Ho" asked for, Ho_autosomal has to go somewhere and goes to Ho
  only_ho <- RADdiversity:::.stat_blocks(df, "Ho", always = "n")
  expect_identical(names(only_ho$Ho), c("population", "n", "Ho", "Ho_autosomal"))
  expect_true(all(vapply(b, function(d) names(d)[1] == "population", logical(1))))
})

test_that(".print_compact() never prints nothing, and .stat_blocks() drops a statistic with no columns", {
  label_only <- data.frame(population = c("a", "b"))
  expect_gt(length(capture.output(RADdiversity:::.print_compact(label_only))), 0L)
  b <- RADdiversity:::.stat_blocks(data.frame(population = "a", Ho = 1, Ho_se = 2), c("Ho", "Fis"),
                                   always = "n")
  expect_identical(names(b), "Ho")                          # no block holding only the label
})

test_that(".print_compact() measures numbers as print() shows them, so no block is left without its label", {
  ## 9.25e-05 prints as 0.0000925: a width taken from as.character() is too small
  ## and lets print() wrap a block by itself, with no row label.
  set.seed(1)
  df <- data.frame(population = c("north", "south"),
                   a_autosomal = c(0.005297, 0.00416), a_se = c(0.0000925, 0.0001053),
                   a_se_ind = c(0.0001397, 0.000253), a_se_combined = c(0.0001676, 0.0002741),
                   a_lo = c(0.00509, 0.003992), a_hi = c(0.005423, 0.004261),
                   b_autosomal = c(0.005858, 0.004996), b_se = c(8.424e-05, 1.08e-04),
                   b_se_ind = c(5.316e-05, 5.462e-05), b_se_combined = c(9.961e-05, 1.21e-04),
                   b_lo = c(0.005695, 0.004763), b_hi = c(0.006005, 0.005187))
  for (width in c(80, 60, 45, 30)) {
    out <- capture.output(RADdiversity:::.print_compact(df, width = width))
    expect_lt(max(nchar(out)), max(width, nchar(out[1]) + 1L), label = paste("width", width))
    headers <- sum(grepl("^ *population", out))
    expect_equal(sum(grepl("^ *north ", out)), headers, info = paste("width", width))
    expect_equal(sum(grepl("^ *south ", out)), headers, info = paste("width", width))
    expect_equal(length(out), 3L * headers, info = paste("width", width))    # header + 2 rows each
  }
})

test_that(".print_compact() repeats the row label on every block when it has to wrap", {
  df <- data.frame(population = c("Sipsey_Fork_2019", "Cahaba_River_2020"), n = c(15L, 10L),
                   `Ho (SE)` = c("0.3134 (0.0099)", "0.2461 (0.0162)"),
                   `He (SE)` = c("0.3466 (0.0059)", "0.2956 (0.0072)"), check.names = FALSE)
  wide <- capture.output(RADdiversity:::.print_compact(df, width = 80))
  expect_length(wide, 3L)
  narrow <- capture.output(RADdiversity:::.print_compact(df, width = 40))
  expect_gt(length(narrow), 3L)
  ## every block starts with its header row and then both labels
  expect_equal(sum(grepl("Sipsey_Fork_2019", narrow)), sum(grepl("population", narrow)))
  expect_equal(sum(grepl("Cahaba_River_2020", narrow)), sum(grepl("population", narrow)))
})

## ---- diversity_stats() ------------------------------------------------------

test_that("the short diversity summary is short, fits 80 columns and labels every estimate (SE)", {
  snps <- div_run("sim.allsnps.vcf.gz", sites = 1e6)
  haps <- div_run("sim.haps.vcf.gz")
  for (res in list(snps, haps)) {
    out <- lines_of(res)
    expect_lt(length(out), 45)
    expect_lte(max(nchar(out)), 80)
    expect_true(any(grepl("each cell is estimate (standard error)", out, fixed = TRUE)))
    expect_true(any(grepl("TAKE FROM THIS RUN", out)))
    expect_true(any(grepl("^CHECKS", out)))
    ## the other uncertainty columns and the support counts are not shown
    expect_false(any(grepl("_lo|_hi|_se_ind|Ar_n|priv_total", out)))
  }

  ## SNP VCF: Ho and He (and per-site values), each with the combined SE.
  ## pct_poly is not in the primary table.
  out <- lines_of(snps)
  header <- grep("^ *population +n ", out, value = TRUE)
  expect_match(header, "Ho \\(SE\\) +He \\(SE\\)$")
  expect_false(grepl("Fis|Ar|private|pct_poly", header))
  expect_false(any(grepl("pct_poly", out)))
  pp <- snps$per_population
  est_se <- RADdiversity:::.est_se
  expect_true(all(vapply(est_se(pp$Ho, pp$Ho_se_combined), function(cell) any(grepl(cell, out, fixed = TRUE)),
                         logical(1))))
  expect_true(all(vapply(est_se(pp$He, pp$He_se_combined), function(cell) any(grepl(cell, out, fixed = TRUE)),
                         logical(1))))
  expect_false(any(grepl(est_se(pp$Ho, pp$Ho_se)[1], out, fixed = TRUE)))    # not the locus SE
  a <- snps$autosomal
  site_header <- grep("Ho_autosomal \\(SE\\)", out, value = TRUE)
  expect_length(site_header, 1L)
  expect_true(any(grepl(est_se(a$He_autosomal, a$He_autosomal_se_combined)[1], out, fixed = TRUE)))

  ## Haplotype VCF: Fis with the combined SE, Ar and the rarefied private
  ## alleles (the count priv_total, not privAr) with the locus SE.
  out <- lines_of(haps)
  header <- grep("^ *population +n ", out, value = TRUE)
  expect_match(header, "Fis \\(SE\\) +Ar \\(SE\\) +rarefied private alleles \\(SE\\)")
  expect_false(grepl("Ho|He|pct_poly", header))
  pp <- haps$per_population
  rich <- haps$richness
  expect_true(any(grepl(est_se(pp$Fis, pp$Fis_se_combined)[1], out, fixed = TRUE)))
  expect_true(any(grepl(est_se(rich$Ar, rich$Ar_se)[1], out, fixed = TRUE)))
  expect_true(any(grepl(est_se(rich$priv_total, rich$priv_total_se)[1], out, fixed = TRUE)))
  expect_false(any(grepl(est_se(rich$privAr, rich$privAr_se)[1], out, fixed = TRUE)))
})

test_that("the short diversity summary says which SE each column uses, and what simulations showed", {
  snps <- lines_of(div_run("sim.allsnps.vcf.gz"))
  expect_true(any(grepl("SE counts variation among RAD loci and among individuals", snps)))
  expect_true(any(grepl(RADdiversity:::.coverage$combined, snps, fixed = TRUE)))
  expect_false(any(grepl("pct_poly", snps)))

  haps <- lines_of(div_run("sim.haps.vcf.gz"))
  expect_true(any(grepl("SEs of Ar and rarefied private alleles count variation among RAD loci only", haps)))
  expect_true(any(grepl(RADdiversity:::.coverage$loci, haps, fixed = TRUE)))
  expect_true(any(grepl("allelic richness", haps)))                # what Ar is
  expect_match(gsub("\\s+", " ", paste(haps, collapse = " ")),
               "found only in this population")                     # what private alleles are
  expect_match(gsub("\\s+", " ", paste(haps, collapse = " ")),
               "It is an expected value")                           # not a whole number
  expect_match(gsub("\\s+", " ", paste(haps, collapse = " ")),
               "added up over the [0-9,]+ loci")                      # and over which loci

  ## Without the individual SEs the cells hold the locus SE, and the header says so.
  off <- div_run("sim.allsnps.vcf.gz", se_individuals = FALSE)
  out <- lines_of(off)
  expect_true(any(grepl("se_individuals = FALSE was set", out)))
  expect_true(any(grepl("SE counts variation among RAD loci only", out)))
  expect_false(any(grepl(RADdiversity:::.coverage$combined, out, fixed = TRUE)))
  expect_false(any(grepl("_se_ind|_se_combined", out)))
  pp <- off$per_population
  expect_true(any(grepl(RADdiversity:::.est_se(pp$Ho, pp$Ho_se)[1], out, fixed = TRUE)))
})

test_that("the short diversity summary lists every check as ok, info or look", {
  out <- lines_of(div_run("sim.haps.vcf.gz"))
  checks <- grep("^  (ok|info|look) ", out, value = TRUE)
  expect_gte(length(checks), 3L)
  expect_true(any(grepl("^  ok +every population had enough genotypes", checks)))
  ## a fired check reads "look"
  set.seed(3)
  res <- div_run("sim.haps.vcf.gz", complete_case = FALSE)
  res$settings$boot <- "individuals"                             # a comparison mode
  expect_true(any(grepl("^  look +boot = ", lines_of(res))))
})

test_that("the short diversity summary says why a population has no per-site value", {
  ## popB's `sites` is smaller than its variant records, so it is ignored for popB only
  res <- div_run("sim.allsnps.vcf.gz", sites = c(popA = 1e6, popB = 5))
  expect_true(is.na(res$autosomal$He_autosomal[res$autosomal$population == "popB"]))
  joined <- gsub("\\s+", " ", paste(lines_of(res), collapse = " "))
  expect_match(joined, "look no per-site value for popB \\(its `sites` value, 5, is smaller than its [0-9,]+ variant records\\)")
  ## the same run with a good value for both has no such line
  ok <- gsub("\\s+", " ", paste(lines_of(div_run("sim.allsnps.vcf.gz", sites = 1e6)), collapse = " "))
  expect_false(grepl("no per-site value", ok))
})

test_that("print() keeps the population name on every block of a wide table", {
  res <- div_run("sim.allsnps.vcf.gz", sites = 1e6)
  old <- options(width = 80); on.exit(options(old), add = TRUE)
  out <- capture.output(print(res))
  ## A SNP run does not print Ar and private alleles; the haplotype run does.
  expect_false(any(grepl("^\\$richness", out)))
  hap_out <- capture.output(print(div_run("sim.haps.vcf.gz")))
  snp_out <- out
  for (tab in c("per_population", "richness", "autosomal")) {
    out <- if (tab == "richness") hap_out else snp_out
    start <- grep(paste0("^\\$", tab, "$"), out)
    stop_at <- min(c(grep("^\\$", out)[grep("^\\$", out) > start], grep("^Take from|^NOTE|^summary", out)))
    block <- out[(start + 1):(stop_at - 1)]
    block <- block[nzchar(trimws(block))]
    headers <- sum(grepl("^ *population ", block))
    expect_gt(headers, 0L, label = tab)
    ## a header row, then one row for each of the 2 populations, for every block
    expect_equal(sum(grepl("^ *popA ", block)), headers, info = tab)
    expect_equal(sum(grepl("^ *popB ", block)), headers, info = tab)
    expect_equal(length(block), 3L * headers, info = tab)
  }
})

test_that("the short diversity summary keeps every row label with long popmap names", {
  pm <- long_popmap()
  on.exit(unlink(pm), add = TRUE)
  for (vcf in c("sim.allsnps.vcf.gz", "sim.haps.vcf.gz")) {
    out <- lines_of(div_run(vcf, popmap = pm, sites = if (vcf == "sim.allsnps.vcf.gz") 1e6))
    expect_lte(max(nchar(out)), 80)
    rows <- grep("^ *(Sipsey_Fork_2019|Cahaba_River_2020) +(15|10) ", out, value = TRUE)
    expect_length(rows, 2L)                                   # both populations in the results table
  }
})

test_that("summary(details = TRUE) of diversity_stats() is the full report, with labelled blocks", {
  res <- div_run("sim.allsnps.vcf.gz", sites = 1e6)
  out <- lines_of(res, details = TRUE)
  brief <- lines_of(res)
  expect_gt(length(out), 2 * length(brief))
  ## the recommended-SE table and the advice come first, then the full tables
  expect_lt(grep("^RESULTS", out)[1], grep("^FULL TABLES", out)[1])
  expect_lt(grep("TAKE FROM THIS RUN", out)[1], grep("^FULL TABLES", out)[1])
  expect_gt(grep("^FULL TABLES", out)[1], 0)
  ## every uncertainty column is there, in one small table per statistic that
  ## starts with the population
  for (col in c("Ho_se", "Ho_se_ind", "Ho_se_combined", "Ho_lo", "Ho_hi", "He_se_combined",
                "Fis_se_combined", "Ar_se", "priv_total_lo", "Ho_autosomal_se_combined",
                "He_autosomal_hi"))
    expect_true(any(grepl(paste0("\\b", col, "\\b"), out)), info = col)
  heads <- grep("^ *population ", out, value = TRUE)
  expect_gte(length(heads), 8L)
  ## no block is split: every header row lists its whole statistic
  expect_true(any(grepl("population +n +Ho +Ho_se +Ho_se_ind +Ho_se_combined +Ho_lo +Ho_hi", heads)))
  ## what to report, with how often that SE held the truth in simulations
  expect_true(any(grepl("Report Ho_se_combined", out)))
  expect_true(any(grepl("Report Ar_se", out)))
})

## ---- differentiation_stats() ------------------------------------------------

dif_run <- function(vcf, popmap = lg("popmap.tsv"), ...)
  suppressMessages(differentiation_stats(lg(vcf), popmap, nboot = 100, seed = 1, ...))

test_that("the short differentiation summary shows FST and D as estimate (SE) and says beta has no SE", {
  skip_if_not_installed("hierfstat")
  res <- dif_run("sim.haps.vcf.gz")
  out <- lines_of(res)
  expect_lt(length(out), 45)
  expect_lte(max(nchar(out)), 80)
  expect_true(any(grepl("^DIFFERENTIATION:", out)))
  header <- grep("^ *pair ", out, value = TRUE)
  expect_match(header, "FST \\(SE\\) +D \\(SE\\) +beta")
  expect_false(grepl("FIS|_lo|_hi", header))
  pw <- res$pairwise
  est_se <- RADdiversity:::.est_se
  expect_true(any(grepl(est_se(pw$FST, pw$FST_se)[1], out, fixed = TRUE)))
  expect_true(any(grepl(est_se(pw$D, pw$D_se)[1], out, fixed = TRUE)))
  expect_true(any(grepl("beta +Weir & Goudet.*it has no SE", out)))
  expect_true(any(grepl("SE counts variation among RAD loci only \\(not checked by simulation\\)", out)))
  ## two populations: one sentence, no separate "global" and "pairwise" D
  expect_true(any(grepl("D uses the [0-9,]+ of [0-9,]+ records typed in both populations", out)))
  expect_false(any(grepl("global D|pairwise D uses", out)))
  expect_false(any(grepl("FST_lo|D_lo|FIS_se|_se_", out)))          # no other uncertainty columns
  ## what to report depends on the file
  expect_true(any(grepl("Report D \\(or beta\\) with FST", out)))
  snp <- lines_of(dif_run("sim.allsnps.vcf.gz"))
  expect_true(any(grepl("Report FST and D, each with its SE", snp)))
  expect_false(any(grepl("prefer D", snp)))
})

test_that("the short differentiation summary without beta says why, and has no beta column", {
  res <- dif_run("sim.allsnps.vcf.gz", beta = FALSE)
  out <- lines_of(res)
  expect_false(grepl("beta", grep("^ *pair ", out, value = TRUE)))
  expect_true(any(grepl("^  info +beta not computed \\(beta = FALSE\\)", out)))
  expect_false(any(grepl("beta +Weir & Goudet", out)))
})

test_that("with 3 or more populations the short summary adds an 'all populations' row", {
  res <- dif_run("sim.allsnps.vcf.gz", popmap = lg("pm_three.tsv"), beta = FALSE)
  out <- lines_of(res)
  expect_true(any(grepl("^ *all populations ", out)))
  expect_equal(sum(grepl("^ *g[1-3] - g[1-3] ", out)), nrow(res$pairwise))     # popmap: g1, g2, g3
  est_se <- RADdiversity:::.est_se
  g <- res$global
  expect_true(any(grepl(est_se(c(g$FST, res$pairwise$FST), c(g$FST_se, res$pairwise$FST_se))[1], out,
                        fixed = TRUE)))
})

test_that("the short differentiation summary keeps every row label with long popmap names", {
  pm <- long_popmap()
  on.exit(unlink(pm), add = TRUE)
  out <- lines_of(dif_run("sim.allsnps.vcf.gz", popmap = pm, beta = FALSE))
  expect_lte(max(nchar(out)), 80)
  expect_true(any(grepl("Sipsey_Fork_2019 - Cahaba_River_2020", out)))
})

test_that("the full differentiation report leaves out the beta table when beta was not computed", {
  res <- dif_run("sim.haps.vcf.gz", beta = FALSE)
  out <- lines_of(res, details = TRUE)
  expect_false(any(grepl("^ *pair +beta$", out)))                     # no empty beta table
  expect_false(any(grepl("beta is Weir & Goudet", out)))              # and no explanation of it
  expect_equal(sum(grepl("beta not computed \\(beta = FALSE\\)", out)), 1L)   # the reason, once
})

test_that("summary(details = TRUE) of differentiation_stats() has every column in labelled blocks", {
  res <- dif_run("sim.haps.vcf.gz", beta = FALSE)
  out <- lines_of(res, details = TRUE)
  expect_gt(length(out), length(lines_of(res)))
  for (col in c("FST_se", "FST_lo", "FST_hi", "FIS_se", "D_lo", "D_hi", "D_records"))
    expect_true(any(grepl(paste0("\\b", col, "\\b"), out)), info = col)
  expect_true(any(grepl("^ *scope +FST +FST_se +FST_lo +FST_hi", out)))
  expect_true(any(grepl("^ *pair +FST +FST_se +FST_lo +FST_hi", out)))
  expect_true(any(grepl("The FIS to report is", out)))
  expect_lt(grep("^RESULTS", out)[1], grep("^FULL TABLES", out)[1])
  expect_true(any(grepl("D uses the [0-9,]+ of [0-9,]+ records typed in both populations", out)))
})

## ---- pi_allsites() ----------------------------------------------------------

ex <- function(name) system.file("extdata", name, package = "RADdiversity")
pi_run <- function(popmap = ex("example_popmap.tsv"))
  suppressMessages(pi_allsites(ex("example.allsites.vcf.gz"), popmap, nboot = 50, seed = 1))
pi_cached <- local({ res <- NULL; function() { if (is.null(res)) res <<- pi_run(); res } })

test_that("the short pi summary shows pi_nc, pi, dxy and da_nc as estimate (SE)", {
  res <- pi_cached()
  out <- lines_of(res)
  expect_lt(length(out), 45)
  expect_lte(max(nchar(out)), 80)
  expect_true(any(grepl("^NUCLEOTIDE DIVERSITY:", out)))
  expect_match(grep("^ *population +sites ", out, value = TRUE), "pi_nc \\(SE\\) +pi \\(SE\\)")
  expect_match(grep("^ *pair ", out, value = TRUE), "dxy \\(SE\\) +da_nc \\(SE\\)")
  est_se <- RADdiversity:::.est_se
  p <- res$pi
  expect_true(all(vapply(est_se(p$pi_nc, p$pi_nc_se), function(cell) any(grepl(cell, out, fixed = TRUE)),
                         logical(1))))
  expect_true(all(vapply(est_se(p$pi, p$pi_se), function(cell) any(grepl(cell, out, fixed = TRUE)),
                         logical(1))))
  expect_true(any(grepl(est_se(res$dxy$dxy, res$dxy$dxy_se)[1], out, fixed = TRUE)))
  expect_true(any(grepl(est_se(res$dxy$da_nc, res$dxy$da_nc_se)[1], out, fixed = TRUE)))
  ## the other uncertainty columns, and da (built on pi), are not shown
  expect_false(any(grepl("_lo|_hi|pi_se|da \\(SE\\)", out)))
  expect_true(any(grepl("pi_nc +REPORT THIS", out)))
  expect_true(any(grepl("SE counts variation among RAD loci only \\(not checked by simulation\\)", out)))
  expect_true(any(grepl("PER-INDIVIDUAL HETEROZYGOSITY PER SITE", out)))
  expect_match(grep("^ *population +mean het", out, value = TRUE), "\\(no SE\\)")
  expect_true(any(grepl("^  info +each population's pi is averaged over its own sites", out)))
})

test_that("print() of a pi_allsites() result keeps the population name on every block", {
  old <- options(width = 80); on.exit(options(old), add = TRUE)
  pi_out <- capture.output(print(pi_cached()))
  start <- grep("^\\$pi$", pi_out)
  block <- pi_out[(start + 1):(grep("^\\$dxy", pi_out) - 1)]
  block <- block[nzchar(trimws(block))]
  headers <- sum(grepl("^ *population ", block))
  expect_gt(headers, 0L)
  expect_equal(sum(grepl("^ *north ", block)), headers)
  expect_equal(sum(grepl("^ *south ", block)), headers)
})

test_that("the short pi summary mentions complete_sites, and keeps row labels with long names", {
  res <- pi_cached()
  res$settings$complete_sites <- TRUE
  expect_true(any(grepl("^  info +complete_sites: each population uses only sites typed in all",
                        lines_of(res))))

  pm <- readLines(ex("example_popmap.tsv"))
  pm <- sub("\tnorth$", "\tSipsey_Fork_2019", pm)
  pm <- sub("\tsouth$", "\tCahaba_River_2020", pm)
  path <- tempfile(fileext = ".tsv")
  writeLines(pm, path)
  on.exit(unlink(path), add = TRUE)
  out <- lines_of(pi_run(popmap = path))
  expect_lte(max(nchar(out)), 80)
  expect_true(any(grepl("^ *Sipsey_Fork_2019 +80,000 ", out)))
  expect_true(any(grepl("Sipsey_Fork_2019 - Cahaba_River_2020", out)))
})

test_that("summary(details = TRUE) of pi_allsites() has every column in labelled blocks", {
  res <- pi_cached()
  out <- lines_of(res, details = TRUE)
  expect_gt(length(out), length(lines_of(res)))
  for (col in c("pi_nc_se", "pi_nc_lo", "pi_nc_hi", "pi_se", "dxy_se", "da", "da_nc_hi"))
    expect_true(any(grepl(paste0("\\b", col, "\\b"), out)), info = col)
  expect_true(any(grepl("^ *population +sites +pi_nc +pi_nc_se +pi_nc_lo +pi_nc_hi", out)))
  expect_true(any(grepl("^ *pair +dxy +dxy_se +dxy_lo +dxy_hi", out)))
  expect_lt(grep("^RESULTS", out)[1], grep("^FULL TABLES", out)[1])
})

## ---- het_between_pops() -----------------------------------------------------

het_run <- function(vcf = "sim.allsnps.vcf.gz", popmap = lg("popmap.tsv"), ...)
  suppressMessages(het_between_pops(lg(vcf), popmap, nboot_g2 = 50, seed = 1, ...))

test_that("the short heterozygosity summary shows the combined test as difference (SE), CI and p", {
  res <- het_run()
  out <- lines_of(res)
  expect_lt(length(out), 50)
  expect_lte(max(nchar(out)), 80)
  expect_true(any(grepl("^HETEROZYGOSITY BETWEEN POPULATIONS:", out)))
  expect_true(any(grepl("^DO THE POPULATIONS DIFFER IN HETEROZYGOSITY\\?", out)))
  expect_true(any(grepl("^DO THEY DIFFER IN INBREEDING \\(F\\)\\?", out)))      # pinned wording
  expect_true(any(grepl("Identity disequilibrium g2", out)))                     # pinned wording
  headers <- grep("^ *pair ", out, value = TRUE)
  expect_length(headers, 2L)
  expect_match(headers, "difference \\(SE\\) +95% CI +p \\(combined\\) +Hedges g")
  ## each cell is the estimate and the combined SE, formatted with the package's rule
  est_se <- RADdiversity:::.est_se
  for (tab in list(res$pairwise_tests, res$pairwise_F_tests))
    expect_true(any(grepl(est_se(tab$diff, tab$se_combined)[1], out, fixed = TRUE)))
  d <- RADdiversity:::.se_decimals(res$pairwise_tests$se_combined)
  expect_true(any(grepl(sprintf("%.*f to %.*f", d, res$pairwise_tests$ci_lo, d, res$pairwise_tests$ci_hi),
                        out, fixed = TRUE)))
  ## Welch's and Wilcoxon's p-values are for comparison and stay in the full report
  expect_false(any(grepl("p_welch|p_wilcox|Wilcoxon p|_BH", out)))
  ## no SE for the population means: comparing by overlapping intervals is the mistake this avoids
  expect_true(any(grepl("mean_het has no SE here", out)))
  expect_false(any(grepl("mean_het \\(SE\\)|\\bse\\b", grep("^ *population", out, value = TRUE))))
  expect_true(any(grepl(RADdiversity:::.coverage$het_test, out, fixed = TRUE)))
  expect_true(any(grepl("^  (ok|info|look) ", out)))
})

test_that("with 3 populations the short heterozygosity summary shows the overall test and BH p-values", {
  res <- het_run(popmap = lg("pm_three.tsv"))
  out <- lines_of(res)
  expect_lte(max(nchar(out)), 80)
  joined <- gsub("\\s+", " ", paste(out, collapse = " "))
  expect_match(joined, "Overall test across all 3 populations \\(Welch's one-way ANOVA\\): p = ")
  expect_true(any(grepl("p \\(combined, BH\\)", out)))
  expect_equal(sum(grepl("^ *g[1-3] - g[1-3] ", out)), 2L * nrow(res$pairwise_tests))  # both questions
  expect_match(joined, "read the pairs only if the overall test is significant")
})

test_that("the short heterozygosity summary flags failed libraries and g2 above 0", {
  res <- suppressWarnings(het_run(vcf = "one_dead_ind.vcf.gz"))
  out <- lines_of(res)
  expect_true(any(grepl("^  look +KEPT, but they look like failed libraries: A04", out)))   # pinned wording
  expect_false(any(grepl("no individual looks like a failed library", out)))

  res <- het_run()
  res$g2$g2_lo <- c(0.01, 0.02)
  joined <- gsub("\\s+", " ", paste(lines_of(res), collapse = " "))
  expect_match(joined, "look individuals differ in inbreeding in popA and popB \\(g2 above 0\\)")
  expect_match(joined, "report _se_combined from diversity_stats\\(\\) \\(computed by default\\)")
  res$g2$g2_lo <- c(-0.1, -0.2)
  expect_match(gsub("\\s+", " ", paste(lines_of(res), collapse = " ")),
               "ok g2 is not above 0 in any population")
})

test_that("the heterozygosity summaries flag call rate tracking heterozygosity either way", {
  res <- het_run()
  res$missingness_confound$no_variation <- FALSE
  res$missingness_confound$too_few <- FALSE
  res$missingness_confound$overall <- data.frame(r = -0.6, p_value = 0.001, n = 20,
                                                 call_rate_gap = 0)
  joined <- gsub("\\s+", " ", paste(lines_of(res), collapse = " "))
  expect_match(joined, "look individuals with more missing data look MORE heterozygous")
  expect_true(any(grepl("NEGATIVE and significant", lines_of(res, details = TRUE))))
  expect_true(any(grepl("MORE heterozygous", capture.output(print(res)))))
  res$missingness_confound$overall$r <- 0.6
  joined <- gsub("\\s+", " ", paste(lines_of(res), collapse = " "))
  expect_match(joined, "look individuals with more missing data look LESS heterozygous")
  res$missingness_confound$overall$p_value <- 0.2           # not significant either way
  joined <- gsub("\\s+", " ", paste(lines_of(res), collapse = " "))
  expect_match(joined, "ok heterozygosity does not track call rate")
})

test_that("the short heterozygosity summary keeps every row label with long popmap names", {
  pm <- long_popmap()
  on.exit(unlink(pm), add = TRUE)
  out <- lines_of(het_run(popmap = pm))
  expect_lte(max(nchar(out)), 80)
  ## The two-question tables are too wide for 80 columns with names this long, so
  ## they wrap; every block still starts with its header and its row label.
  expect_gte(sum(grepl("^ *Sipsey_Fork_2019 - Cahaba_River_2020", out)), 2L)
  expect_equal(sum(grepl("^ *Sipsey_Fork_2019 - Cahaba_River_2020", out)),
               sum(grepl("^ *pair ", out)))
  expect_true(any(grepl("^ *Sipsey_Fork_2019 +15 ", out)))
})

test_that("the short heterozygosity summary stays short with 5 populations (10 pairs)", {
  pm <- readLines(lg("popmap.tsv"))
  path <- tempfile(fileext = ".tsv")
  writeLines(paste0(sub("\t.*", "", pm), "\tp", rep(1:5, length.out = length(pm))), path)
  on.exit(unlink(path), add = TRUE)
  res <- het_run(popmap = path)
  expect_equal(nrow(res$pairwise_tests), 10L)
  out <- lines_of(res)
  expect_lte(max(nchar(out)), 80)
  expect_lt(length(out), 65)                                # 69 lines when all 10 pairs were shown
  ## at most 5 pairs per question, the most significant, and where to find the rest
  expect_equal(sum(grepl("^ *p[1-5] - p[1-5] ", out)), 10L)   # 5 rows for each of the 2 questions
  joined <- gsub("\\s+", " ", paste(out, collapse = " "))
  expect_match(joined, "Showing the 5 most significant of 10 pairs; the full table is x\\$pairwise_tests")
  expect_match(joined, "Showing the 5 most significant of 10 pairs; the full table is x\\$pairwise_F_tests")
  ## the full report still shows every pair
  full <- lines_of(res, details = TRUE)
  expect_gte(sum(grepl("^ *p[1-5] - p[1-5] ", full)), 20L)
})

test_that("summary(details = TRUE) of het_between_pops() has every test column in labelled blocks", {
  res <- het_run(popmap = lg("pm_three.tsv"))
  out <- lines_of(res, details = TRUE)
  expect_gt(length(out), length(lines_of(res)))
  for (col in c("se_combined", "ci_lo", "p_combined_BH", "p_welch", "p_wilcox_BH", "hedges_g", "g2_hi"))
    expect_true(any(grepl(paste0("\\b", col, "\\b"), out)), info = col)
  expect_true(any(grepl("^ *pair +diff +se_combined +ci_lo +ci_hi +df +p_combined +p_combined_BH", out)))
  expect_true(any(grepl("^ *pair +hedges_g +p_welch +p_welch_BH +p_wilcox +p_wilcox_BH", out)))
  expect_true(any(grepl("Interpretation notes", out)))
  expect_lt(grep("^DO THE POPULATIONS", out)[1], grep("^FULL TABLES", out)[1])
  expect_true(any(grepl("^ *g[1-3] - g[1-3] ", out)))
})

## ---- summary(x)$tables and the manuscript table -------------------------------

## Every cell of every table in `tables` is on the screen in the short view.
## This assumes nothing is cut: the short views show at most 10 pairs (5 per
## question for het_between_pops()), so call it only with runs that have fewer.
expect_tables_on_screen <- function(res, info = NULL) {
  out <- lines_of(res)
  for (nm in names(summary(res)$tables)) {
    tab <- summary(res)$tables[[nm]]
    for (col in setdiff(names(tab), c("population", "pair")))
      for (cell in unique(tab[[col]]))
        expect_true(any(grepl(cell, out, fixed = TRUE)), info = paste(info, nm, col, cell))
  }
}

test_that("summary(x)$tables holds what the short view prints, for all four functions", {
  snps <- div_run("sim.allsnps.vcf.gz", sites = 1e6)
  haps <- div_run("sim.haps.vcf.gz")
  expect_named(summary(snps)$tables, c("results", "per_site"))
  expect_named(summary(haps)$tables, "results")
  expect_identical(names(summary(snps)$tables$results),
                   c("population", "n", "Ho (SE)", "He (SE)"))
  expect_identical(names(summary(snps)$tables$per_site),
                   c("population", "Ho_autosomal (SE)", "He_autosomal (SE)"))
  expect_identical(names(summary(haps)$tables$results),
                   c("population", "n", "Fis (SE)", "Ar (SE)", "rarefied private alleles (SE)"))
  expect_named(summary(div_run("sim.allsnps.vcf.gz"))$tables, "results")   # no `sites`: no per_site
  expect_tables_on_screen(snps, "snps")
  expect_tables_on_screen(haps, "haps")

  dif <- dif_run("sim.haps.vcf.gz", beta = FALSE)
  expect_named(summary(dif)$tables, "results")
  expect_identical(names(summary(dif)$tables$results), c("pair", "FST (SE)", "D (SE)"))
  expect_tables_on_screen(dif, "dif")

  het <- het_run()
  expect_named(summary(het)$tables, c("heterozygosity", "inbreeding", "populations"))
  expect_identical(names(summary(het)$tables$heterozygosity),
                   c("pair", "difference (SE)", "95% CI", "p (combined)", "Hedges g"))
  expect_tables_on_screen(het, "het")

  pi_res <- pi_cached()
  expect_named(summary(pi_res)$tables, c("results", "between", "individual_het"))
  expect_tables_on_screen(pi_res, "pi")

  ## details = TRUE does not change what is held
  expect_identical(summary(snps, details = TRUE)$tables, summary(snps)$tables)
})

test_that("the tables hold every pair, even where the printed summary cuts the list", {
  ## six populations: 15 pairs for differentiation, 15 pairs per question for the het tests
  pm <- readLines(lg("popmap.tsv"))
  path <- tempfile(fileext = ".tsv")
  writeLines(paste0(sub("\t.*", "", pm), "\tp", rep(1:6, length.out = length(pm))), path)
  on.exit(unlink(path), add = TRUE)

  dif <- dif_run("sim.allsnps.vcf.gz", popmap = path, beta = FALSE)
  expect_equal(nrow(dif$pairwise), 15L)
  expect_equal(nrow(summary(dif)$tables$results), 16L)               # 15 pairs + "all populations"
  expect_equal(sum(grepl("^ *p[1-6] - p[1-6] ", lines_of(dif))), 10L)   # the screen shows 10

  het <- het_run(popmap = path)
  expect_equal(nrow(summary(het)$tables$heterozygosity), 15L)
  expect_equal(nrow(summary(het)$tables$inbreeding), 15L)
  expect_equal(nrow(summary(het)$tables$populations), 6L)
  expect_equal(sum(grepl("^ *p[1-6] - p[1-6] ", lines_of(het))), 10L)     # 5 per question
})

## diversity_table() replaced a 5-line recipe the docs used to repeat. The
## recipe is kept here as the reference the function must match.
plain_df <- function(x) {
  class(x) <- "data.frame"
  attr(x, "caption") <- NULL
  x
}
manuscript_table1 <- function(div_snps, div_haps) {
  snp <- summary(div_snps)$tables
  hap <- summary(div_haps)$tables
  parts  <- list(snp$results, snp$per_site, hap$results[names(hap$results) != "n"])
  table1 <- Reduce(function(a, b) merge(a, b, by = "population", sort = FALSE),
                   Filter(Negate(is.null), parts))
  table1[match(snp$results$population, table1$population), ]
}

test_that("diversity_table() holds exactly the columns the two default summaries show", {
  snps <- div_run("sim.allsnps.vcf.gz", sites = 1e6)
  haps <- div_run("sim.haps.vcf.gz")
  table1 <- diversity_table(snps, haps)
  expect_s3_class(table1, c("raddiv_table", "data.frame"))
  expect_identical(names(table1),
                   c("population", "n", "Ho (SE)", "He (SE)", "Ho_autosomal (SE)",
                     "He_autosomal (SE)", "Fis (SE)", "Ar (SE)", "rarefied private alleles (SE)"))
  ## the same cells as the old recipe, in popmap order
  reference <- manuscript_table1(snps, haps)
  expect_equal(plain_df(table1)[names(reference)], reference, ignore_attr = TRUE)
  expect_identical(table1$population, snps$per_population$population)
  ## every cell is printed by the summary of the run it came from
  snp_out <- lines_of(snps); hap_out <- lines_of(haps)
  for (col in c("Ho (SE)", "He (SE)", "Ho_autosomal (SE)", "He_autosomal (SE)"))
    for (cell in table1[[col]]) expect_true(any(grepl(cell, snp_out, fixed = TRUE)), info = col)
  for (col in c("Fis (SE)", "Ar (SE)", "rarefied private alleles (SE)"))
    for (cell in table1[[col]]) expect_true(any(grepl(cell, hap_out, fixed = TRUE)), info = col)
  ## none of the columns the summary leaves out is in the table
  expect_false(any(grepl("_se|_lo|_hi|_n$|priv_total|privAr|pct_poly|sites_used", names(table1))))

  ## without `sites` there are no per-site columns
  expect_identical(names(diversity_table(div_run("sim.allsnps.vcf.gz"), haps)),
                   c("population", "n", "Ho (SE)", "He (SE)", "Fis (SE)", "Ar (SE)",
                     "rarefied private alleles (SE)"))
  ## one run alone works
  expect_identical(names(diversity_table(snps = snps)),
                   c("population", "n", "Ho (SE)", "He (SE)", "Ho_autosomal (SE)", "He_autosomal (SE)"))
  expect_identical(names(diversity_table(haps = haps)),
                   c("population", "n", "Fis (SE)", "Ar (SE)", "rarefied private alleles (SE)"))

  ## rows are joined by population name: reversing the popmap order of the
  ## haplotype run does not scramble the table
  pm <- readLines(lg("popmap.tsv"))
  rev_path <- tempfile(fileext = ".tsv")
  writeLines(c(pm[grepl("popB$", pm)], pm[grepl("popA$", pm)]), rev_path)
  on.exit(unlink(rev_path), add = TRUE)
  reversed <- diversity_table(snps, div_run("sim.haps.vcf.gz", popmap = rev_path))
  expect_identical(plain_df(reversed), plain_df(table1))
})

test_that("diversity_table() stops on swapped or mismatched runs, and warns on a different n", {
  snps <- div_run("sim.allsnps.vcf.gz")
  haps <- div_run("sim.haps.vcf.gz")
  expect_error(diversity_table(), "at least one")
  expect_error(diversity_table(haps, snps), "Did you swap")
  expect_error(diversity_table(snps = list(a = 1)), "must be the result of diversity_stats")
  other <- haps
  other$per_population$population[1] <- "elsewhere"
  expect_error(diversity_table(snps, other), "Only in the haplotype run: elsewhere")
  fewer <- haps
  fewer$per_population$n[1] <- fewer$per_population$n[1] - 1L
  expect_warning(diversity_table(snps, fewer), "different individuals")
})

test_that("diversity_table()'s caption says where each column came from and which SE it holds", {
  snps <- div_run("sim.allsnps.vcf.gz", sites = 1e6)
  haps <- div_run("sim.haps.vcf.gz")
  cap <- attr(diversity_table(snps, haps), "caption")
  expect_match(cap, "per-site Ho and He are from the SNP VCF")
  expect_match(cap, "Fis, Ar and rarefied private alleles are from the haplotype VCF")
  expect_match(cap, sprintf("g = %d gene copies", haps$settings$g))
  expect_match(cap, "among RAD loci and among individuals")
  expect_match(cap, "SEs for Ar and rarefied private alleles count variation among RAD loci only")
  ## the printed table ends with the caption
  out <- capture.output(print(diversity_table(snps, haps)))
  expect_true(any(out == "Caption:"))
  expect_false(any(grepl("^ *1 ", out)))                          # no row numbers
  ## without the individual SEs, the caption says loci only
  off <- attr(diversity_table(div_run("sim.allsnps.vcf.gz", se_individuals = FALSE)), "caption")
  expect_match(off, "SEs for Ho and He count variation among RAD loci only")
  expect_false(grepl("per-site", off))
})

test_that("diversity_table() writes to CSV like a plain data frame", {
  table1 <- diversity_table(div_run("sim.allsnps.vcf.gz"), div_run("sim.haps.vcf.gz"))
  a <- tempfile(fileext = ".csv"); b <- tempfile(fileext = ".csv")
  on.exit(unlink(c(a, b)), add = TRUE)
  utils::write.csv(table1, a, row.names = FALSE)
  plain <- table1; class(plain) <- "data.frame"; attr(plain, "caption") <- NULL
  utils::write.csv(plain, b, row.names = FALSE)
  expect_identical(readLines(a), readLines(b))
})
