fx <- function(name) test_path("fixtures", name)

test_that("diversity_stats() returns the expected structure on a small fixture", {
  outdir <- tempfile("raddiversity-test-")
  dir.create(outdir)
  on.exit(unlink(outdir, recursive = TRUE), add = TRUE)
  res <- suppressMessages(capture.output(
    result <- diversity_stats(fx("small.haps.vcf"), fx("small_popmap.tsv"),
                               g = 4, nboot = 50, outdir = outdir)
  ))
  expect_named(result, c("per_population", "richness", "autosomal", "estimator_comparison",
                         "he_difference", "fis_by_call_rate", "settings"))
  expect_setequal(result$per_population$population, c("popA", "popB"))
  expect_true(all(c("Ho", "He", "Fis", "pct_poly") %in% names(result$per_population)))
  expect_true(all(c("Ar", "privAr", "priv_total") %in% names(result$richness)))
  expect_true(file.exists(file.path(outdir, "diversity_per_population.haps.tsv")))
  expect_true(file.exists(file.path(outdir, "diversity_richness.haps.tsv")))
})

test_that("results are quiet objects: no files without outdir, tables on print(), the report on summary()", {
  vcf <- normalizePath(fx("small.haps.vcf")); pm <- normalizePath(fx("small_popmap.tsv"))
  wd <- tempfile("nowrite-"); dir.create(wd)
  old <- setwd(wd)
  on.exit({ setwd(old); unlink(wd, recursive = TRUE) }, add = TRUE)

  out <- capture.output(res <- suppressMessages(diversity_stats(vcf, pm, g = 4, nboot = 0)))
  expect_length(out, 0)                                   # computing prints nothing
  expect_s3_class(res, "raddiv_diversity")
  short <- capture.output(p <- print(res))
  expect_identical(p, res)                                # print() returns its input
  expect_true(any(grepl("Take from this haplotype VCF", short)))
  expect_lt(length(short), 40)                            # print() stays short
  full <- capture.output(print(summary(res)))
  expect_true(any(grepl("TAKE FROM THIS RUN", full)))
  expect_gt(length(full), length(short))

  het <- suppressMessages(het_between_pops(vcf, pm, min_call = 0.5))
  expect_s3_class(het, "raddiv_het")
  expect_true(any(grepl("pairwise_tests", capture.output(print(het)))))
  expect_true(any(grepl("Overdispersion check", capture.output(print(summary(het))))))

  dif <- suppressMessages(differentiation_stats(vcf, pm, nboot = 0))
  expect_s3_class(dif, "raddiv_differentiation")
  expect_true(any(grepl("DIFFERENTIATION", capture.output(print(summary(dif))))))

  expect_length(list.files(wd), 0)                        # and nothing was written
})

test_that("diversity_stats(hierfstat_check = TRUE) agrees with hierfstat", {
  skip_if_not_installed("hierfstat")
  outdir <- tempfile("raddiversity-test-")
  dir.create(outdir)
  on.exit(unlink(outdir, recursive = TRUE), add = TRUE)
  msgs <- capture_messages(expect_no_warning(invisible(capture.output(
    diversity_stats(fx("small.haps.vcf"), fx("small_popmap.tsv"), g = 4, nboot = 0,
                    outdir = outdir, hierfstat_check = TRUE)))))
  expect_true(any(grepl("cross-check vs hierfstat::basic.stats", msgs)))
  expect_true(any(grepl("cross-check vs hierfstat::allelic.richness", msgs)))
})

test_that("diversity_stats() requires g and validates its arguments", {
  expect_error(diversity_stats(fx("small.haps.vcf"), fx("small_popmap.tsv")),
               "missing")
  expect_error(
    suppressMessages(diversity_stats(fx("small.haps.vcf"), fx("small_popmap.tsv"),
                                      g = 4, min_n = 1)),
    "min_n"
  )
})

test_that("diversity_stats() rejects an inexact boot value instead of partial-matching it", {
  expect_error(
    diversity_stats(fx("small.haps.vcf"), fx("small_popmap.tsv"),
                     g = 4, nboot = 5, boot = "individual"),
    "must be one of"
  )
})

test_that("with a seed, diversity_stats() and het_between_pops() are reproducible and restore the caller's RNG state", {
  set.seed(999)
  expected <- runif(1)

  set.seed(999)
  a <- suppressMessages(diversity_stats(fx("small.haps.vcf"), fx("small_popmap.tsv"),
                                        g = 4, nboot = 20, seed = 1))
  expect_equal(runif(1), expected)
  b <- suppressMessages(diversity_stats(fx("small.haps.vcf"), fx("small_popmap.tsv"),
                                        g = 4, nboot = 20, seed = 1))
  expect_identical(a$per_population, b$per_population)

  set.seed(999)
  invisible(suppressMessages(het_between_pops(fx("small.haps.vcf"), fx("small_popmap.tsv"),
                                              min_call = 0.5, seed = 1)))
  expect_equal(runif(1), expected)
})

test_that("with seed = NULL (the default), set.seed() before the call makes the result reproducible", {
  run <- function() suppressMessages(diversity_stats(fx("small.haps.vcf"), fx("small_popmap.tsv"),
                                                     g = 4, nboot = 20))
  set.seed(5)
  a <- run()
  set.seed(5)
  b <- run()
  expect_identical(a$per_population, b$per_population)
  set.seed(6)
  c <- run()
  expect_false(identical(a$per_population$He_lo, c$per_population$He_lo))
})

test_that("diversity_stats() sites= accepts a scalar, a named vector, a data frame, or a sumstats_summary path -- all equivalent when the value is the same for every population", {
  outdir <- tempfile("raddiversity-test-")
  dir.create(outdir)
  on.exit(unlink(outdir, recursive = TRUE), add = TRUE)
  go <- function(sites) {
    out <- NULL
    invisible(capture.output(suppressMessages(
      out <- diversity_stats(fx("small.snps.vcf"), fx("small_popmap.tsv"),
                              g = 4, nboot = 0, outdir = outdir, sites = sites)
    )))
    out
  }

  r_scalar <- go(100000)
  r_vector <- go(c(popA = 100000, popB = 100000))
  r_df     <- go(data.frame(population = c("popA", "popB"), sites = c(100000, 100000)))
  r_path   <- go(fx("sumstats_summary_small.tsv"))

  expect_false(is.null(r_scalar$autosomal))
  expect_identical(r_scalar$autosomal, r_vector$autosomal)
  expect_identical(r_scalar$autosomal, r_df$autosomal)
  expect_identical(r_scalar$autosomal, r_path$autosomal)
  expect_true(all(c("sites_used", "Ho_autosomal", "He_autosomal") %in% names(r_scalar$autosomal)))
  expect_true(file.exists(file.path(outdir, "diversity_autosomal.snps.tsv")))
})

test_that("diversity_stats() sites= as a per-population vector gives each population its own denominator", {
  outdir <- tempfile("raddiversity-test-")
  dir.create(outdir)
  on.exit(unlink(outdir, recursive = TRUE), add = TRUE)
  out <- NULL
  invisible(capture.output(suppressMessages(
    out <- diversity_stats(fx("small.snps.vcf"), fx("small_popmap.tsv"), g = 4, nboot = 0,
                            outdir = outdir, sites = c(popA = 100000, popB = 50000))
  )))
  aut <- out$autosomal
  n_rec <- 15L  # tests/testthat/fixtures/small.snps.vcf has 15 records
  ## Each population's own He (per marker, from per_population) scaled by ITS
  ## OWN sites value -- not a shared scalar, which would give both
  ## populations the same denominator regardless of their true sites.
  for (p in c("popA", "popB")) {
    he_marker <- out$per_population$He[out$per_population$population == p]
    sites_p   <- c(popA = 100000, popB = 50000)[[p]]
    expect_equal(aut$He_autosomal[aut$population == p], he_marker * n_rec / sites_p)
  }
  ## And the two populations' own values are NOT equal to each other (unlike
  ## the equal-sites case above) -- confirms the per-population denominator
  ## actually took effect rather than silently falling back to one shared value.
  expect_false(isTRUE(all.equal(aut$He_autosomal[aut$population == "popA"],
                                 aut$He_autosomal[aut$population == "popB"])))
})

test_that("diversity_stats() sites= flags an implausible per-population value as NA rather than erroring the whole run", {
  outdir <- tempfile("raddiversity-test-")
  dir.create(outdir)
  on.exit(unlink(outdir, recursive = TRUE), add = TRUE)
  out <- NULL
  invisible(capture.output(suppressMessages(
    out <- diversity_stats(fx("small.snps.vcf"), fx("small_popmap.tsv"), g = 4, nboot = 0,
                            outdir = outdir, sites = c(popA = 100000, popB = 5))
  )))
  aut <- out$autosomal
  expect_false(is.na(aut$He_autosomal[aut$population == "popA"]))
  expect_true(is.na(aut$He_autosomal[aut$population == "popB"]))
})

test_that("diversity_stats() sites= errors when a population has no value, naming it", {
  expect_error(
    suppressMessages(diversity_stats(fx("small.snps.vcf"), fx("small_popmap.tsv"),
                                      g = 4, nboot = 0, sites = c(popA = 100000))),
    "popB"
  )
})

test_that("diversity_stats() sites= rejects a bad shape before the VCF is even read", {
  expect_error(diversity_stats("nonexistent.vcf", fx("small_popmap.tsv"), g = 4, sites = -5),
               "non-negative")
  expect_error(diversity_stats("nonexistent.vcf", fx("small_popmap.tsv"), g = 4,
                                sites = "no/such/file.tsv"),
               "not found")
})

test_that("diversity_stats() ignores sites= entirely (autosomal NULL) on a haplotype VCF", {
  outdir <- tempfile("raddiversity-test-")
  dir.create(outdir)
  on.exit(unlink(outdir, recursive = TRUE), add = TRUE)
  out <- NULL
  invisible(capture.output(suppressMessages(
    out <- diversity_stats(fx("small.haps.vcf"), fx("small_popmap.tsv"), g = 4, nboot = 0,
                            outdir = outdir, sites = 100000)
  )))
  expect_null(out$autosomal)
})

test_that("per-site values scale each population by its own variant records (Stacks' -r pruning)", {
  ## Three populations; population C is blanked at 30% of records, as Stacks'
  ## -r does with -p below the number of populations. Stacks' `Sites` for C
  ## then excludes those records, so C's numerator must exclude them too.
  set.seed(42)
  L <- 600
  ids <- paste0("i", 1:30)
  pops <- list(A = ids[1:10], B = ids[11:20], C = ids[21:30])
  g <- sim_genotypes(stats::runif(L, 0.05, 0.95), 30)
  A1 <- g$A1; A2 <- g$A2
  colnames(A1) <- colnames(A2) <- ids
  pruned <- sample(L, 0.3 * L)
  A1[pruned, pops$C] <- NA
  A2[pruned, pops$C] <- NA
  H <- sim_H(A1, A2, records_per_locus = 3L)
  sites <- 5e4 + c(A = L, B = L, C = L - length(pruned))

  res <- diversity_stats(H, pops, g = 16, nboot = 0, sites = sites, verbose = FALSE)
  counts <- RADdiversity:::.pop_counts(H, pops)
  typed <- RADdiversity:::.typed_by_pop(H, pops)
  hets <- RADdiversity:::.het_by_pop(H, pops, seq_len(L))
  sum_hs <- vapply(1:3, function(k) {
    hs <- hs_nei_chesser(rowSums((counts[[k]] / rowSums(counts[[k]]))^2),
                         hets[, k] / typed[, k], typed[, k])
    sum(hs, na.rm = TRUE)
  }, numeric(1))
  expect_equal(res$autosomal$variant_records, c(L, L, L - length(pruned)))
  expect_equal(res$autosomal$He_autosomal, unname(sum_hs / sites), tolerance = 1e-12)
})

test_that("per-site values warn when records were removed after reading", {
  snps <- fx("small.snps.vcf")
  pm <- fx("small_popmap.tsv")
  sumstats <- fx("sumstats_summary_small.tsv")
  H <- read_stacks_vcf(snps, verbose = FALSE)
  expect_identical(H$n_records_read, nrow(H$A1))
  H_f <- filter_call_rate(H, min_call = 1, verbose = FALSE)
  skip_if(nrow(H_f$A1) == nrow(H$A1), "fixture has no incomplete record to filter")
  expect_identical(H_f$n_records_read, nrow(H$A1))
  expect_warning(diversity_stats(H_f, pm, g = 4, nboot = 0, sites = 1e5, verbose = FALSE),
                 "removed after reading")
  ## With the sumstats file, its unfiltered Variant_Sites is used, and the
  ## mismatch with the VCF is reported.
  expect_warning(res <- diversity_stats(H_f, pm, g = 4, nboot = 0, sites = sumstats,
                                        verbose = FALSE),
                 "Variant_Sites")
  expect_equal(res$autosomal$variant_records,
               read_sumstats_summary(sumstats)$all_positions$variant_sites)
  ## The unfiltered VCF with its own sumstats file: no warning.
  expect_no_warning(diversity_stats(H, pm, g = 4, nboot = 0, sites = sumstats, verbose = FALSE))
})

test_that("fis_by_call_rate is flat without dropout and rises with it", {
  sim <- function(dropout, seed) {
    set.seed(seed)
    L <- 3000; n <- 30
    g <- sim_genotypes(stats::runif(L, 0.1, 0.9), n)
    A1 <- g$A1; A2 <- g$A2
    samp <- c(paste0("a", 1:15), paste0("b", 1:15))
    dimnames(A1) <- dimnames(A2) <- list(NULL, samp)
    if (dropout) {
      ## Half the loci carry a null allele: heterozygotes there are called
      ## homozygous half the time, and 30% of genotypes go missing.
      bad <- matrix(rep(stats::runif(L) < 0.5, n), L, n)
      het <- A1 != A2
      to_hom <- bad & het & matrix(stats::runif(L * n) < 0.5, L, n)
      A2[to_hom] <- A1[to_hom]
      gone <- bad & matrix(stats::runif(L * n) < 0.3, L, n)
    } else {
      gone <- matrix(stats::runif(L * n) < 0.1, L, n)      # missing at random
    }
    A1[gone] <- NA; A2[gone] <- NA
    pops <- list(popA = samp[1:15], popB = samp[16:30])
    diversity_stats(sim_H(A1, A2), pops, g = 20, nboot = 0, verbose = FALSE)$fis_by_call_rate
  }
  clean <- sim(FALSE, 1)
  expect_true(all(c("population", "call_rate", "n_records", "He", "Fis", "Fis_se") %in% names(clean)))
  expect_lt(max(abs(clean$Fis[clean$n_records > 100])), 0.05)

  dirty <- sim(TRUE, 2)
  a <- dirty[dirty$population == "popA", ]
  expect_gt(a$Fis[a$call_rate == "<75%"], a$Fis[a$call_rate == "100%"] + 0.1)
})

test_that("diversity_stats() notices data that were filtered by minor allele count", {
  set.seed(5)
  L <- 2000; n <- 30
  ## Neutral-like frequencies: many rare variants.
  g <- sim_genotypes(pmin(stats::rbeta(L, 0.2, 2), 0.5), n)
  samp <- c(paste0("a", 1:15), paste0("b", 1:15))
  dimnames(g$A1) <- dimnames(g$A2) <- list(NULL, samp)
  H <- sim_H(g$A1, g$A2)
  H <- filter_mac(H, min_mac = 1, verbose = FALSE)          # variable records only
  pops <- list(popA = samp[1:15], popB = samp[16:30])

  raw <- diversity_stats(H, pops, g = 20, nboot = 0, verbose = FALSE)$settings$prior_filters
  expect_false(raw$looks_mac_filtered)
  expect_gt(raw$rare_share, 0.1)

  filtered_H <- filter_mac(H, min_mac = 3, verbose = FALSE)
  res <- diversity_stats(filtered_H, pops, g = 20, nboot = 0, verbose = FALSE)
  pf <- res$settings$prior_filters
  expect_true(pf$looks_mac_filtered)
  expect_equal(pf$min_allele_count, 3L)
  ## Described in summary() as information, not raised as a warning.
  expect_true(any(grepl("filtered by minor allele count", capture.output(summary(res)))))
})

test_that("het_between_pops() returns the expected structure on a small fixture", {
  outdir <- tempfile("raddiversity-test-")
  dir.create(outdir)
  on.exit(unlink(outdir, recursive = TRUE), add = TRUE)
  invisible(capture.output(
    result <- suppressMessages(
      het_between_pops(fx("small.haps.vcf"), fx("small_popmap.tsv"),
                        min_call = 0.5, outdir = outdir)
    )
  ))
  expect_named(result, c("individual_heterozygosity", "population_summary", "omnibus",
                         "pairwise_tests", "pairwise_F_tests", "overdispersion", "g2",
                         "missingness_confound", "settings"))
  expect_null(result$omnibus)                      # two populations: no overall test
  expect_equal(nrow(result$individual_heterozygosity), 7L)
  expect_true(all(c("p_welch", "p_wilcox", "hedges_g") %in% names(result$pairwise_tests)))
  expect_true(file.exists(file.path(outdir, "individual_heterozygosity.haps.tsv")))
  expect_true(file.exists(file.path(outdir, "het_between_pops_tests.haps.tsv")))
})

test_that("het_between_pops() runs on the SNP and haplotype VCFs into one outdir without overwriting", {
  ## Regression: the output names used to be fixed, so the second run
  ## silently replaced the first run's files.
  outdir <- tempfile("raddiversity-test-")
  dir.create(outdir)
  on.exit(unlink(outdir, recursive = TRUE), add = TRUE)
  ## The same fixture under the two Stacks file names -- only the name matters.
  vcfs <- file.path(outdir, c("populations.snps.vcf", "populations.haps.vcf"))
  file.copy(fx("small.haps.vcf"), vcfs)
  for (v in vcfs)
    invisible(capture.output(suppressMessages(
      het_between_pops(v, fx("small_popmap.tsv"), min_call = 0.5, outdir = outdir))))
  expect_true(all(file.exists(file.path(outdir, c(
    "individual_heterozygosity.snps.tsv", "individual_heterozygosity.haps.tsv",
    "het_between_pops_tests.snps.tsv", "het_between_pops_tests.haps.tsv")))))
})
