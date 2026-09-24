fx <- function(name) test_path("fixtures", name)

## The individual standard errors (se_individuals) are on by default, but the
## toy fixtures have only 3-4 individuals per population, so they warn there.
## Most tests here are not about those SEs, so this wrapper runs them with the
## option off; a test that needs them passes se_individuals = TRUE. The default
## itself is tested in test-jackknife-individuals.R.
diversity_stats <- function(..., se_individuals = FALSE)
  RADdiversity::diversity_stats(..., se_individuals = se_individuals)

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
  brief <- capture.output(print(summary(res)))            # summary() is short by default ...
  expect_true(any(grepl("TAKE FROM THIS RUN", brief)))
  expect_lt(length(brief), 45)
  full <- capture.output(print(summary(res, details = TRUE)))   # ... details = TRUE has it all
  expect_true(any(grepl("TAKE FROM THIS RUN", full)))
  expect_gt(length(full), length(short))
  expect_gt(length(full), 2 * length(brief))

  het <- suppressMessages(het_between_pops(vcf, pm, min_call = 0.5))
  expect_s3_class(het, "raddiv_het")
  expect_true(any(grepl("pairwise_tests", capture.output(print(het)))))
  expect_true(any(grepl("Identity disequilibrium g2", capture.output(print(summary(het))))))

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

test_that("per-site values carry the same uncertainty columns as Ho and He, scaled by records / sites", {
  ## Ho_autosomal = Ho * records / sites, a fixed multiplier per population, so
  ## every SE and bootstrap bound is the Ho (or He) one times that multiplier.
  set.seed(11)
  g <- sim_genotypes(stats::runif(600, 0.1, 0.9), 20)
  samp <- c(paste0("a", 1:10), paste0("b", 1:10))
  A1 <- g$A1; A2 <- g$A2; dimnames(A1) <- dimnames(A2) <- list(NULL, samp)
  H <- sim_H(A1, A2, records_per_locus = 3L)
  pops <- list(A = samp[1:10], B = samp[11:20])
  run <- function(...) diversity_stats(H, pops, g = 10, sites = c(A = 40000, B = 50000),
                                       verbose = FALSE, ...)
  res <- run(nboot = 100, se_individuals = TRUE)
  pp <- res$per_population
  aut <- res$autosomal
  multiplier <- aut$variant_records / aut$sites_used
  suffixes <- c("", "_se", "_se_ind", "_se_combined", "_lo", "_hi")
  for (stat in c("Ho", "He"))
    for (suffix in suffixes) {
      col <- paste0(stat, "_autosomal", suffix)
      expect_true(col %in% names(aut), info = col)
      expect_equal(aut[[col]], pp[[paste0(stat, suffix)]] * multiplier, info = col)
    }
  expect_identical(names(aut)[4:9], paste0("Ho_autosomal", suffixes))   # estimate, then its uncertainty
  expect_equal(aut$He_autosomal, autosomal_het(pp$He, aut$variant_records, aut$sites_used))

  ## Without the individual SEs there are no _se_ind / _se_combined columns;
  ## without the bootstrap the bounds are NA, as they are for Ho and He.
  plain <- run(nboot = 0, se_individuals = FALSE)$autosomal
  expect_false(any(grepl("_se_ind|_se_combined", names(plain))))
  expect_true(all(is.na(plain$Ho_autosomal_lo)) && all(is.na(plain$He_autosomal_hi)))
  expect_true(all(is.finite(plain$Ho_autosomal_se)))

  ## A population whose sites value is implausible (fewer sites than variant
  ## records) gets NA in every per-site column and does not touch the others.
  odd <- diversity_stats(H, pops, g = 10, sites = c(A = 40000, B = 5), nboot = 0,
                         se_individuals = TRUE, verbose = FALSE)$autosomal
  per_site <- grep("_autosomal", names(odd), value = TRUE)
  expect_true(all(is.na(unlist(odd[odd$population == "B", per_site]))))
  expect_true(all(is.finite(unlist(odd[odd$population == "A",
                                       setdiff(per_site, c("Ho_autosomal_lo", "Ho_autosomal_hi",
                                                           "He_autosomal_lo", "He_autosomal_hi"))]))))
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

test_that("a multi-allelic SNP does not turn a SNP VCF into haplotype data", {
  ## GATK and freebayes can emit sites such as A -> C,T. Every allele is still
  ## one base, so every record is one site and per-site values are given.
  set.seed(3)
  n <- 6
  samp <- c(paste0("a", 1:n), paste0("b", 1:n))
  gt <- function() sample(c("0/0", "0/1", "1/1"), 2 * n, TRUE, prob = c(0.4, 0.4, 0.2))
  lines <- c("##fileformat=VCFv4.2",
             paste(c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", samp),
                   collapse = "\t"))
  for (i in 1:60)
    lines <- c(lines, paste(c("chr1", i * 2000, ".", "A", "C", ".", "PASS", ".", "GT", gt()),
                            collapse = "\t"))
  lines <- c(lines, paste(c("chr1", 999999, ".", "A", "C,T", ".", "PASS", ".", "GT",
                            "0/2", gt()[-1]), collapse = "\t"))
  vcf <- tempfile(fileext = ".vcf"); writeLines(lines, vcf)
  pm <- list(A = samp[1:n], B = samp[n + 1:n])
  H <- read_stacks_vcf(vcf, verbose = FALSE)
  expect_false(RADdiversity:::.is_haplotype_H(H))
  res <- diversity_stats(H, pm, g = 4, nboot = 0, sites = 1e5, verbose = FALSE)
  expect_false(res$settings$is_haplotype)
  expect_false(is.null(res$autosomal))
  ## The same record with multi-base (haplotype) alleles is haplotype data.
  H$alleles[[61]] <- c("AG", "CT", "TC")
  expect_true(RADdiversity:::.is_haplotype_H(H))
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
  w <- testthat::capture_warnings(
    diversity_stats(H_f, pm, g = 4, nboot = 0, sites = 1e5, verbose = FALSE))
  expect_true(any(grepl("removed after reading", w)))
  ## With the sumstats file, only its Sites is used: the SNP count comes from
  ## the data, so the result equals passing those Sites as numbers.
  all_pos <- read_sumstats_summary(sumstats)$all_positions
  res_file <- suppressWarnings(diversity_stats(H_f, pm, g = 4, nboot = 0, sites = sumstats,
                                               verbose = FALSE))
  res_num <- suppressWarnings(diversity_stats(H_f, pm, g = 4, nboot = 0, verbose = FALSE,
                                              sites = stats::setNames(all_pos$sites,
                                                                      all_pos$population)))
  expect_equal(res_file$autosomal, res_num$autosomal)
  expect_equal(res_file$autosomal$variant_records,
               unname(colSums(RADdiversity:::.typed_by_pop(H_f, read_popmap(pm, verbose = FALSE)) > 0)))
  ## The unfiltered VCF with its own sumstats file: no warning.
  expect_no_warning(diversity_stats(H, pm, g = 4, nboot = 0, sites = sumstats, verbose = FALSE))
})

test_that("per-site values do not inflate after a MAC filter (Variant_Sites is not used)", {
  ## Rare SNPs have low He. Scaling the He of the SNPs kept by the unfiltered
  ## SNP count would treat the removed SNPs as average and inflate per-site He.
  set.seed(3)
  L <- 2000; n <- 10
  ids <- paste0("i", 1:(2 * n))
  pops <- list(A = ids[1:n], B = ids[n + 1:n])
  p <- c(stats::runif(L / 2, 0.2, 0.8), stats::runif(L / 2, 0.01, 0.04))
  g <- sim_genotypes(p, 2 * n)
  colnames(g$A1) <- colnames(g$A2) <- ids
  H <- sim_H(g$A1, g$A2)
  H$n_records_read <- L
  H_f <- filter_mac(H, min_mac = 3, verbose = FALSE)
  full <- diversity_stats(H, pops, g = 4, nboot = 0, sites = 1e6, verbose = FALSE)
  filtered <- suppressWarnings(diversity_stats(H_f, pops, g = 4, nboot = 0, sites = 1e6,
                                               verbose = FALSE))
  ## The filter can only remove He, never add it.
  expect_true(all(filtered$autosomal$He_autosomal <= full$autosomal$He_autosomal))
})

test_that("per-site values warn when whole RAD loci were removed, or the VCF was thinned", {
  snps <- fx("small.snps.vcf")
  pm <- fx("small_popmap.tsv")
  sumstats <- fx("sumstats_summary_small.tsv")
  H <- read_stacks_vcf(snps, verbose = FALSE)
  loci <- unique(H$locus_raw)
  H_drop <- RADdiversity:::.subset_H(H, H$locus_raw != loci[1])
  w <- testthat::capture_warnings(
    diversity_stats(H_drop, pm, g = 4, nboot = 0, sites = sumstats, verbose = FALSE))
  expect_true(any(grepl("whole RAD loci were removed", w)))
  ## Removed without a filter_*() function, so the reason is not known.
  expect_true(any(grepl("not by a filter_*() function", w, fixed = TRUE)))

  ## A VCF thinned to one SNP per RAD locus BEFORE reading (as with Stacks'
  ## --write-single-snp), paired with the unthinned Stacks summary file.
  set.seed(8)
  L <- 300
  ids <- paste0("i", 1:12)
  pops <- list(A = ids[1:6], B = ids[7:12])
  g <- sim_genotypes(stats::runif(L, 0.2, 0.8), 12)
  colnames(g$A1) <- colnames(g$A2) <- ids
  H_thin <- filter_thin_one_snp(sim_H(g$A1, g$A2, records_per_locus = 3L), verbose = FALSE)
  H_thin$n_records_read <- nrow(H_thin$A1)
  H_thin$n_loci_read <- length(unique(H_thin$locus_raw))
  row_all <- function(p) paste(c(p, 0, 1e5, L, L, 1, rep(0, 24)), collapse = "\t")
  row_var <- function(p) paste(c(p, 0, rep(0, 24)), collapse = "\t")
  trip <- paste(as.vector(rbind(c("Num_Indv", "P", "Obs_Het", "Obs_Hom", "Exp_Het", "Exp_Hom",
                                  "Pi", "Fis"), "Var", "StdErr")), collapse = "\t")
  f <- tempfile(fileext = ".tsv")
  writeLines(c("# Variant positions", paste0("# Pop ID\tPrivate\t", trip), row_var("A"), row_var("B"),
               "# All positions (variant and fixed)",
               paste0("# Pop ID\tPrivate\tSites\tVariant_Sites\tPolymorphic_Sites\t%Polymorphic_Loci\t", trip),
               row_all("A"), row_all("B")), f)
  w <- testthat::capture_warnings(
    diversity_stats(H_thin, pops, g = 4, nboot = 0, sites = f, verbose = FALSE))
  expect_true(any(grepl("looks thinned", w)))
})

test_that("RAD loci emptied by a MAC filter give no loci warning; loci removed by call rate do", {
  ex <- function(f) system.file("extdata", f, package = "RADdiversity")
  pm <- ex("example_popmap.tsv")
  sumstats <- ex("example.sumstats_summary.tsv")
  H <- read_stacks_vcf(ex("example.snps.vcf.gz"), verbose = FALSE)
  H_mac <- filter_mac(H, min_mac = 3, verbose = FALSE)
  expect_gt(H_mac$filter_log$loci_removed, 0L)        # the filter did empty some loci
  w <- testthat::capture_warnings(
    diversity_stats(H_mac, pm, g = 10, nboot = 0, sites = sumstats, verbose = FALSE))
  expect_true(any(grepl("records were removed after reading (by filter_mac)", w, fixed = TRUE)))
  ## Their sites were still sequenced, so Stacks' Sites is right: no loci warning.
  expect_false(any(grepl("whole RAD loci", w)))

  H_call <- filter_call_rate(H_mac, min_call = 0.95, verbose = FALSE)
  w <- testthat::capture_warnings(
    res <- diversity_stats(H_call, pm, g = 10, nboot = 0, sites = sumstats, verbose = FALSE))
  loci_warning <- grep("whole RAD loci", w, value = TRUE)
  expect_length(loci_warning, 1L)
  expect_true(grepl(sprintf("%d whole RAD loci were removed after reading by filter_call_rate.",
                            H_call$filter_log$loci_removed[2]), loci_warning, fixed = TRUE))
  ## The full report lists the filters.
  report <- capture.output(print(summary(res, details = TRUE)))
  expect_true(any(grepl("2. filter_call_rate (min_call = 0.95, pooled)", report, fixed = TRUE)))
})

test_that("per-site values from data thinned in R warn to use the unthinned data", {
  set.seed(9)
  ids <- paste0("i", 1:12)
  g <- sim_genotypes(stats::runif(300, 0.2, 0.8), 12)
  colnames(g$A1) <- colnames(g$A2) <- ids
  H <- sim_H(g$A1, g$A2, records_per_locus = 3L)
  H$n_records_read <- nrow(H$A1)
  H_thin <- filter_thin_one_snp(H, verbose = FALSE)
  w <- testthat::capture_warnings(
    diversity_stats(H_thin, list(A = ids[1:6], B = ids[7:12]), g = 4, nboot = 0, sites = 1e5,
                    verbose = FALSE))
  expect_true(any(grepl("compute them from the unthinned SNP data", w, fixed = TRUE)))
})

test_that("a population with no record at min_n stops, naming it", {
  expect_error(diversity_stats(fx("small.snps.vcf"), fx("small_popmap.tsv"), g = 4, nboot = 0,
                               min_n = 4, verbose = FALSE),
               "population\\(s\\): popB \\(at most 3 genotyped at any record, of 3 individuals\\)")
})

test_that("`sites` on a haplotype VCF is reported as needing the SNP VCF, whatever its value", {
  for (s in c(5, 1e6)) {
    res <- diversity_stats(fx("small.haps.vcf"), fx("small_popmap.tsv"), g = 4, nboot = 0,
                           sites = s, verbose = FALSE)
    expect_null(res$autosomal)
    printed <- capture.output(print(res))
    expect_true(any(grepl("per-site values need the SNP VCF", printed)))
    expect_false(any(grepl("smaller than the number of variant records", printed)))
    expect_true(any(grepl("HAPLOTYPE VCF, so the autosomal",
                          capture.output(print(summary(res, details = TRUE))))))
  }
})

test_that("the 'Ar is capped at 2' note is for SNP VCFs, not haplotype VCFs", {
  ## Both toy files have mean Ar below 2. On the SNP VCF that is the biallelic
  ## ceiling. On the haplotype VCF it is only a low mean, and the note would
  ## tell the user to run the file they already ran.
  capped <- function(vcf) {
    res <- diversity_stats(fx(vcf), fx("small_popmap.tsv"), g = 4, nboot = 0, verbose = FALSE)
    expect_lt(max(res$richness$Ar), 2.001)
    any(grepl("Ar is capped at 2", capture.output(print(summary(res, details = TRUE)))))
  }
  expect_true(capped("small.snps.vcf"))
  expect_false(capped("small.haps.vcf"))
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

test_that("print() and summary() note when Ar rests on uneven records between populations", {
  set.seed(12)
  L <- 400
  g <- sim_genotypes(stats::runif(L, 0.1, 0.9), 16)
  samp <- c(paste0("a", 1:8), paste0("b", 1:8))
  A1 <- g$A1; A2 <- g$A2
  dimnames(A1) <- dimnames(A2) <- list(NULL, samp)
  pops <- list(popA = samp[1:8], popB = samp[9:16])
  ## print() shows the Ar notes for haplotype data only (Ar is not taken from a
  ## SNP VCF), so the records get multi-base alleles here.
  as_haps <- function(H) { H$alleles <- rep(list(c("AT", "CG")), nrow(H$A1)); H }
  even <- diversity_stats(as_haps(sim_H(A1, A2)), pops, g = 12, nboot = 0, verbose = FALSE)
  expect_false(any(grepl("Ar rests on under 90%", capture.output(print(even)))))

  ## popB loses 3 of its 8 individuals at a quarter of the records: 10 gene
  ## copies there, below g = 12, so its Ar rests on ~75% of popA's records.
  gone <- which(stats::runif(L) < 0.25)
  A1[gone, samp[9:11]] <- NA; A2[gone, samp[9:11]] <- NA
  uneven <- diversity_stats(as_haps(sim_H(A1, A2)), pops, g = 12, nboot = 0, verbose = FALSE)
  expect_lt(uneven$richness$Ar_n[2], 0.9 * uneven$richness$Ar_n[1])
  expect_true(any(grepl("Ar rests on under 90%", capture.output(print(uneven)))))
  snp_run <- diversity_stats(sim_H(A1, A2), pops, g = 12, nboot = 0, verbose = FALSE)
  expect_false(any(grepl("Ar rests on under 90%", capture.output(print(snp_run)))))
  expect_true(any(grepl("Ar_n differs by more than 10%",
                        capture.output(print(summary(uneven, details = TRUE))))))

  ## The short summary shows it as a `look` check (for haplotype data only: on a
  ## SNP VCF Ar is not shown, so its coverage is not either).
  Hh <- sim_H(A1, A2)
  Hh$alleles <- rep(list(c("AT", "CG")), L)               # multi-base alleles: haplotype data
  brief <- capture.output(print(summary(diversity_stats(Hh, pops, g = 12, nboot = 0,
                                                        verbose = FALSE))))
  brief <- gsub("\\s+", " ", paste(brief, collapse = " "))
  expect_match(brief, "look Ar is averaged over [0-9,]+ loci in popA but [0-9,]+ in popB")
  expect_match(brief, "populations are compared over partly different loci")
})

test_that("fis_by_call_rate's SE uses only the loci in each call-rate group", {
  ## 600 loci, all fully typed except 5 with one missing individual in popA:
  ## the "90-99%" group of popA has 5 loci, and its SE must be the jackknife
  ## over those 5 alone, not over all 600 (which would inflate it by ~12%).
  set.seed(11)
  L <- 600; n <- 10
  g <- sim_genotypes(stats::runif(L, 0.2, 0.8), 2 * n)
  samp <- c(paste0("a", 1:n), paste0("b", 1:n))
  dimnames(g$A1) <- dimnames(g$A2) <- list(NULL, samp)
  g$A1[1:5, 1] <- NA; g$A2[1:5, 1] <- NA
  pops <- list(popA = samp[1:n], popB = samp[n + 1:n])
  tab <- diversity_stats(sim_H(g$A1, g$A2), pops, g = 10, nboot = 0,
                         verbose = FALSE)$fis_by_call_rate
  row <- tab[tab$population == "popA" & tab$call_rate == "90-99%", ]
  expect_equal(row$n_records, 5L)
  ## The same group computed on its own 5 loci:
  few <- sim_H(g$A1[1:5, ], g$A2[1:5, ])
  alone <- diversity_stats(few, pops, g = 10, nboot = 0, verbose = FALSE)$per_population
  expect_equal(row$Fis_se, alone$Fis_se[1], tolerance = 1e-10)
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

  ## A VCF sample outside the popmap (an outgroup) carrying singletons at
  ## every record does not hide the filter: the check uses popmap samples.
  out1 <- matrix(1L, nrow(filtered_H$A1), 1, dimnames = list(NULL, "outgroup"))
  out2 <- out1
  out2[, 1] <- ifelse(filtered_H$A1[, 1] == 1L, 2L, 1L)   # heterozygous everywhere
  with_outgroup <- filtered_H
  with_outgroup$A1 <- cbind(filtered_H$A1, out1)
  with_outgroup$A2 <- cbind(filtered_H$A2, out2)
  with_outgroup$samples <- c(filtered_H$samples, "outgroup")
  pf2 <- diversity_stats(with_outgroup, pops, g = 20, nboot = 0,
                         verbose = FALSE)$settings$prior_filters
  expect_true(pf2$looks_mac_filtered)
  expect_equal(pf2$n_samples, 30L)
  expect_true(any(grepl("over the 30 individuals in the popmap",
                        capture.output(summary(res, details = TRUE)))))
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
                         "pairwise_tests", "pairwise_F_tests", "g2",
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

test_that("the rarefied private alleles' SE and interval are privAr's times privAr_n", {
  ## priv_total is privAr added up over the privAr_n loci, the same loci for
  ## every population, so its uncertainty is privAr's scaled by that count
  ## (the SE the simulations checked), not the resampled sum's.
  ex <- function(f) system.file("extdata", f, package = "RADdiversity")
  res <- diversity_stats(ex("example.haps.vcf.gz"), ex("example_popmap.tsv"), g = 20,
                         nboot = 100, seed = 1, verbose = FALSE)
  r <- res$richness
  expect_lt(max(r$privAr_n), res$settings$n_records_used)   # some loci fall below g
  expect_equal(r$priv_total, r$privAr * r$privAr_n)
  expect_equal(r$priv_total_se, r$privAr_se * r$privAr_n)
  expect_equal(r$priv_total_lo, r$privAr_lo * r$privAr_n)
  expect_equal(r$priv_total_hi, r$privAr_hi * r$privAr_n)
})
