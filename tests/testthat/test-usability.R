## Behaviour that matters for interactive use: argument checks with clear
## messages, `verbose = FALSE` meaning silence, results that keep full
## precision and expose every table, and filters that chain.

fx <- function(name) test_path("fixtures", name)
haps <- function() read_stacks_vcf(fx("small.haps.vcf"), verbose = FALSE)

test_that("scalar arguments are checked with a message naming the argument", {
  vcf <- fx("small.haps.vcf")
  pm <- fx("small_popmap.tsv")
  expect_error(diversity_stats(vcf, pm, g = c(4, 6)), "`g` must be a single whole number")
  expect_error(diversity_stats(vcf, pm, g = 4.9), "`g` must be a single whole number")
  expect_error(diversity_stats(vcf, pm, g = "20"), "`g` must be a single whole number")
  expect_error(diversity_stats(vcf, pm, g = 4, boot = c("loci", "both")), "`boot` must be one of")
  expect_error(diversity_stats(vcf, pm, g = 4, nboot = -1), "`nboot`")
  expect_error(diversity_stats(vcf, pm, g = 4, se_individuals = "yes"), "TRUE or FALSE")
  expect_error(diversity_stats(vcf, pm, g = 4, seed = "a"), "`seed`")
  expect_error(identity_disequilibrium(vcf, pm, min_call = 2), "`min_call`")
  expect_error(filter_maf(haps(), min_maf = c(0.1, 0.2)), "`min_maf`")
  expect_error(filter_maf(haps(), min_maf = -3), "`min_maf`")
  expect_error(filter_max_het(haps(), max_ho = 2), "`max_ho`")
  expect_error(kinship_check(haps(), threshold = "high"), "`threshold`")
  expect_error(hwe_test(haps(), stop_after = 0), "`stop_after`")
})

test_that("verbose = FALSE silences every progress message", {
  vcf <- fx("small.haps.vcf")
  pm <- fx("small_popmap.tsv")
  expect_no_message(diversity_stats(vcf, pm, g = 4, nboot = 10, se_individuals = TRUE,
                                    verbose = FALSE))
  expect_no_message(het_between_pops(vcf, pm, min_call = 0.5, verbose = FALSE))
  expect_no_message(differentiation_stats(vcf, pm, nboot = 10, verbose = FALSE))
  expect_no_message(individual_inbreeding(vcf, pm, min_call = 0.5, verbose = FALSE))
  expect_no_message(identity_disequilibrium(vcf, pm, nboot = 10, nperm = 10, min_call = 0.5,
                                            verbose = FALSE))
  expect_no_message(hwe_test(vcf, pm, n_draws = 100, verbose = FALSE))
  expect_no_message(kinship_check(vcf, verbose = FALSE))
  expect_no_message(read_stacks_vcf(vcf, verbose = FALSE))
  expect_no_message(filter_call_rate(haps(), 0.5, popmap = pm, verbose = FALSE))
})

test_that("results keep full precision; only print() and the files round", {
  res <- diversity_stats(fx("small.haps.vcf"), fx("small_popmap.tsv"), g = 4, nboot = 0,
                         verbose = FALSE)
  he <- res$per_population$He
  expect_false(isTRUE(all.equal(he, round(he, 4), tolerance = 0)))
  printed <- capture.output(print(res))
  expect_true(any(grepl(sprintf("%.4f", round(he[1], 4)), printed, fixed = TRUE)))
})

test_that("every table of a het_between_pops() result is a list element, not hidden in attributes", {
  res <- het_between_pops(fx("small.haps.vcf"), fx("small_popmap.tsv"), min_call = 0.5,
                          verbose = FALSE)
  expect_null(attributes(res)$report)
  for (tab in c("population_summary", "overdispersion", "g2"))
    expect_s3_class(res[[tab]], "data.frame")
  expect_true(is.list(res$missingness_confound))
  expect_true(is.list(res$settings))
})

test_that("filters chain, including filter_low_conf_alt()", {
  H <- read_stacks_vcf(fx("small_ad.haps.vcf"), verbose = FALSE)
  out <- H |>
    filter_low_conf_alt(min_alt_reads = 2, verbose = FALSE) |>
    filter_max_het(max_ho = 1, verbose = FALSE) |>
    filter_thin_one_snp(verbose = FALSE)
  expect_s3_class(out, "raddiv_vcf")
  expect_true(all(c("A1", "A2", "locus", "locus_raw", "fields") %in% names(out)))
})

test_that("filter_low_conf_alt() explains that it needs the raw VCF columns", {
  H <- haps()
  H$fields <- NULL
  expect_error(filter_low_conf_alt(H, verbose = FALSE), "H\\$fields")
  expect_error(low_conf_alt_calls(H), "H\\$fields")
})

test_that("a sample in the data but not in the popmap is named when writing FSTAT", {
  H <- haps()
  expect_error(write_fstat(H, tempfile(), popmap = list(popA = H$samples[1:3]), verbose = FALSE),
               "not in `popmap`")
})
