## het_between_pops(): which individuals and loci each part uses, the
## failed-individual warning, and the populations it cannot analyze.

## `ns` individuals per population, `n_rec` records under HWE, no missing data.
het_data <- function(ns, n_rec = 300) {
  gs <- lapply(ns, function(n) sim_genotypes(stats::runif(n_rec, 0.1, 0.9), n))
  A1 <- do.call(cbind, lapply(gs, `[[`, "A1"))
  A2 <- do.call(cbind, lapply(gs, `[[`, "A2"))
  ids <- unlist(Map(function(p, n) paste0(p, seq_len(n)), letters[seq_along(ns)], ns))
  colnames(A1) <- colnames(A2) <- ids
  list(A1 = A1, A2 = A2,
       pops = split(ids, factor(rep(toupper(letters[seq_along(ns)]), ns))))
}
blank <- function(d, rows, ids) {
  d$A1[rows, ids] <- NA
  d$A2[rows, ids] <- NA
  d
}

test_that("populations absent from many records keep every individual", {
  ## Each population is absent from a different 40% of records, so few loci
  ## pass the POOLED call-rate filter. Individuals used to be excluded on that
  ## pooled set, which removed all of them.
  set.seed(4)
  d <- het_data(c(10, 10, 10))
  for (p in d$pops) d <- blank(d, sample(300, 120), p)
  expect_no_warning(res <- het_between_pops(sim_H(d$A1, d$A2), d$pops, min_loci = 100,
                                            nboot_g2 = 10, verbose = FALSE))
  expect_equal(res$population_summary$n, c(10L, 10L, 10L))
  expect_true(all(res$population_summary$n_loci >= 100))
  expect_true(all(is.finite(res$pairwise_tests$p_welch)))
})

test_that("an individual with no call at its population's loci is left out, without NaN", {
  set.seed(2)
  d <- het_data(c(10, 10))
  ## b1 is missing wherever population B clears 90%, and typed elsewhere.
  d <- blank(d, 1:150, "b1")
  d <- blank(d, 151:300, c("b2", "b3"))
  res <- suppressWarnings(het_between_pops(sim_H(d$A1, d$A2), d$pops, min_call = 0.9,
                                           min_loci = 20, nboot_g2 = 10, verbose = FALSE))
  summary_B <- res$population_summary[res$population_summary$population == "B", ]
  expect_equal(summary_B$n, 9L)
  expect_true(all(is.finite(unlist(summary_B[, c("mean_het", "sd", "se", "min", "max")]))))
  expect_true(is.finite(res$g2$g2[res$g2$population == "B"]))
  expect_equal(res$pairwise_tests$n2, 9L)
})

test_that("a population it cannot analyze stops with a message that says why", {
  set.seed(1)
  d <- het_data(c(8, 1))
  expect_error(het_between_pops(sim_H(d$A1, d$A2), d$pops, verbose = FALSE),
               "fewer than 2 individuals: B")

  d <- het_data(c(10, 10))
  for (j in 1:300) d <- blank(d, j, paste0("b", (j + 0:1) %% 10 + 1))   # B's call rate: 80%
  expect_error(het_between_pops(sim_H(d$A1, d$A2), d$pops, min_call = 0.85, verbose = FALSE),
               "B \\(highest call rate 80%\\)")
})

test_that("a failed library is kept and named in a warning", {
  set.seed(5)
  d <- het_data(c(20, 5), n_rec = 2000)
  d <- blank(d, which(stats::runif(2000) < 0.98), "b5")
  H <- sim_H(d$A1, d$A2)
  expect_warning(res <- het_between_pops(H, d$pops, nboot_g2 = 10, verbose = FALSE),
                 "b5 \\(B: .*They are KEPT")
  expect_identical(res$settings$flagged_individuals, "b5")
  expect_true("b5" %in% res$individual_heterozygosity$sample)
  expect_output(print(res), "failed libraries and were kept: b5")
  expect_output(print(summary(res)), "KEPT, but they look like failed libraries: b5")
  ## The same check in the other functions that report per-individual values.
  expect_warning(individual_inbreeding(H, d$pops, verbose = FALSE), "b5")
  expect_warning(identity_disequilibrium(H, d$pops, nboot = 0, nperm = 0, verbose = FALSE), "b5")
})

test_that("the failed-library call rate is judged on the records where its population has data", {
  set.seed(6)
  d <- het_data(c(10, 10), n_rec = 400)
  d <- blank(d, 1:200, unlist(d$pops$B))          # B absent from half the records: not a failure
  d <- blank(d, 201:320, "b1")                     # b1: 60% missing among B's 200 records
  d <- blank(d, 201:280, "b2")                     # b2: 40% missing
  expect_warning(
    flagged <- RADdiversity:::.warn_failed_individuals(sim_H(d$A1, d$A2), d$pops, min_loci = 50),
    "b1")
  expect_identical(flagged, "b1")

  ## A dataset with fewer records than min_loci does not make everyone "failed".
  small <- het_data(c(5, 5), n_rec = 30)
  expect_no_warning(RADdiversity:::.warn_failed_individuals(sim_H(small$A1, small$A2), small$pops,
                                                             min_loci = 50))
})
