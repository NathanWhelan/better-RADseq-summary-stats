## Simulated diploid genotypes: n individuals x L biallelic loci with allele
## frequencies p, each individual inbred with probability F at every locus
## (its two alleles identical by descent).
sim_H <- function(n, L, F = 0, seed = 1, prefix = "s") {
  set.seed(seed)
  p <- stats::runif(L, 0.1, 0.9)
  a1 <- matrix(1L + (stats::runif(L * n) < p), L, n)
  ibd <- matrix(stats::runif(L * n) < F, L, n)
  a2 <- ifelse(ibd, a1, 1L + (stats::runif(L * n) < p))
  storage.mode(a2) <- "integer"
  samp <- paste0(prefix, seq_len(n))
  dimnames(a1) <- dimnames(a2) <- list(NULL, samp)
  list(A1 = a1, A2 = a2, locus = paste0("l", seq_len(L)), locus_raw = paste0("l", seq_len(L)),
       alleles = replicate(L, c("A", "C"), simplify = FALSE), n_alleles = rep(2L, L),
       samples = samp)
}

test_that("F is ~0 under random mating and ~F under inbreeding", {
  H0 <- sim_H(30, 3000, F = 0)
  f0 <- individual_inbreeding(H0, list(pop = H0$samples))$F
  expect_lt(abs(mean(f0)), 0.02)
  H2 <- sim_H(30, 3000, F = 0.2, seed = 2)
  f2 <- individual_inbreeding(H2, list(pop = H2$samples))$F
  expect_lt(abs(mean(f2) - 0.2), 0.03)
})

test_that("an individual heterozygous everywhere has F < 0", {
  H <- sim_H(20, 500)
  H$A1[, 1] <- 1L; H$A2[, 1] <- 2L
  f <- individual_inbreeding(H, list(pop = H$samples))
  expect_lt(f$F[f$sample == "s1"], 0)
})

test_that("with complete data, the population mean of F equals diversity_stats() FIS", {
  ## E_i is then the same for everyone, so mean(1 - O_i/E) = 1 - sum(Ho)/sum(Hs).
  A <- sim_H(12, 400, F = 0.1, seed = 3, prefix = "a")
  B <- sim_H(9, 400, F = 0.0, seed = 4, prefix = "b")
  H <- A; H$A1 <- cbind(A$A1, B$A1); H$A2 <- cbind(A$A2, B$A2); H$samples <- c(A$samples, B$samples)
  pops <- list(popA = A$samples, popB = B$samples)
  pm <- tempfile(); writeLines(paste0(H$samples, "\t", rep(names(pops), lengths(pops))), pm)
  f <- individual_inbreeding(H, pops)
  fis <- suppressMessages(diversity_stats(H, pm, g = 4, nboot = 0, stem = "f"))$per_population
  ## diversity_stats() rounds FIS to 4 decimals, so compare on an absolute
  ## scale (a relative tolerance fails for the population whose FIS is ~0).
  for (p in names(pops))
    expect_lt(abs(mean(f$F[f$population == p]) - fis$Fis[fis$population == p]), 6e-5)
})

test_that("het_between_pops() tests F between populations and finds an inbreeding difference", {
  A <- sim_H(15, 800, F = 0.3, seed = 5, prefix = "a")
  B <- sim_H(15, 800, F = 0.0, seed = 6, prefix = "b")
  H <- A; H$A1 <- cbind(A$A1, B$A1); H$A2 <- cbind(A$A2, B$A2); H$samples <- c(A$samples, B$samples)
  pm <- tempfile(); writeLines(paste0(H$samples, "\t", rep(c("inbred", "outbred"), each = 15)), pm)
  res <- suppressMessages(het_between_pops(H, pm, stem = "f"))
  expect_true(all(c("sample", "heterozygosity", "F") %in% names(res$individual_heterozygosity)))
  ft <- res$pairwise_F_tests
  expect_identical(names(ft), names(res$pairwise_tests))
  expect_gt(ft$diff, 0.15)                 # inbred minus outbred F
  expect_lt(ft$p_welch, 1e-4)
  expect_true(any(grepl("INBREEDING", capture.output(print(res)))))
})

test_that("individual_inbreeding() validates its inputs", {
  H <- sim_H(5, 50)
  expect_error(individual_inbreeding(H, list(H$samples)), "named list")
  expect_error(individual_inbreeding(H, list(pop = H$samples), min_call = 2), "min_call")
})
