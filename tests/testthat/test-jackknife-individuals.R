## Two populations of simulated genotypes; F gives each individual's
## probability of carrying two identical-by-descent alleles at a locus.
sim_two_pops <- function(F, L, seed) {
  set.seed(seed)
  n <- length(F)
  p <- stats::runif(L, 0.1, 0.9)
  a1 <- matrix(1L + (stats::runif(L * n) < p), L, n)
  ibd <- matrix(stats::runif(L * n) < rep(F, each = L), L, n)
  a2 <- ifelse(ibd, a1, 1L + (stats::runif(L * n) < p))
  samp <- c(paste0("a", seq_len(n / 2)), paste0("b", seq_len(n / 2)))
  dimnames(a1) <- dimnames(a2) <- list(NULL, samp)
  pm <- tempfile()
  writeLines(paste0(samp, "\t", rep(c("popA", "popB"), each = n / 2)), pm)
  list(H = list(A1 = a1, A2 = a2, locus = paste0("l", seq_len(L)),
                locus_raw = paste0("l", seq_len(L)),
                alleles = replicate(L, c("A", "C"), simplify = FALSE),
                n_alleles = rep(2L, L), samples = samp),
       popmap = pm)
}
run_ind <- function(d, g = 4)
  suppressMessages(diversity_stats(d$H, d$popmap, g = g, nboot = 0, stem = "j",
                                   se_individuals = TRUE))

test_that("with complete data, the individual SE of Ho is exactly sd(het)/sqrt(n)", {
  ## Ho is then the mean of the individuals' heterozygosities, and the
  ## delete-one jackknife SE of a mean is sd/sqrt(n).
  d <- sim_two_pops(rep(0.1, 12), 300, seed = 1)
  res <- expect_no_warning(run_ind(d))     # no mismatch with the point estimates
  hA <- colMeans(d$H$A1[, 1:6] != d$H$A2[, 1:6])
  expect_equal(res$per_population$Ho_se_ind[1], round(sd(hA) / sqrt(6), 4))
  expect_true(all(c("Ho_se_ind", "He_se_ind", "Fis_se_ind") %in% names(res$per_population)))
  expect_true(all(c("Ar_se_ind", "privAr_se_ind") %in% names(res$richness)))
  ## each _se_ind sits right after its locus-based _se
  nm <- names(res$per_population)
  expect_equal(match("Ho_se_ind", nm), match("Ho_se", nm) + 1L)
})

test_that("individual SEs grow relative to locus SEs when individuals differ in inbreeding", {
  ratio <- function(d) { pp <- run_ind(d)$per_population; pp$Ho_se_ind[1] / pp$Ho_se[1] }
  same   <- sim_two_pops(rep(0.1, 20), 1000, seed = 2)
  varied <- sim_two_pops(rep(c(0, 0.6), 10), 1000, seed = 3)
  expect_gt(ratio(varied), 2 * ratio(same))
})

test_that("Ar/privAr individual SEs are NA when removing one individual leaves fewer than g copies", {
  d <- sim_two_pops(rep(0.1, 8), 200, seed = 4)      # 4 per population: 2(n-1) = 6
  expect_message(res <- diversity_stats(d$H, d$popmap, g = 8, nboot = 0, stem = "j",
                                        se_individuals = TRUE), "g <= 6")
  expect_true(all(is.na(res$richness$Ar_se_ind)))
  expect_false(anyNA(res$per_population$He_se_ind))
})

test_that("the report explains the _se_ind columns only when they are there", {
  d <- sim_two_pops(rep(0.1, 12), 200, seed = 5)
  expect_true(any(grepl("_se_ind columns", capture.output(print(run_ind(d))))))
  plain <- suppressMessages(diversity_stats(d$H, d$popmap, g = 4, nboot = 0, stem = "j"))
  expect_false(any(grepl("_se_ind", capture.output(print(plain)))))
  expect_false("Ho_se_ind" %in% names(plain$per_population))
})
