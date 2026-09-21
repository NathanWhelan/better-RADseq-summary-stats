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
  expect_equal(res$per_population$Ho_se_ind[1], sd(hA) / sqrt(6))
  expect_true(all(c("Ho_se_ind", "He_se_ind", "Fis_se_ind") %in% names(res$per_population)))
  ## Ar and privAr have no individual SE: their locus SE covers both sources.
  expect_false(any(c("Ar_se_ind", "privAr_se_ind") %in% names(res$richness)))
  ## _se, _se_ind and _se_combined sit side by side, and combined adds variances.
  nm <- names(res$per_population)
  expect_equal(match(c("Ho_se_ind", "Ho_se_combined"), nm), match("Ho_se", nm) + 1:2)
  pp <- res$per_population
  expect_equal(pp$He_se_combined, sqrt(pp$He_se^2 + pp$He_se_ind^2))
  expect_equal(pp$Fis_se_combined, sqrt(pp$Fis_se^2 + pp$Fis_se_ind^2))
})

test_that("individual SEs grow relative to locus SEs when individuals differ in inbreeding", {
  ratio <- function(d) { pp <- run_ind(d)$per_population; pp$Ho_se_ind[1] / pp$Ho_se[1] }
  same   <- sim_two_pops(rep(0.1, 20), 1000, seed = 2)
  varied <- sim_two_pops(rep(c(0, 0.6), 10), 1000, seed = 3)
  expect_gt(ratio(varied), 2 * ratio(same))
})

test_that("the individual jackknife uses the same records in every replicate", {
  ## Keep only records where both populations have >= 5 of their 10
  ## individuals genotyped. min_n = 5 and min_n = 2 then choose the same
  ## records, so with the record set fixed in the jackknife every result must
  ## be identical. (Re-deciding records per replicate would drop the records
  ## with exactly 5 genotyped individuals from some replicates under min_n = 5.)
  d <- sim_two_pops(rep(0.1, 20), 600, seed = 7)
  set.seed(7)
  miss <- matrix(stats::runif(length(d$H$A1)) < 0.35, nrow(d$H$A1))
  d$H$A1[miss] <- NA; d$H$A2[miss] <- NA
  typed <- cbind(rowSums(!is.na(d$H$A1[, 1:10])), rowSums(!is.na(d$H$A1[, 11:20])))
  keep <- rowSums(typed >= 5) == 2
  expect_true(any(typed[keep, ] == 5))                  # the case being tested exists
  H <- RADdiversity:::.subset_H(d$H, keep)
  run <- function(min_n) suppressMessages(diversity_stats(H, d$popmap, g = 4, nboot = 0, stem = "j",
                                                          min_n = min_n, se_individuals = TRUE))
  expect_equal(run(5)$per_population, run(2)$per_population)
})

test_that("a warning flags individual SEs inflated by records with only 2 genotyped individuals", {
  d <- sim_two_pops(rep(0.1, 12), 500, seed = 6)       # 6 per population
  set.seed(6)
  miss <- matrix(stats::runif(length(d$H$A1)) < 0.6, nrow(d$H$A1))
  d$H$A1[miss] <- NA; d$H$A2[miss] <- NA
  expect_warning(res <- suppressMessages(diversity_stats(d$H, d$popmap, g = 4, nboot = 0,
                                                         stem = "j", se_individuals = TRUE)),
                 "probably too large")
  expect_true(all(res$settings$jackknife_boundary > RADdiversity:::.jackknife_boundary_share))
  expect_true(any(grepl("probably too large", capture.output(print(res)))))
  ## min_n = 3: no used record has only 2 genotyped individuals.
  expect_no_warning(res3 <- suppressMessages(
    diversity_stats(d$H, d$popmap, g = 4, nboot = 0, stem = "j", min_n = 3, se_individuals = TRUE)))
  expect_true(all(res3$settings$jackknife_boundary == 0))
})

test_that("the report explains the _se_ind columns only when they are there", {
  d <- sim_two_pops(rep(0.1, 12), 200, seed = 5)
  expect_true(any(grepl("_se_ind columns", capture.output(print(summary(run_ind(d), details = TRUE))))))
  plain <- suppressMessages(diversity_stats(d$H, d$popmap, g = 4, nboot = 0, stem = "j",
                                            se_individuals = FALSE))
  expect_false(any(grepl("_se_ind", capture.output(print(summary(plain, details = TRUE))))))
  expect_false(any(grepl("_se_ind|_se_combined", capture.output(print(summary(plain))))))
  expect_false("Ho_se_ind" %in% names(plain$per_population))
})

test_that("the individual SEs are on by default, and se_individuals = FALSE turns them off", {
  d <- sim_two_pops(rep(0.1, 12), 200, seed = 5)
  on <- suppressMessages(diversity_stats(d$H, d$popmap, g = 4, nboot = 0, stem = "j"))
  expect_true(on$settings$se_individuals)
  expect_true(all(c("Ho_se_ind", "Ho_se_combined", "He_se_combined", "Fis_se_combined") %in%
                    names(on$per_population)))
  off <- suppressMessages(diversity_stats(d$H, d$popmap, g = 4, nboot = 0, stem = "j",
                                          se_individuals = FALSE))
  expect_false(off$settings$se_individuals)
  expect_false(any(grepl("_se_ind|_se_combined", names(off$per_population))))
  ## Turning them off brings back the note that the SEs hold individuals fixed.
  expect_true(any(grepl("se_individuals = TRUE", capture.output(print(off)))))
  expect_false(any(grepl("se_individuals = TRUE", capture.output(print(on)))))
})
