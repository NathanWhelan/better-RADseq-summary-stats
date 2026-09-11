## Simulated genotypes: individual i is inbred with probability F[i] at each
## locus (its two alleles identical by descent), so a spread of F across
## individuals creates identity disequilibrium and a constant F does not.
sim_F <- function(F, L, seed) {
  set.seed(seed)
  n <- length(F)
  p <- stats::runif(L, 0.1, 0.9)
  a1 <- matrix(1L + (stats::runif(L * n) < p), L, n)
  ibd <- matrix(stats::runif(L * n) < rep(F, each = L), L, n)
  a2 <- ifelse(ibd, a1, 1L + (stats::runif(L * n) < p))
  samp <- paste0("s", seq_len(n))
  dimnames(a1) <- dimnames(a2) <- list(NULL, samp)
  list(A1 = a1, A2 = a2, locus = paste0("l", seq_len(L)), locus_raw = paste0("l", seq_len(L)),
       alleles = replicate(L, c("A", "C"), simplify = FALSE), n_alleles = rep(2L, L),
       samples = samp)
}

## Brute force straight from the definition: mean of h_il h_il' over pairs of
## distinct loci typed in the same individual, divided by the mean of
## h_il h_jl' over pairs of distinct loci typed in two different individuals.
g2_brute <- function(a1, a2) {
  m <- !is.na(a1); h <- m & (a1 != a2)
  L <- nrow(h); n <- ncol(h)
  s_num <- s_den <- d_num <- d_den <- 0
  for (i in 1:n) for (j in 1:n) for (l in 1:L) for (k in 1:L) {
    if (l == k || !m[l, i] || !m[k, j]) next
    v <- h[l, i] * h[k, j]
    if (i == j) { s_num <- s_num + v; s_den <- s_den + 1 }
    else        { d_num <- d_num + v; d_den <- d_den + 1 }
  }
  (s_num / s_den) / (d_num / d_den) - 1
}

test_that("g2 matches a brute-force evaluation of its definition, with missing data", {
  H <- sim_F(runif(6, 0, 0.6), 8, seed = 1)
  H$A1[cbind(c(1, 3, 5, 8), c(2, 2, 4, 6))] <- NA
  H$A2[cbind(c(1, 3, 5, 8), c(2, 2, 4, 6))] <- NA
  got <- RADdiversity:::.g2_summary(H$A1, H$A2, nboot = 0, nperm = 0)$g2
  expect_equal(got, unname(g2_brute(H$A1, H$A2)), tolerance = 1e-12)
})

test_that("with complete data g2 is the closed form (n - 1) P_same / P_diff - 1", {
  H <- sim_F(rep(0.1, 10), 50, seed = 2)
  h <- (H$A1 != H$A2) * 1
  r <- colSums(h); cl <- rowSums(h); n <- ncol(h)
  Ps <- sum(r^2 - r); Pd <- sum(cl)^2 - sum(cl^2) - Ps
  expect_equal(RADdiversity:::.g2_summary(H$A1, H$A2, 0, 0)$g2, (n - 1) * Ps / Pd - 1,
               tolerance = 1e-12)
})

test_that("g2 is ~0 when individuals share one F, and > 0 when F varies among them", {
  same <- sim_F(rep(0.1, 30), 2000, seed = 3)
  varied <- sim_F(rep(c(0, 0.5), 15), 2000, seed = 4)
  res <- identity_disequilibrium(same, list(pop = same$samples), nboot = 200, nperm = 200)
  expect_true(res$g2_lo < 0 && res$g2_hi > 0)            # CI covers 0
  res2 <- identity_disequilibrium(varied, list(pop = varied$samples), nboot = 200, nperm = 200)
  expect_gt(res2$g2_lo, 0)
  expect_lt(res2$p_value, 0.01)
})

test_that("identity_disequilibrium() restores the caller's RNG state", {
  H <- sim_F(rep(0.1, 8), 100, seed = 5)
  set.seed(9); before <- .Random.seed
  identity_disequilibrium(H, list(pop = H$samples), nboot = 20, nperm = 20)
  expect_identical(.Random.seed, before)
})

test_that("with complete data g2 agrees with inbreedR::g2_snps()", {
  skip_if_not_installed("inbreedR")
  H <- sim_F(rep(c(0, 0.3), 10), 300, seed = 6)
  geno <- t((H$A1 != H$A2) * 1)                   # individuals x loci, 1 = heterozygous
  ref <- inbreedR::g2_snps(geno, nperm = 0, nboot = 0, verbose = FALSE)$g2
  expect_equal(RADdiversity:::.g2_summary(H$A1, H$A2, 0, 0)$g2, ref, tolerance = 1e-8)
})
