###############################################################################
#
#  R/selftest.R -- self-tests for the estimators (diversity_core_selftest())
#  and for the individual-level test's type I error calibration
#  (het_between_pops_selftest()). Both need no data files and no packages.
#
###############################################################################

#' Check the estimators in R/estimators.R against brute-force Monte Carlo
#'
#' Reproduces the checks run by `Rscript diversity_core.R --selftest`: every
#' estimator here against brute-force Monte Carlo, and the Stacks formulas
#' against real published Stacks output. Prints a report to the console.
#'
#' @return Invisibly, `TRUE` if every check passed, `FALSE` otherwise. The
#'   caller's own RNG state is restored when this function returns, so
#'   calling it interactively does not affect subsequent random draws in
#'   the caller's session.
#' @examples
#' \donttest{
#' diversity_core_selftest()   # about 20 seconds of Monte Carlo
#' }
#' @export
diversity_core_selftest <- function() {

  ## Restore the caller's RNG state on exit -- see the matching comment in
  ## diversity_stats(). This function calls set.seed() repeatedly for its
  ## own Monte Carlo checks below.
  restore_rng <- .save_rng_state()
  on.exit(restore_rng(), add = TRUE)

  ok <- TRUE
  chk <- function(label, got, want, tol = 1e-9) {
    good <- isTRUE(all(abs(got - want) < tol))
    ok <<- ok && good
    cat(sprintf("  %-58s %-12s %s\n", label,
                paste(signif(got, 6), collapse = " "),
                if (good) "PASS" else sprintf("FAIL (expected %s)",
                                              paste(signif(want, 6), collapse = " "))))
  }

  cat("\n--- 1. Nei-Chesser Hs is UNBIASED at any FIS; Stacks Pi is not ---\n")
  set.seed(11)
  simulate_bias <- function(n, F0, p0, reps = 40000) {
    ## draw n diploids with inbreeding coefficient F0 at a biallelic locus
    ib  <- stats::runif(reps * n) < F0
    a1  <- stats::runif(reps * n) < p0
    a2  <- ifelse(ib, a1, stats::runif(reps * n) < p0)
    A1 <- matrix(a1, reps, n); A2 <- matrix(a2, reps, n)
    p  <- (rowSums(A1) + rowSums(A2)) / (2 * n)
    ho <- rowMeans(A1 != A2)
    c(true = 2 * p0 * (1 - p0),
      nei_chesser = mean(hs_biallelic(p, ho, n)),
      stacks_pi   = mean(hs_stacks_pi(p, n)))
  }
  for (cfg in list(c(n = 10, F = 0.00), c(n = 10, F = 0.20),
                   c(n = 15, F = 0.40), c(n = 30, F = 0.20))) {
    v <- simulate_bias(cfg[1], cfg[2], 0.35)
    cat(sprintf("  n=%2d F=%.2f  true %.5f | Nei-Chesser %.5f (%+.2f%%) | Stacks Pi %.5f (%+.2f%%)\n",
                cfg[1], cfg[2], v["true"], v["nei_chesser"],
                100 * (v["nei_chesser"] / v["true"] - 1), v["stacks_pi"],
                100 * (v["stacks_pi"] / v["true"] - 1)))
  }
  cat("  Monte-Carlo error is about +/-0.15%. Nei-Chesser sits on zero at every\n")
  cat("  F; Stacks Pi is low by roughly F/(2n-1), growing as n shrinks.\n")

  cat("\n--- 2. FIS reaches its theoretical bounds ---\n")
  n <- 10
  ## every individual heterozygous  -> true FIS = -1
  chk("all heterozygous, Nei-Chesser FIS (must be exactly -1)",
      1 - 1 / hs_biallelic(0.5, 1, n), -1)
  chk("all heterozygous, Stacks-Pi FIS (wrong: -(n-1)/n)",
      1 - 1 / hs_stacks_pi(0.5, n), -(n - 1) / n)
  ## no heterozygotes at p = 0.5    -> FIS = +1
  chk("no heterozygotes at p=0.5, Nei-Chesser FIS", 1 - 0 / hs_biallelic(0.5, 0, n), 1)

  cat("\n--- 3. single minor-allele copy forces FIS = 0 (both estimators) ---\n")
  chk("MAC=1: Nei-Chesser Hs equals Ho exactly", hs_biallelic(1 / (2 * n), 1 / n, n), 1 / n)

  cat("\n--- 4. ratio of sums is invariant to monomorphic sites ---\n")
  set.seed(2); L <- 500
  pv <- stats::runif(L, 0.05, 0.95); nn <- rep(12, L)
  hov <- pmax(0, pmin(1, 2 * pv * (1 - pv) * 0.93))
  hev <- hs_biallelic(pv, hov, nn)
  f_poly <- fis_ratio_of_sums(hov, hev)
  f_pad  <- fis_ratio_of_sums(c(hov, rep(0, 800)), c(hev, rep(0, 800)))
  chk("FIS with vs without 800 monomorphic sites appended", f_pad, f_poly)

  cat("\n--- 5. rarefied allelic richness vs brute-force Monte Carlo ---\n")
  set.seed(5); cts <- c(20, 12, 6, 3, 2, 1, 1); pool <- rep(seq_along(cts), cts)
  for (g in c(5, 10, 20)) {
    an <- rare_richness(cts, g)
    mc <- mean(replicate(60000, length(unique(sample(pool, g)))))
    cat(sprintf("  g=%2d  analytic %.4f  simulated %.4f  diff %.4f (MC err ~0.004)\n",
                g, an, mc, abs(an - mc)))
    ok <- ok && abs(an - mc) < 0.02
  }
  chk("at g = N, rarefied richness equals observed richness",
      rare_richness(c(7, 3, 2, 1), 13), 4)

  cat("\n--- 6. rarefied PRIVATE allelic richness vs brute-force Monte Carlo ---\n")
  set.seed(6)
  m <- rbind(c(10, 8, 4, 2, 0, 0), c(12, 0, 0, 3, 5, 4))
  pa <- rep(seq_len(ncol(m)), m[1, ]); pb <- rep(seq_len(ncol(m)), m[2, ])
  for (g in c(4, 8, 12)) {
    mc <- rowMeans(replicate(40000, {
      sa <- unique(sample(pa, g)); sb <- unique(sample(pb, g))
      c(length(setdiff(sa, sb)), length(setdiff(sb, sa))) }))
    an <- rare_private_all(m, g)
    cat(sprintf("  g=%2d  pop1 analytic %.4f sim %.4f | pop2 analytic %.4f sim %.4f\n",
                g, an[1], mc[1], an[2], mc[2]))
    ok <- ok && max(abs(an - mc)) < 0.03
  }

  cat("\n--- 7. three-population private-allele logic, hand-checked ---\n")
  ## allele 1 only in pop1; allele 2 in all; allele 3 only in pop3
  m3 <- rbind(c(4, 4, 0), c(0, 8, 0), c(0, 4, 4))
  pr <- rare_private_all(m3, 8)   # g = N, so exact presence/absence
  chk("private counts at g = N (expect 1, 0, 1)", pr, c(1, 0, 1))

  cat("\n--- 8. autosomal conversion (Schmidt et al. 2021) ---\n")
  chk("He 0.2911 over 10,348 SNPs of 1,000,000 sequenced sites",
      autosomal_het(0.2911, 10348, 1e6), 0.2911 * 10348 / 1e6, tol = 1e-12)

  cat("\n--- 9. the two Stacks-Pi forms agree, and Pi = Exp_Het * 2n/(2n-1) ---\n")
  for (n in c(2L, 10L, 15L)) for (p in c(0.1, 0.35, 0.5)) {
    cnt <- c(round(2 * n * p), 2 * n - round(2 * n * p))
    ph <- cnt[1] / sum(cnt)
    chk(sprintf("n=%2d p=%.2f  counts form == biallelic form", n, ph),
        gene_div_2n_counts(cnt), hs_stacks_pi(ph, n), tol = 1e-12)
  }
  ## The relationship a user can check on their own populations.sumstats.tsv,
  ## anchored on a real Stacks summary file: Exp_Het 0.40000 at Num_Indv 2
  ## is reported as Pi 0.53333.
  chk("published Stacks row: Exp_Het 0.4 at n=2 -> Pi",
      0.40000 * (2 * 2) / (2 * 2 - 1), 0.53333, tol = 1e-5)
  cat("  Pi / Exp_Het ratio:")
  for (n in c(2L, 10L, 15L, 25L))
    cat(sprintf("  n=%d %.3f", n, (2 * n) / (2 * n - 1)))
  cat("\n  So compare this package's He against Stacks' Pi, never Exp_Het.\n")

  cat("\n--- 10. Stacks formulas, against real published Stacks output ---\n")
  ## Rows copied verbatim from populations.sumstats.tsv in the Galaxy stacks2
  ## test data (N, P, Obs Het, Exp Het, Pi, Fis). Catchen, on stacks-users:
  ## "The calculation for Fis uses pi for expected heterozygosity:
  ##  Fis = (pi - obs_het)/pi. Pi is calculated as
  ##  pi = 1 - Sum_i( (n_i choose 2) ) / (n choose 2)."  (Hohenlohe et al. 2010)
  st <- data.frame(N = c(2, 2), P = c(0.50, 0.75), oh = c(1.0, 0.5),
                   eh = c(0.50000, 0.37500), pi = c(0.66667, 0.50000),
                   fis = c(-0.50000, 0.00000))
  for (i in seq_len(nrow(st))) {
    row <- st[i, ]; n <- 2 * row$N; sp2 <- row$P^2 + (1 - row$P)^2
    ni <- c(round(n * row$P), n - round(n * row$P))
    chk(sprintf("N=%d p=%.2f  Exp Het = 1 - sum p^2", row$N, row$P), 1 - sp2, row$eh, 1e-4)
    chk("            Pi = (n/(n-1))(1 - sum p^2)", (n / (n - 1)) * (1 - sp2), row$pi, 1e-4)
    chk("            Pi = 1 - sum C(n_i,2)/C(n,2)",
        1 - sum(choose(ni, 2)) / choose(n, 2), row$pi, 1e-4)
    chk("            Fis = (Pi - Obs Het)/Pi", (row$pi - row$oh) / row$pi, row$fis, 1e-4)
  }
  chk("summary row: Exp_Het 0.4 at N=2 -> Pi", 0.4 * 4 / 3, 0.53333, 1e-5)
  cat("  So: Stacks divides FIS by Pi (which assumes FIS = 0), and averages\n")
  cat("  per-locus ratios. Two independent errors. See README.md.\n")

  cat("\n--- 11. FIS ratio-of-sums stays unbiased under LOCUS-DRIVEN missingness ---\n")
  ## Real RAD data drops loci in clusters (a locus fails broadly) rather than
  ## dropping individuals independently at random. Simulate that: each locus
  ## has probability q of being "degraded", in which case a random subset of
  ## individuals is missing there, so per-locus n varies a lot from locus to
  ## locus even though every individual's own overall call rate stays high.
  ## fis_ratio_of_sums() should stay close to the true FIS regardless, because
  ## it weights each locus by its own information (via the sums), unlike a
  ## naive mean of per-locus ratios.
  set.seed(31)
  n <- 20; Ltest <- 3000; F0 <- 0.15; p0 <- 0.35; q <- 0.30
  sim_available <- function() {
    degraded <- stats::runif(Ltest) < q
    ho <- he <- rep(NA_real_, Ltest)
    for (j in seq_len(Ltest)) {
      nj <- if (degraded[j]) sample(2:(n - 1), 1) else n
      ib <- stats::runif(nj) < F0
      a1 <- stats::runif(nj) < p0
      a2 <- ifelse(ib, a1, stats::runif(nj) < p0)
      p  <- (sum(a1) + sum(a2)) / (2 * nj)
      hoj <- mean(a1 != a2)
      ho[j] <- hoj; he[j] <- hs_biallelic(p, hoj, nj)
    }
    c(ros = fis_ratio_of_sums(ho, he),
      mor = mean(1 - ho / he, na.rm = TRUE))     # naive mean of per-locus ratios
  }
  v <- rowMeans(replicate(300, sim_available()))
  cat(sprintf("  true FIS = %.3f | ratio-of-sums %.4f (%+.1f%% off) | mean-of-ratios %.4f (%+.1f%% off)\n",
              F0, v["ros"], 100 * (v["ros"] / F0 - 1), v["mor"], 100 * (v["mor"] / F0 - 1)))
  chk("ratio-of-sums FIS within 5% of true FIS under locus-driven missingness",
      abs(v["ros"] / F0 - 1) < 0.05, TRUE)
  cat("  (mean-of-ratios is shown for contrast, not required to pass -- it's the\n")
  cat("  aggregation this package deliberately does not use.)\n")

  cat("\n", if (ok) "ALL CHECKS PASSED\n" else "*** SOME CHECKS FAILED ***\n", sep = "")
  invisible(ok)
}

#' Calibration self-test for het_between_pops()
#'
#' Reproduces the type I error table run by `Rscript het_between_pops.R
#' --selftest`: 400 simulations of two populations with IDENTICAL true
#' heterozygosity, comparing the rejection rate of Welch's t / Wilcoxon
#' (individual-level) against a locus bootstrap, for a nominal 5% test.
#' Prints a report to the console.
#'
#' @return Invisibly, `TRUE`. The caller's own RNG state is restored when
#'   this function returns, so calling it interactively does not affect
#'   subsequent random draws in the caller's session.
#' @examples
#' \donttest{
#' het_between_pops_selftest()   # a few seconds of simulation
#' }
#' @export
het_between_pops_selftest <- function() {
  ## Restore the caller's RNG state on exit -- see the matching comment in
  ## diversity_stats(). This function calls set.seed() repeatedly for its
  ## own Monte Carlo checks below.
  restore_rng <- .save_rng_state()
  on.exit(restore_rng(), add = TRUE)

  cat("--- calibration of the individual-level tests, 400 simulations ---\n")
  cat("    two populations with IDENTICAL true Ho; target rejection rate 5%\n\n")
  set.seed(404); nsim <- 400; L <- 1500; n1 <- 15; n2 <- 10
  sim <- function(f_sd) {
    p <- stats::runif(L, 0.05, 0.95)
    F1 <- pmin(pmax(stats::rnorm(n1, .1, f_sd), 0), .9)
    F2 <- pmin(pmax(stats::rnorm(n2, .1, f_sd), 0), .9)
    ih <- function(Fv, n) vapply(seq_len(n), function(j)
            mean(stats::runif(L) < 2*p*(1-p)*(1-Fv[j])), numeric(1))
    list(a = ih(F1, n1), b = ih(F2, n2))
  }
  cat(sprintf("  %-24s %10s %10s\n", "individual F variation", "Welch t", "Wilcoxon"))
  for (sd in c(0.02, 0.10, 0.25)) {
    rr <- replicate(nsim, { x <- sim(sd)
      c(stats::t.test(x$a, x$b)$p.value,
        suppressWarnings(stats::wilcox.test(x$a, x$b))$p.value) })
    cat(sprintf("  SD(F) = %.2f %-13s %9.1f%% %9.1f%%\n", sd, "",
                100 * mean(rr[1, ] < 0.05), 100 * mean(rr[2, ] < 0.05)))
  }
  cat("\n--- and what the LOCUS bootstrap does on the same data, for contrast ---\n")
  set.seed(202); L <- 1500
  for (sd in c(0.00, 0.10)) {
    rej <- replicate(120, {
      p <- stats::runif(L, .05, .95)
      F1 <- pmin(pmax(stats::rnorm(n1, .1, sd), 0), .9)
      F2 <- pmin(pmax(stats::rnorm(n2, .1, sd), 0), .9)
      mk <- function(Fv, n) { h <- matrix(FALSE, L, n)
        for (j in seq_len(n)) h[, j] <- stats::runif(L) < 2*p*(1-p)*(1-Fv[j]); h }
      h1 <- mk(F1, n1); h2 <- mk(F2, n2)
      a <- rowMeans(h1); b <- rowMeans(h2)
      bo <- replicate(300, { i <- sample.int(L, L, TRUE); mean(a[i]) - mean(b[i]) })
      !(0 >= stats::quantile(bo, .025) && 0 <= stats::quantile(bo, .975)) })
    cat(sprintf("  SD(F) = %.2f -> locus bootstrap rejects %.0f%% of the time\n",
                sd, 100 * mean(rej)))
  }
  cat("\n  The individual-level tests hold their nominal rate; the locus bootstrap\n")
  cat("  does not, as soon as individuals differ from one another.\n")
  invisible(TRUE)
}
