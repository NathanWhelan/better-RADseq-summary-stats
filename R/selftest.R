###############################################################################
#
#  R/selftest.R -- self-tests a user (or a reviewer) can run on their own
#  computer without any data files or extra packages:
#
#    diversity_core_selftest()    every estimator in R/estimators.R against
#                                 brute-force Monte Carlo and published
#                                 Stacks output
#    het_between_pops_selftest()  the type I error of het_between_pops()'s
#                                 combined test, of Welch's t and Wilcoxon
#                                 alone, and of a locus bootstrap
#
#  Both use fixed seeds (so they give the same numbers every time) and put
#  the caller's own random-number state back when they finish. The same
#  checks, and many more, run in tests/testthat/.
#
###############################################################################

## Not exported. Collects self-test results. `check(label, got, want, tol)`
## records one row (and prints it when verbose); `results()` returns them all
## as a data frame. `rows` lives in this function, and `<<-` inside check()
## adds to that shared `rows` rather than creating a new local copy.
.selftest_recorder <- function(verbose) {
  rows <- list()
  check <- function(label, got, want, tol = 1e-9) {
    pass <- isTRUE(all(abs(got - want) < tol))
    rows[[length(rows) + 1L]] <<- data.frame(
      check = label, value = paste(signif(got, 6), collapse = " "),
      expected = paste(signif(want, 6), collapse = " "), tolerance = tol, pass = pass)
    if (verbose)
      cat(sprintf("  %-58s %-12s %s\n", label, paste(signif(got, 6), collapse = " "),
                  if (pass) "PASS" else sprintf("FAIL (expected %s)",
                                                paste(signif(want, 6), collapse = " "))))
    invisible(pass)
  }
  results <- function() do.call(rbind, rows)
  list(check = check, results = results)
}

#' Check the estimators against brute-force Monte Carlo
#'
#' Runs the checks behind `Rscript diversity_core.R --selftest`: every
#' estimator in [estimators] against brute-force Monte Carlo, and the Stacks
#' formulas against published Stacks output.
#'
#' @param verbose Print the report. Default `TRUE`.
#' @return Invisibly, a data frame with one row per check: `check`, `value`,
#'   `expected`, `tolerance` and `pass`. `all(result$pass)` is `TRUE` when
#'   every check passed. Your own random-number state is left as it was.
#' @examples
#' \donttest{
#' result <- diversity_core_selftest()   # about 20 seconds of Monte Carlo
#' all(result$pass)
#' }
#' @export
diversity_core_selftest <- function(verbose = TRUE) {
  .check_flag(verbose, "verbose")
  restore_rng <- .save_rng_state()
  on.exit(restore_rng(), add = TRUE)
  rec <- .selftest_recorder(verbose)
  check <- rec$check
  say <- function(...) if (verbose) cat(...)

  say("\n--- 1. Nei-Chesser Hs is UNBIASED at any FIS; Stacks Pi is not ---\n")
  set.seed(11)
  simulate_bias <- function(n, F0, p0, reps = 40000) {
    ## n diploids with inbreeding coefficient F0 at a biallelic locus
    ibd <- stats::runif(reps * n) < F0
    a1 <- stats::runif(reps * n) < p0
    a2 <- ifelse(ibd, a1, stats::runif(reps * n) < p0)
    A1 <- matrix(a1, reps, n)
    A2 <- matrix(a2, reps, n)
    p <- (rowSums(A1) + rowSums(A2)) / (2 * n)
    ho <- rowMeans(A1 != A2)
    c(true = 2 * p0 * (1 - p0),
      nei_chesser = mean(hs_biallelic(p, ho, n)),
      stacks_pi = mean(hs_stacks_pi(p, n)))
  }
  for (setting in list(c(n = 10, F = 0.00), c(n = 10, F = 0.20),
                       c(n = 15, F = 0.40), c(n = 30, F = 0.20))) {
    v <- simulate_bias(setting[["n"]], setting[["F"]], 0.35)
    say(sprintf("  n=%2d F=%.2f  true %.5f | Nei-Chesser %.5f (%+.2f%%) | Stacks Pi %.5f (%+.2f%%)\n",
                setting[["n"]], setting[["F"]], v[["true"]], v[["nei_chesser"]],
                100 * (v[["nei_chesser"]] / v[["true"]] - 1), v[["stacks_pi"]],
                100 * (v[["stacks_pi"]] / v[["true"]] - 1)))
  }
  say("  Monte-Carlo error is about +/-0.15%. Nei-Chesser sits on zero at every\n")
  say("  F; Stacks Pi is low by roughly F/(2n-1), growing as n shrinks.\n")

  say("\n--- 2. FIS reaches its theoretical bounds ---\n")
  n <- 10
  check("all heterozygous, Nei-Chesser FIS (must be exactly -1)",
        1 - 1 / hs_biallelic(0.5, 1, n), -1)
  check("all heterozygous, Stacks-Pi FIS (wrong: -(n-1)/n)",
        1 - 1 / hs_stacks_pi(0.5, n), -(n - 1) / n)
  check("no heterozygotes at p=0.5, Nei-Chesser FIS", 1 - 0 / hs_biallelic(0.5, 0, n), 1)

  say("\n--- 3. single minor-allele copy forces FIS = 0 (both estimators) ---\n")
  check("MAC=1: Nei-Chesser Hs equals Ho exactly", hs_biallelic(1 / (2 * n), 1 / n, n), 1 / n)

  say("\n--- 4. ratio of sums is invariant to monomorphic sites ---\n")
  set.seed(2)
  n_loci <- 500
  p_values <- stats::runif(n_loci, 0.05, 0.95)
  ho_values <- pmax(0, pmin(1, 2 * p_values * (1 - p_values) * 0.93))
  he_values <- hs_biallelic(p_values, ho_values, rep(12, n_loci))
  check("FIS with vs without 800 monomorphic sites appended",
        fis_ratio_of_sums(c(ho_values, rep(0, 800)), c(he_values, rep(0, 800))),
        fis_ratio_of_sums(ho_values, he_values))

  say("\n--- 5. rarefied allelic richness vs brute-force Monte Carlo ---\n")
  set.seed(5)
  counts <- c(20, 12, 6, 3, 2, 1, 1)
  gene_pool <- rep(seq_along(counts), counts)
  for (g in c(5, 10, 20)) {
    simulated <- mean(replicate(60000, length(unique(sample(gene_pool, g)))))
    check(sprintf("g=%2d analytic vs simulated (MC error ~0.004)", g),
          rare_richness(counts, g), simulated, tol = 0.02)
  }
  check("at g = N, rarefied richness equals observed richness",
        rare_richness(c(7, 3, 2, 1), 13), 4)

  say("\n--- 6. rarefied PRIVATE allelic richness vs brute-force Monte Carlo ---\n")
  set.seed(6)
  count_mat <- rbind(c(10, 8, 4, 2, 0, 0), c(12, 0, 0, 3, 5, 4))
  pool_a <- rep(seq_len(ncol(count_mat)), count_mat[1, ])
  pool_b <- rep(seq_len(ncol(count_mat)), count_mat[2, ])
  for (g in c(4, 8, 12)) {
    simulated <- rowMeans(replicate(40000, {
      drawn_a <- unique(sample(pool_a, g))
      drawn_b <- unique(sample(pool_b, g))
      c(length(setdiff(drawn_a, drawn_b)), length(setdiff(drawn_b, drawn_a)))
    }))
    check(sprintf("g=%2d both populations, analytic vs simulated", g),
          rare_private_all(count_mat, g), simulated, tol = 0.03)
  }

  say("\n--- 7. three-population private-allele logic, hand-checked ---\n")
  ## allele 1 only in pop1; allele 2 in all; allele 3 only in pop3
  three_pops <- rbind(c(4, 4, 0), c(0, 8, 0), c(0, 4, 4))
  check("private counts at g = N (expect 1, 0, 1)", rare_private_all(three_pops, 8), c(1, 0, 1))

  say("\n--- 8. autosomal conversion (Schmidt et al. 2021) ---\n")
  check("He 0.2911 over 10,348 SNPs of 1,000,000 sequenced sites",
        autosomal_het(0.2911, 10348, 1e6), 0.2911 * 10348 / 1e6, tol = 1e-12)

  say("\n--- 9. the two Stacks-Pi forms agree, and Pi = Exp_Het * 2n/(2n-1) ---\n")
  for (n in c(2L, 10L, 15L)) {
    for (p in c(0.1, 0.35, 0.5)) {
      allele_counts <- c(round(2 * n * p), 2 * n - round(2 * n * p))
      p_hat <- allele_counts[1] / sum(allele_counts)
      check(sprintf("n=%2d p=%.2f  counts form == biallelic form", n, p_hat),
            gene_div_2n_counts(allele_counts), hs_stacks_pi(p_hat, n), tol = 1e-12)
    }
  }
  ## Anchored on a real Stacks summary file: Exp_Het 0.40000 at Num_Indv 2 is
  ## reported as Pi 0.53333.
  check("published Stacks row: Exp_Het 0.4 at n=2 -> Pi", 0.40000 * (2 * 2) / (2 * 2 - 1),
        0.53333, tol = 1e-5)
  say("  Pi / Exp_Het ratio:")
  for (n in c(2L, 10L, 15L, 25L)) say(sprintf("  n=%d %.3f", n, (2 * n) / (2 * n - 1)))
  say("\n  So compare this package's He against Stacks' Pi, never Exp_Het.\n")

  say("\n--- 10. Stacks formulas, against real published Stacks output ---\n")
  ## Rows copied from populations.sumstats.tsv in the Galaxy stacks2 test data
  ## (N, P, Obs Het, Exp Het, Pi, Fis). Catchen, on stacks-users: "The
  ## calculation for Fis uses pi for expected heterozygosity:
  ## Fis = (pi - obs_het)/pi. Pi is calculated as
  ## pi = 1 - Sum_i( (n_i choose 2) ) / (n choose 2)." (Hohenlohe et al. 2010)
  stacks_rows <- data.frame(N = c(2, 2), P = c(0.50, 0.75), obs_het = c(1.0, 0.5),
                            exp_het = c(0.50000, 0.37500), pi = c(0.66667, 0.50000),
                            fis = c(-0.50000, 0.00000))
  for (i in seq_len(nrow(stacks_rows))) {
    row <- stacks_rows[i, ]
    n_copies <- 2 * row$N
    sum_p2 <- row$P^2 + (1 - row$P)^2
    copies <- c(round(n_copies * row$P), n_copies - round(n_copies * row$P))
    check(sprintf("N=%d p=%.2f  Exp Het = 1 - sum p^2", row$N, row$P), 1 - sum_p2, row$exp_het, 1e-4)
    check("            Pi = (n/(n-1))(1 - sum p^2)", (n_copies / (n_copies - 1)) * (1 - sum_p2),
          row$pi, 1e-4)
    check("            Pi = 1 - sum C(n_i,2)/C(n,2)",
          1 - sum(choose(copies, 2)) / choose(n_copies, 2), row$pi, 1e-4)
    check("            Fis = (Pi - Obs Het)/Pi", (row$pi - row$obs_het) / row$pi, row$fis, 1e-4)
  }
  check("summary row: Exp_Het 0.4 at N=2 -> Pi", 0.4 * 4 / 3, 0.53333, 1e-5)
  say("  So: Stacks divides FIS by Pi (which assumes FIS = 0), and averages\n")
  say("  per-locus ratios; see vignette(\"rationale\"), section 2.\n")

  say("\n--- 11. FIS ratio-of-sums stays unbiased under LOCUS-DRIVEN missingness ---\n")
  ## Real RAD data lose loci in clusters (a locus fails broadly), not
  ## individuals independently at random. Each simulated locus is "degraded"
  ## with probability q, and then a random subset of individuals is missing
  ## there, so per-locus n varies a lot. The ratio of sums should stay close
  ## to the true FIS; a mean of per-locus ratios is shown for contrast.
  set.seed(31)
  n <- 20
  n_test_loci <- 3000
  F0 <- 0.15
  p0 <- 0.35
  q <- 0.30
  simulate_missing <- function() {
    degraded <- stats::runif(n_test_loci) < q
    ho <- he <- rep(NA_real_, n_test_loci)
    for (j in seq_len(n_test_loci)) {
      n_j <- if (degraded[j]) sample(2:(n - 1), 1) else n
      ibd <- stats::runif(n_j) < F0
      a1 <- stats::runif(n_j) < p0
      a2 <- ifelse(ibd, a1, stats::runif(n_j) < p0)
      p <- (sum(a1) + sum(a2)) / (2 * n_j)
      ho[j] <- mean(a1 != a2)
      he[j] <- hs_biallelic(p, ho[j], n_j)
    }
    c(ratio_of_sums = fis_ratio_of_sums(ho, he),
      mean_of_ratios = mean(1 - ho / he, na.rm = TRUE))
  }
  v <- rowMeans(replicate(300, simulate_missing()))
  say(sprintf("  true FIS = %.3f | ratio-of-sums %.4f (%+.1f%% off) | mean-of-ratios %.4f (%+.1f%% off)\n",
              F0, v[["ratio_of_sums"]], 100 * (v[["ratio_of_sums"]] / F0 - 1),
              v[["mean_of_ratios"]], 100 * (v[["mean_of_ratios"]] / F0 - 1)))
  check("ratio-of-sums FIS within 5% of true FIS under locus-driven missingness",
        abs(v[["ratio_of_sums"]] / F0 - 1) < 0.05, TRUE)
  say("  (mean-of-ratios is shown for contrast, not required to pass -- it is the\n")
  say("  aggregation this package deliberately does not use.)\n")

  results <- rec$results()
  say("\n", if (all(results$pass)) "ALL CHECKS PASSED\n" else "*** SOME CHECKS FAILED ***\n", sep = "")
  invisible(results)
}

#' Calibration self-test for het_between_pops()
#'
#' Runs the type I error simulation behind `Rscript het_between_pops.R
#' --selftest`: simulations of two populations with IDENTICAL true
#' heterozygosity, each with its own allele frequencies, comparing how often
#' each test rejects at the 5% level. The tests are the combined test that
#' [het_between_pops()] reports (individuals and loci), Welch's t and
#' Wilcoxon (individuals only), and a bootstrap over loci (loci only).
#'
#' Three settings: individuals differ in inbreeding but the populations are
#' not differentiated (a locus bootstrap goes wrong); the populations are
#' differentiated (FST = 0.1) but individuals are alike (Welch's t goes
#' wrong); and both. The full grid is in `inst/sims/het_test_null.R`.
#'
#' @param verbose Print the report. Default `TRUE`.
#' @return Invisibly, a data frame with one row per test and setting:
#'   `method`, `fst` (differentiation between the populations), `sd_F` (the
#'   spread of inbreeding among individuals) and `rejection_rate` (the target
#'   is 0.05). Your own random-number state is left as it was.
#' @examples
#' \donttest{
#' het_between_pops_selftest()   # a few seconds of simulation
#' }
#' @export
het_between_pops_selftest <- function(verbose = TRUE) {
  .check_flag(verbose, "verbose")
  restore_rng <- .save_rng_state()
  on.exit(restore_rng(), add = TRUE)
  say <- function(...) if (verbose) cat(...)

  n_sim <- 300
  n_loci <- 1500
  n1 <- 15
  n2 <- 10
  ## Each locus has an ancestral allele frequency; each population draws its
  ## own (Balding-Nichols, same FST for both), so both have the same expected
  ## heterozygosity. Each individual's F is drawn around 0.1 with standard
  ## deviation sd_F; its two alleles are identical by descent with that
  ## probability at each locus.
  draw_F <- function(n, sd_F) pmin(pmax(stats::rnorm(n, 0.1, sd_F), 0), 0.9)
  population <- function(p, n, fst, sd_F) {
    if (fst > 0) p <- stats::rbeta(n_loci, p * (1 - fst) / fst, (1 - p) * (1 - fst) / fst)
    a1 <- matrix(stats::runif(n_loci * n) < p, n_loci)
    a2 <- matrix(stats::runif(n_loci * n) < p, n_loci)
    ibd <- matrix(stats::runif(n_loci * n) < rep(draw_F(n, sd_F), each = n_loci), n_loci)
    a2[ibd] <- a1[ibd]
    (a1 != a2) * 1                                   # loci x individuals, heterozygous
  }
  one_dataset <- function(fst, sd_F) {
    p <- stats::runif(n_loci, 0.05, 0.95)
    h1 <- population(p, n1, fst, sd_F)
    h2 <- population(p, n2, fst, sd_F)
    typed1 <- matrix(1, n_loci, n1)
    typed2 <- matrix(1, n_loci, n2)
    locus <- .locus_variance_of_difference(h1, typed1, h2, typed2, seq_len(n_loci))
    test <- .two_sample(colMeans(h1), colMeans(h2), "a", "b", n_loci, "heterozygosity",
                        verbose = FALSE, locus = locus)
    ## Locus bootstrap of the difference in mean heterozygosity.
    locus_diff <- rowMeans(h1) - rowMeans(h2)
    boot <- vapply(seq_len(200), function(i)
      mean(locus_diff[sample.int(n_loci, n_loci, TRUE)]), numeric(1))
    boot_ci <- stats::quantile(boot, c(0.025, 0.975), names = FALSE)
    c(combined = test$p_combined < 0.05, welch = test$p_welch < 0.05,
      wilcoxon = test$p_wilcox < 0.05, locus_bootstrap = boot_ci[1] > 0 || boot_ci[2] < 0)
  }

  settings <- list(c(fst = 0, sd_F = 0.10), c(fst = 0.10, sd_F = 0), c(fst = 0.10, sd_F = 0.10))
  say(sprintf("--- two populations with IDENTICAL true heterozygosity, n = %d and %d, %s loci ---\n",
              n1, n2, format(n_loci, big.mark = ",")))
  say(sprintf("    %d simulations per setting; %% rejected at the 5%% level (target 5%%)\n\n", n_sim))
  say(sprintf("  %-34s %9s %9s %9s %16s\n", "setting", "combined", "Welch t", "Wilcoxon",
              "locus bootstrap"))
  set.seed(404)
  results <- lapply(settings, function(s) {
    rates <- rowMeans(replicate(n_sim, one_dataset(s[["fst"]], s[["sd_F"]])))
    say(sprintf("  FST = %.2f, SD of F = %.2f %11s %8.1f%% %8.1f%% %8.1f%% %15.1f%%\n",
                s[["fst"]], s[["sd_F"]], "", 100 * rates[["combined"]], 100 * rates[["welch"]],
                100 * rates[["wilcoxon"]], 100 * rates[["locus_bootstrap"]]))
    data.frame(method = names(rates), fst = s[["fst"]], sd_F = s[["sd_F"]],
               rejection_rate = unname(rates))
  })
  say("\n  Loci alone go wrong when individuals differ in inbreeding; individuals alone\n")
  say("  go wrong when populations are differentiated. The combined test counts both.\n")
  out <- do.call(rbind, results)
  rownames(out) <- NULL
  invisible(out)
}
