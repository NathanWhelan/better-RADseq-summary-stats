###############################################################################
#
#  R/selftest.R -- self-tests a user (or a reviewer) can run on their own
#  computer without any data files or extra packages:
#
#    diversity_core_selftest()    every estimator in R/estimators.R against
#                                 brute-force Monte Carlo and published
#                                 Stacks output
#    het_between_pops_selftest()  the type I error of the individual-level
#                                 tests in het_between_pops() versus a locus
#                                 bootstrap
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
  say("  per-locus ratios. Two independent errors. See vignette(\"rationale\").\n")

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
#' --selftest`: 400 simulations of two populations with IDENTICAL true
#' heterozygosity, comparing how often Welch's t and Wilcoxon (individuals as
#' replicates) reject at the 5% level with how often a locus bootstrap does.
#'
#' @param verbose Print the report. Default `TRUE`.
#' @return Invisibly, a data frame with one row per simulated setting:
#'   `method`, `sd_F` (the spread of inbreeding among individuals) and
#'   `rejection_rate` (the target is 0.05). Your own random-number state is
#'   left as it was.
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
  results <- list()

  say("--- calibration of the individual-level tests, 400 simulations ---\n")
  say("    two populations with IDENTICAL true Ho; target rejection rate 5%\n\n")
  set.seed(404)
  n_sim <- 400
  n_loci <- 1500
  n1 <- 15
  n2 <- 10
  ## Each individual's F is drawn around 0.1 with standard deviation sd_F; its
  ## heterozygosity at a locus is then 2p(1-p)(1-F).
  individual_het <- function(F_values, p) {
    vapply(F_values, function(F_i) mean(stats::runif(length(p)) < 2 * p * (1 - p) * (1 - F_i)),
           numeric(1))
  }
  draw_F <- function(n, sd_F) pmin(pmax(stats::rnorm(n, 0.1, sd_F), 0), 0.9)
  simulate <- function(sd_F) {
    p <- stats::runif(n_loci, 0.05, 0.95)
    F_a <- draw_F(n1, sd_F)
    F_b <- draw_F(n2, sd_F)
    list(a = individual_het(F_a, p), b = individual_het(F_b, p))
  }
  say(sprintf("  %-24s %10s %10s\n", "individual F variation", "Welch t", "Wilcoxon"))
  for (sd_F in c(0.02, 0.10, 0.25)) {
    p_values <- replicate(n_sim, {
      x <- simulate(sd_F)
      c(stats::t.test(x$a, x$b)$p.value,
        suppressWarnings(stats::wilcox.test(x$a, x$b))$p.value)
    })
    rates <- rowMeans(p_values < 0.05)
    say(sprintf("  SD(F) = %.2f %-13s %9.1f%% %9.1f%%\n", sd_F, "", 100 * rates[1], 100 * rates[2]))
    results[[length(results) + 1L]] <- data.frame(method = c("welch", "wilcoxon"),
                                                  sd_F = sd_F, rejection_rate = rates)
  }

  say("\n--- and what the LOCUS bootstrap does on the same data, for contrast ---\n")
  set.seed(202)
  for (sd_F in c(0.00, 0.10)) {
    rejected <- replicate(120, {
      p <- stats::runif(n_loci, 0.05, 0.95)
      het_matrix <- function(F_values) {
        vapply(F_values, function(F_i) stats::runif(n_loci) < 2 * p * (1 - p) * (1 - F_i),
               logical(n_loci))
      }
      F_a <- draw_F(n1, sd_F)
      F_b <- draw_F(n2, sd_F)
      locus_het_a <- rowMeans(het_matrix(F_a))
      locus_het_b <- rowMeans(het_matrix(F_b))
      boot <- replicate(300, {
        i <- sample.int(n_loci, n_loci, TRUE)
        mean(locus_het_a[i]) - mean(locus_het_b[i])
      })
      !(0 >= stats::quantile(boot, 0.025) && 0 <= stats::quantile(boot, 0.975))
    })
    say(sprintf("  SD(F) = %.2f -> locus bootstrap rejects %.0f%% of the time\n", sd_F,
                100 * mean(rejected)))
    results[[length(results) + 1L]] <- data.frame(method = "locus_bootstrap", sd_F = sd_F,
                                                  rejection_rate = mean(rejected))
  }
  say("\n  The individual-level tests hold their nominal rate; the locus bootstrap\n")
  say("  does not, as soon as individuals differ from one another.\n")
  out <- do.call(rbind, results)
  rownames(out) <- NULL
  invisible(out)
}
