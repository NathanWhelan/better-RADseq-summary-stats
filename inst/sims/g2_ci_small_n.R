###############################################################################
#
#  inst/sims/g2_ci_small_n.R -- how well do g2's bootstrap interval and
#  permutation test behave with few individuals?
#
#  THE QUESTION. identity_disequilibrium() gives g2 a 95% interval from a
#  bootstrap over individuals, and a p-value from permutations (David et al.
#  2007; inbreedR does the same). A bootstrap over individuals has only the
#  sampled individuals to work with, and RAD-seq studies often have 8-15 per
#  population.
#
#  HOW IT IS ANSWERED. Random mating, no inbreeding (true g2 = 0), unlinked
#  biallelic loci with allele frequency uniform on 0.05-0.5, no missing data.
#  8, 15 and 30 individuals; 400 and 2,000 loci; 300 repeats each. g2 comes
#  from the package's own function, with 200 bootstrap replicates and 199
#  permutations. For each setting: how often the 95% interval misses the
#  true 0 (should be 5%), on which side, how often the permutation p is
#  below 0.05 (should be 5%), and the bootstrap SE against the real spread
#  of g2 over repeats.
#
#  Results are quoted in ?identity_disequilibrium and vignette("rationale"),
#  section 5.
#
#  Run from the package source:
#    Rscript inst/sims/g2_ci_small_n.R [reps]
#  About 5 minutes on 4 cores (unix). Needs the package installed, or
#  devtools::load_all() and then source() this file.
#
#  Result (seeds 20260923 + ..., 300 repeats per setting):
#     n  loci  interval misses 0  below 0  above 0  permutation p < 0.05  boot SE / true SD
#     8   400                17%      16%       1%                    2%               0.78
#     8  2000                19%      19%       0%                    3%               0.72
#    15   400                15%      15%       1%                    4%               0.82
#    15  2000                10%      10%       1%                    4%               0.89
#    30   400                 7%       7%       0%                    3%               0.99
#    30  2000                 7%       6%       1%                    5%               0.99
#  So the interval is too narrow below about 30 individuals, and misses on
#  the low side; the permutation test holds its level at every size. (A 5%
#  rate over 300 repeats has a Monte Carlo SE of about 1.3%.)
#
###############################################################################

if (!requireNamespace("RADdiversity", quietly = TRUE)) stop("Install RADdiversity first.")
suppressMessages(library(RADdiversity))
args <- commandArgs(trailingOnly = TRUE)
reps <- if (length(args) && grepl("^[0-9]+$", args[1])) as.integer(args[1]) else 300L
cores <- if (.Platform$OS.type == "unix") min(8L, parallel::detectCores()) else 1L
g2_summary <- RADdiversity:::.g2_summary
options(width = 200)

one_run <- function(n, loci) {
  p <- stats::runif(loci, 0.05, 0.5)
  a1 <- matrix(1L + (stats::runif(loci * n) < p), loci)
  a2 <- matrix(1L + (stats::runif(loci * n) < p), loci)
  s <- g2_summary(a1, a2, nboot = 200L, nperm = 199L)
  c(g2 = s$g2, se = s$g2_se, lo = s$g2_lo, hi = s$g2_hi, p = s$p_value)
}

## Each repeat sets its own seed, so the results are the same on any number
## of cores (a seed set once before mclapply() does not reach the workers).
settings <- expand.grid(loci = c(400L, 2000L), n = c(8L, 15L, 30L))
rows <- lapply(seq_len(nrow(settings)), function(i) {
  s <- settings[i, ]
  runs <- parallel::mclapply(seq_len(reps), function(r) {
    set.seed(20260923 + 1000L * i + r)
    one_run(s$n, s$loci)
  }, mc.cores = cores)
  runs <- do.call(rbind, runs)
  pct <- function(x) sprintf("%.0f%%", 100 * mean(x))
  data.frame(n = s$n, loci = s$loci, reps = reps,
             misses_0 = pct(runs[, "lo"] > 0 | runs[, "hi"] < 0),
             below_0 = pct(runs[, "hi"] < 0), above_0 = pct(runs[, "lo"] > 0),
             permutation_p_below_05 = pct(runs[, "p"] < 0.05),
             boot_se_over_true_sd = round(mean(runs[, "se"]) / stats::sd(runs[, "g2"]), 2))
})
cat("g2 with no inbreeding (true g2 = 0): 95% bootstrap interval over individuals",
    "and permutation test.\n\n")
print(do.call(rbind, rows), row.names = FALSE)
