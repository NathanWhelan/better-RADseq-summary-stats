###############################################################################
#
#  inst/sims/g2_linked_snps.R -- does g2 find "variance in inbreeding" that
#  is not there, when several SNPs sit on one RAD locus?
#
#  THE QUESTION. g2 (identity disequilibrium; David et al. 2007) asks whether
#  an individual that is heterozygous at one locus is more often heterozygous
#  at another than two different individuals are. It assumes the two loci of
#  each pair are inherited independently. SNPs on one RAD tag are not: an
#  individual carrying two different haplotypes of the tag is heterozygous at
#  several of its SNPs at once. Counting those pairs looks exactly like
#  variance in inbreeding.
#
#  HOW IT IS ANSWERED. Random mating, no inbreeding (the null: true g2 = 0).
#  Each RAD tag has 1, 2 or 3 SNPs (probabilities 0.4, 0.4, 0.2; mean 1.8)
#  and 2 to 4 haplotypes with random frequencies. Each individual draws two
#  haplotypes, so SNPs on a tag are linked. 5% of (tag, individual) cells are
#  missing, the whole tag at once, as in RAD data. SNPs that do not vary in
#  the sample are dropped, as Stacks does. 20 individuals; 500, 2,000 and
#  8,000 tags. g2 is computed two ways, with the package's own function:
#    old  every pair of SNPs counts, and the permutation shuffles each SNP on
#         its own (what the package did before 2026-09-23);
#    new  pairs of SNPs on the same tag are left out, and the permutation
#         moves a tag's SNPs together (what identity_disequilibrium() and
#         het_between_pops() do now).
#  For each, the share of runs with permutation p < 0.05 and with the 95%
#  bootstrap interval above 0. Both should be near 5% and 2.5%.
#
#  The positive control: the same, but half the individuals have F = 0.3
#  (their two haplotypes are identical by descent at a tag with probability
#  0.3). True g2 = Var(F) / (1 - mean F)^2 = 0.031. g2 must still be found.
#
#  Results are quoted in vignette("rationale"), section 5 ("g2"), in the
#  header of R/identity_disequilibrium.R and in NEWS.md.
#
#  Run from the package source:
#    Rscript inst/sims/g2_linked_snps.R [reps]
#  Default 100 repeats per setting (a rate of 5% is then +/- 2.2%), on up to
#  8 cores (unix). About 15 minutes on 4 cores. Needs the package installed,
#  or devtools::load_all() and then source() this file.
#
#  Result (seeds 20260923 + ..., set in each repeat, 100 repeats):
#    null      old: p < 0.05 in 54-56% of runs, CI above 0 in 29-33%
#              new: p < 0.05 in 2-4%, CI above 0 in 0-1%
#    positive  old and new: found in 100% (mean g2 0.033-0.035, true 0.031)
#
###############################################################################

if (!requireNamespace("RADdiversity", quietly = TRUE)) stop("Install RADdiversity first.")
suppressMessages(library(RADdiversity))
args <- commandArgs(trailingOnly = TRUE)
reps <- if (length(args) && grepl("^[0-9]+$", args[1])) as.integer(args[1]) else 100L
cores <- if (.Platform$OS.type == "unix") min(8L, parallel::detectCores()) else 1L
options(width = 200)
g2_summary <- RADdiversity:::.g2_summary

## One data set: records x individuals allele matrices (1 = REF, 2 = ALT)
## and the RAD tag of each record. `F_ind` is each individual's inbreeding
## coefficient.
simulate <- function(n_tags, F_ind, missing = 0.05) {
  n <- length(F_ind)
  a1 <- a2 <- tag_of <- vector("list", n_tags)
  for (t in seq_len(n_tags)) {
    s <- sample(1:3, 1L, prob = c(0.4, 0.4, 0.2))           # SNPs on the tag
    k <- sample(2:4, 1L)                                    # its haplotypes
    haps <- matrix(stats::rbinom(k * s, 1L, 0.5), k, s)
    haps[1L, ] <- 0L
    haps[2L, ] <- 1L                                        # every SNP varies
    freq <- stats::rgamma(k, 1)
    h1 <- sample.int(k, n, replace = TRUE, prob = freq)
    h2 <- sample.int(k, n, replace = TRUE, prob = freq)
    ibd <- stats::runif(n) < F_ind
    h2[ibd] <- h1[ibd]
    g1 <- 1L + t(haps[h1, , drop = FALSE])                  # SNPs x individuals
    g2 <- 1L + t(haps[h2, , drop = FALSE])
    gone <- stats::runif(n) < missing                       # the whole tag
    g1[, gone] <- NA
    g2[, gone] <- NA
    a1[[t]] <- g1
    a2[[t]] <- g2
    tag_of[[t]] <- rep(t, s)
  }
  A1 <- do.call(rbind, a1)
  A2 <- do.call(rbind, a2)
  tag <- unlist(tag_of)
  varies <- vapply(seq_len(nrow(A1)), function(j) {
    x <- c(A1[j, ], A2[j, ])
    length(unique(x[!is.na(x)])) > 1L
  }, logical(1))
  list(A1 = A1[varies, , drop = FALSE], A2 = A2[varies, , drop = FALSE], tag = tag[varies])
}

one_run <- function(n_tags, F_ind) {
  d <- simulate(n_tags, F_ind)
  old <- g2_summary(d$A1, d$A2, nboot = 200L, nperm = 199L)
  new <- g2_summary(d$A1, d$A2, nboot = 200L, nperm = 199L, block = d$tag)
  c(snps_per_tag = nrow(d$A1) / length(unique(d$tag)),
    old_g2 = old$g2, old_p = old$p_value, old_lo = old$g2_lo,
    new_g2 = new$g2, new_p = new$p_value, new_lo = new$g2_lo)
}

## Each repeat sets its own seed, so the results are the same on any number
## of cores (a seed set once before mclapply() does not reach the workers).
n_ind <- 20L
settings <- expand.grid(tags = c(500L, 2000L, 8000L), model = c("null", "positive"),
                        stringsAsFactors = FALSE)
rows <- lapply(seq_len(nrow(settings)), function(i) {
  s <- settings[i, ]
  F_ind <- if (s$model == "null") rep(0, n_ind) else rep(c(0, 0.3), each = n_ind / 2)
  runs <- parallel::mclapply(seq_len(reps), function(r) {
    set.seed(20260923 + 1000L * i + r)
    one_run(s$tags, F_ind)
  }, mc.cores = cores)
  runs <- do.call(rbind, runs)
  data.frame(model = s$model, tags = s$tags, reps = reps,
             snps_per_tag = round(mean(runs[, "snps_per_tag"]), 2),
             old_mean_g2 = signif(mean(runs[, "old_g2"]), 2),
             old_p_below_05 = sprintf("%.0f%%", 100 * mean(runs[, "old_p"] < 0.05)),
             old_ci_above_0 = sprintf("%.0f%%", 100 * mean(runs[, "old_lo"] > 0)),
             new_mean_g2 = signif(mean(runs[, "new_g2"]), 2),
             new_p_below_05 = sprintf("%.0f%%", 100 * mean(runs[, "new_p"] < 0.05)),
             new_ci_above_0 = sprintf("%.0f%%", 100 * mean(runs[, "new_lo"] > 0)))
})
cat("g2 with linked SNPs: 20 individuals. null = no inbreeding (true g2 = 0);",
    "positive = half the individuals at F = 0.3 (true g2 = 0.031).\n")
cat("old = every pair of SNPs; new = pairs on the same RAD tag left out.\n\n")
print(do.call(rbind, rows), row.names = FALSE)
