###############################################################################
#
#  inst/sims/coverage_se.R -- does each standard error cover the truth?
#
#  Two populations are sampled again and again from one fixed set of loci
#  with known allele frequencies p. Each individual is inbred with
#  probability F_i (its two alleles identical by descent at a locus), with
#  F_i drawn per individual from a normal distribution with mean 0.1 and
#  SD 0, 0.1 or 0.25, truncated to [0, 0.9]. diversity_stats() then reports
#  Ho, He and FIS with the locus-jackknife SE (`_se`) and the
#  individual-jackknife SE (`_se_ind`); a 95% interval is estimate +/- 1.96 SE.
#
#  The truth is the population value for these loci:
#    He = mean over loci of 2p(1-p)      FIS = E[F]      Ho = He (1 - E[F])
#  and the question is how often each interval contains it.
#
#  The locus SE treats the sampled individuals as fixed. With exchangeable
#  individuals (SD 0) that costs little, because with many loci the error
#  from which individuals were sampled averages away; when individuals differ
#  in inbreeding it does not average away, and only the individual SE sees it.
#
#  Run from the package source:   Rscript inst/sims/coverage_se.R [reps]
#  (default 200 replicates per setting; a few minutes). Needs the package
#  installed, or run under devtools::load_all().
#
###############################################################################

if (!requireNamespace("RADdiversity", quietly = TRUE)) stop("Install RADdiversity first.")
suppressMessages(library(RADdiversity))
args <- commandArgs(trailingOnly = TRUE)
reps <- if (length(args)) as.integer(args[1]) else 200L

set.seed(20260911)
L <- 1000; n <- 15                           # loci; individuals per population
p <- stats::runif(L, 0.1, 0.9)
truth_He <- mean(2 * p * (1 - p))
samp <- c(paste0("a", seq_len(n)), paste0("b", seq_len(n)))
popmap <- tempfile()
writeLines(paste0(samp, "\t", rep(c("popA", "popB"), each = n)), popmap)

draw_F <- function(k, sdF) pmin(pmax(stats::rnorm(k, 0.1, sdF), 0), 0.9)

one_rep <- function(sdF) {
  F <- draw_F(2 * n, sdF)
  a1 <- matrix(1L + (stats::runif(L * 2 * n) < p), L)
  ibd <- matrix(stats::runif(L * 2 * n) < rep(F, each = L), L)
  a2 <- ifelse(ibd, a1, 1L + (stats::runif(L * 2 * n) < p))
  dimnames(a1) <- dimnames(a2) <- list(NULL, samp)
  H <- list(A1 = a1, A2 = a2, locus = paste0("l", seq_len(L)),
            locus_raw = paste0("l", seq_len(L)),
            alleles = replicate(L, c("A", "C"), simplify = FALSE),
            n_alleles = rep(2L, L), samples = samp)
  suppressMessages(diversity_stats(H, popmap, g = 4, nboot = 0, stem = "sim",
                                   se_individuals = TRUE))$per_population
}

rows <- list()
for (sdF in c(0, 0.1, 0.25)) {
  Fbar <- mean(draw_F(1e6, sdF))                # E[F] after truncation
  truth <- c(Ho = truth_He * (1 - Fbar), He = truth_He, Fis = Fbar)
  hit <- list()
  for (rep in seq_len(reps)) {
    pp <- one_rep(sdF)
    for (s in names(truth)) {
      est <- pp[[s]]; se_l <- pp[[paste0(s, "_se")]]; se_i <- pp[[paste0(s, "_se_ind")]]
      hit[[length(hit) + 1]] <- data.frame(
        stat = s, loci = abs(est - truth[[s]]) <= 1.96 * se_l,
        individuals = abs(est - truth[[s]]) <= 1.96 * se_i, ratio = se_i / se_l)
    }
  }
  h <- do.call(rbind, hit)
  for (s in names(truth))
    rows[[length(rows) + 1]] <- data.frame(
      SD_F = sdF, statistic = s,
      coverage_locus_SE = sprintf("%.1f%%", 100 * mean(h$loci[h$stat == s])),
      coverage_individual_SE = sprintf("%.1f%%", 100 * mean(h$individuals[h$stat == s])),
      SE_ratio_ind_to_locus = round(mean(h$ratio[h$stat == s]), 2))
}
cat(sprintf("Coverage of 95%% intervals (estimate +/- 1.96 SE), %d replicates x 2 populations,\n", reps))
cat(sprintf("%d loci, %d individuals per population; mean F = 0.1\n\n", L, n))
print(do.call(rbind, rows), row.names = FALSE)
