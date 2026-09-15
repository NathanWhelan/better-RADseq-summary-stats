###############################################################################
#
#  inst/sims/het_test_null.R -- does het_between_pops() keep its 5% false
#  positive rate when the two populations have their OWN allele frequencies?
#
#  THE QUESTION. het_between_pops() gives each individual one heterozygosity
#  value. Welch's t on those values counts how much individuals vary. It does
#  not count that the typed loci are a sample of the genome. Two populations
#  with the same genome-wide diversity still differ a little at any given set
#  of loci, because drift acts on each locus separately. A simulation in which
#  both populations share the same allele frequency at every locus cannot
#  show this.
#
#  HOW IT IS ANSWERED. Balding-Nichols model: each locus has an ancestral
#  allele frequency p; each population draws its own frequency from
#  Beta(p (1 - FST) / FST, (1 - p) (1 - FST) / FST). Both populations then
#  have the same expected He, 2p(1 - p)(1 - FST), and the same distribution
#  of individual inbreeding F, so the null hypothesis (equal mean observed
#  heterozygosity) is true. FST = 0 means both use p itself. Every "5%" test
#  that rejects is a false positive.
#
#  het_between_pops() reports, for each pair, the combined test (p_combined:
#  individuals plus a jackknife over RAD loci, without counting genotype
#  noise twice; see .locus_variance_of_difference()) and, for comparison,
#  Welch's t and Wilcoxon on individuals only. This script checks all of them.
#  (An earlier run of this script, before the combined test was added to the
#  package, computed it by hand and found the same rates.)
#
#  Settings: FST 0, 0.05, 0.2; SD of individual F 0, 0.05, 0.1 (mean 0.1,
#  truncated at 0); 15 + 10 and 30 + 30 individuals; 1,000 and 10,000 loci.
#
#  Run from the package source:
#    Rscript inst/sims/het_test_null.R [reps]
#  Default 400 repeats per setting (a rate of 5% is then +/- 1.1%), on up to
#  8 cores (unix). About 10 minutes on 8 cores. Needs the package installed.
#
###############################################################################

if (!requireNamespace("RADdiversity", quietly = TRUE)) stop("Install RADdiversity first.")
suppressMessages(library(RADdiversity))
args <- commandArgs(trailingOnly = TRUE)
reps <- if (length(args) && grepl("^[0-9]+$", args[1])) as.integer(args[1]) else 400L
cores <- if (.Platform$OS.type == "unix") min(8L, parallel::detectCores()) else 1L
options(width = 200)

## Individual inbreeding coefficients: normal around 0.1, truncated to [0, 0.9].
draw_F <- function(m, sdF) pmin(pmax(stats::rnorm(m, 0.1, sdF), 0), 0.9)

## Genotypes of n individuals at loci with allele-2 frequencies `p`:
## loci x individuals matrices of allele numbers 1/2. Each individual's two
## alleles are identical by descent at a locus with probability F.
genotypes <- function(p, n, sdF) {
  L <- length(p)
  F_ind <- draw_F(n, sdF)
  a1 <- 1L + (matrix(stats::runif(L * n), L) < p)
  a2 <- 1L + (matrix(stats::runif(L * n), L) < p)
  ibd <- matrix(stats::runif(L * n) < rep(F_ind, each = L), L)
  a2[ibd] <- a1[ibd]
  list(a1 = a1, a2 = a2)
}

## One dataset of two populations as a raddiv_vcf object and popmap.
make_data <- function(L, n1, n2, fst, sdF) {
  p <- stats::runif(L, 0.05, 0.95)
  pop_freq <- function() {
    if (fst == 0) return(p)
    s <- (1 - fst) / fst
    stats::rbeta(L, p * s, (1 - p) * s)
  }
  A <- genotypes(pop_freq(), n1, sdF)
  B <- genotypes(pop_freq(), n2, sdF)
  samples <- c(paste0("a", seq_len(n1)), paste0("b", seq_len(n2)))
  a1 <- cbind(A$a1, B$a1)
  a2 <- cbind(A$a2, B$a2)
  storage.mode(a1) <- storage.mode(a2) <- "integer"
  dimnames(a1) <- dimnames(a2) <- list(NULL, samples)
  loci <- paste0("l", seq_len(L))
  H <- structure(list(A1 = a1, A2 = a2, locus = loci, locus_raw = loci,
                      alleles = rep(list(c("A", "G")), L), n_alleles = rep(2L, L),
                      samples = samples), class = "raddiv_vcf")
  list(H = H, pops = list(popA = samples[seq_len(n1)], popB = samples[n1 + seq_len(n2)]))
}

## p-values from one dataset, all from het_between_pops() itself: the
## combined test (primary), Welch's t and Wilcoxon on heterozygosity, and the
## combined and Welch tests on F. `locus_share` is the locus part of the
## combined variance as a share of Welch's variance.
one_study <- function(L, n1, n2, fst, sdF) {
  d <- make_data(L, n1, n2, fst, sdF)
  r <- suppressWarnings(het_between_pops(d$H, d$pops, nboot_g2 = 0, verbose = FALSE))
  het <- r$pairwise_tests
  h <- d$H$A1 != d$H$A2
  welch_var <- stats::var(colMeans(h[, d$pops$popA])) / n1 +
    stats::var(colMeans(h[, d$pops$popB])) / n2
  c(combined = het$p_combined, welch = het$p_welch, wilcoxon = het$p_wilcox,
    F_combined = r$pairwise_F_tests$p_combined, F_welch = r$pairwise_F_tests$p_welch,
    locus_share = het$se_combined^2 / welch_var - 1)
}

grid <- expand.grid(L = c(1000L, 10000L), n = c("15+10", "30+30"), sdF = c(0, 0.05, 0.1),
                    fst = c(0, 0.05, 0.2), stringsAsFactors = FALSE)
run_setting <- function(i) {
  s <- grid[i, ]
  n1 <- if (s$n == "15+10") 15L else 30L
  n2 <- if (s$n == "15+10") 10L else 30L
  set.seed(5000 + i)
  out <- t(replicate(reps, one_study(s$L, n1, n2, s$fst, s$sdF)))
  rate <- function(col) round(100 * mean(out[, col] < 0.05), 1)
  data.frame(FST = s$fst, sd_F = s$sdF, individuals = s$n, loci = s$L,
             combined = rate("combined"), welch = rate("welch"), wilcoxon = rate("wilcoxon"),
             F_combined = rate("F_combined"), F_welch = rate("F_welch"),
             locus_var_share = round(mean(out[, "locus_share"]), 2))
}
res <- do.call(rbind, parallel::mclapply(seq_len(nrow(grid)), run_setting, mc.cores = cores))
res <- res[order(res$FST, res$sd_F, res$individuals, res$loci), ]
cat(sprintf(paste0("\n== het_between_pops() under the null, populations with their own ",
                   "frequencies: %% of %d repeats with p < 0.05 (target 5) ==\n"), reps))
cat("   combined = p_combined (primary); welch, wilcoxon = individuals only\n")
cat("   locus_var_share = locus part of the combined variance / Welch variance (mean over repeats)\n\n")
print(res, row.names = FALSE)
