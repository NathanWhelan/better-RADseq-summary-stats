###############################################################################
#
#  inst/sims/kinship_checks.R -- the numbers quoted in ?kinship_check, in the
#  header of R/kinship.R and in vignette("rationale"), section 7.
#
#  PARTS (name them on the command line, or run both):
#    dropout      KING's two estimators when one individual's heterozygotes
#                 are called homozygous (allele dropout at low depth)
#    pooled_beta  Goudet's beta on two differentiated populations pooled,
#                 against each population as its own reference
#
#    Rscript inst/sims/kinship_checks.R               # both
#    Rscript inst/sims/kinship_checks.R dropout
#
#  THE TWO KING ESTIMATORS (Manichaikul et al. 2010), for individuals i, j
#  over the loci typed in both:
#    within-family (kinship_check())   (N11 - 2 N20) / (N1_i + N1_j)
#    between-family (PLINK 2's --make-king)
#                         1/2 - (4 N20 + N10 + N01) / (4 min(N1_i, N1_j))
#  N11: both heterozygous. N20: opposite homozygotes. N10, N01: one
#  heterozygous, the other homozygous. N1_i = N11 + N10, i's heterozygous
#  loci. The between-family form is written here from Manichaikul et al.
#  (2010). PLINK 2's source (ComputeKinship() in plink2_matrix_calc.cc,
#  GPL-3) was read only to confirm that --make-king computes this form; no
#  code was copied.
#
#  dropout. 20,000 unlinked biallelic loci, allele frequency uniform on
#  0.05-0.95, 10 repeats. A parent, its offspring (one allele from the parent,
#  one from the population) and an unrelated individual. Then 20% of the
#  offspring's (or the unrelated individual's) heterozygotes are called
#  homozygous for one of its two alleles, at random. For the parent-offspring
#  pair the expected values are exact, whatever the allele frequencies. With
#  H heterozygous loci per individual before dropout: N11 = 0.4 H,
#  N20 = 0.05 H, N10 = 0.6 H, N01 = 0.4 H, N1 = H and 0.8 H. So within-family
#  (0.4 - 0.1)/1.8 = 1/6 = 0.167 and between-family 1/2 - 1.2/3.2 = 1/8 =
#  0.125, against the true 0.25.
#
#  pooled_beta. Two populations at FST = 0.1 (Balding-Nichols), 20
#  individuals each, 5,000 unlinked loci, 20 repeats. Population A holds 16
#  unrelated individuals, a full-sib pair and a half-sib pair. Beta
#  (kinship_check(method = "beta"), i.e. hierfstat::beta.dosage()) is
#  computed on both populations pooled (no popmap) and within each
#  population (with popmap). Counted: the 378 unrelated same-population pairs
#  above 0.0442 (the third-degree cutpoint), and the planted pairs' values
#  (true kinship 0.25 and 0.125).
#
#  Needs the package installed (or devtools::load_all()) and, for
#  pooled_beta, hierfstat.
#
#  Result (seed 20260923):
#    dropout      parent-offspring, no dropout      0.250 within, 0.249 between
#                 parent-offspring, 20% dropout     0.165 within, 0.121 between
#                 unrelated, 20% dropout            -0.107 within, -0.187 between
#    pooled_beta  unrelated same-population pairs above 0.0442:
#                 pooled 80% (74-90% over repeats); within populations 0%
#                 full sibs 0.254, half sibs 0.124 (within populations)
#
###############################################################################

suppressMessages(library(RADdiversity))
args <- commandArgs(trailingOnly = TRUE)
parts <- if (length(args)) args else c("dropout", "pooled_beta")

## A data object kinship_check() accepts, from records x individuals matrices
## of ALT allele copies (0/1) for each individual's two alleles.
as_H <- function(first, second, samples) {
  n_rec <- nrow(first)
  A1 <- 1L + first
  A2 <- 1L + second
  dimnames(A1) <- dimnames(A2) <- list(NULL, samples)
  list(A1 = A1, A2 = A2, locus = paste0("L", seq_len(n_rec)),
       locus_raw = paste0("L", seq_len(n_rec)),
       alleles = replicate(n_rec, c("A", "C"), simplify = FALSE),
       n_alleles = rep(2L, n_rec), samples = samples)
}

## One individual drawn from allele frequencies p: its two alleles.
draw <- function(p) cbind(stats::runif(length(p)) < p, stats::runif(length(p)) < p) * 1L

## A child of two parents (each a two-column allele matrix): one allele from
## each, chosen at random at every locus.
child_of <- function(mother, father) {
  L <- nrow(mother)
  pick <- function(parent) parent[cbind(seq_len(L), sample.int(2L, L, replace = TRUE))]
  cbind(pick(mother), pick(father))
}

if ("dropout" %in% parts) {
  set.seed(20260923)
  L <- 20000L
  ## 20% of an individual's heterozygotes called homozygous for one of its
  ## two alleles, at random.
  drop_hets <- function(g, rate = 0.2) {
    het <- which(g[, 1] != g[, 2])
    hit <- het[stats::runif(length(het)) < rate]
    keep <- sample.int(2L, length(hit), replace = TRUE)
    g[cbind(hit, 3L - keep)] <- g[cbind(hit, keep)]
    g
  }
  between_family <- function(gi, gj) {
    di <- rowSums(gi)
    dj <- rowSums(gj)
    n20 <- sum(abs(di - dj) == 2L)
    n10 <- sum(di == 1L & dj != 1L)
    n01 <- sum(dj == 1L & di != 1L)
    0.5 - (4 * n20 + n10 + n01) / (4 * min(sum(di == 1L), sum(dj == 1L)))
  }
  within_family <- function(gi, gj) {
    H <- as_H(cbind(gi[, 1], gj[, 1]), cbind(gi[, 2], gj[, 2]), c("i", "j"))
    kinship_check(H, threshold = NULL, verbose = FALSE)$pairwise$kinship
  }
  runs <- t(replicate(10L, {
    p <- stats::runif(L, 0.05, 0.95)
    parent <- draw(p)
    offspring <- child_of(parent, draw(p))
    unrelated <- draw(p)
    offspring_d <- drop_hets(offspring)
    unrelated_d <- drop_hets(unrelated)
    c(po_within = within_family(parent, offspring),
      po_between = between_family(parent, offspring),
      po_drop_within = within_family(parent, offspring_d),
      po_drop_between = between_family(parent, offspring_d),
      un_drop_within = within_family(parent, unrelated_d),
      un_drop_between = between_family(parent, unrelated_d))
  }))
  m <- colMeans(runs)
  cat("dropout: KING within-family (kinship_check) vs between-family (PLINK 2), 20,000 loci, 10 repeats\n")
  cat(sprintf("  parent-offspring, no dropout:   within %.3f, between %.3f (true 0.25)\n",
              m[["po_within"]], m[["po_between"]]))
  cat(sprintf("  parent-offspring, 20%% dropout:  within %.3f, between %.3f (exact: 0.167, 0.125)\n",
              m[["po_drop_within"]], m[["po_drop_between"]]))
  cat(sprintf("  unrelated, 20%% dropout:         within %.3f, between %.3f (true 0)\n\n",
              m[["un_drop_within"]], m[["un_drop_between"]]))
}

if ("pooled_beta" %in% parts) {
  if (!requireNamespace("hierfstat", quietly = TRUE)) stop("pooled_beta needs hierfstat.")
  set.seed(20260923)
  L <- 5000L
  fst <- 0.1
  one_run <- function() {
    p0 <- stats::runif(L, 0.05, 0.95)
    own_p <- function() stats::rbeta(L, p0 * (1 - fst) / fst, (1 - p0) * (1 - fst) / fst)
    pa <- own_p()
    pb <- own_p()
    ## Population A: 16 unrelated, a full-sib pair and a half-sib pair.
    mother <- draw(pa)
    father <- draw(pa)
    shared <- draw(pa)
    genomes <- c(lapply(1:16, function(i) draw(pa)),
                 list(child_of(mother, father), child_of(mother, father),
                      child_of(shared, draw(pa)), child_of(shared, draw(pa))),
                 lapply(1:20, function(i) draw(pb)))
    samples <- c(sprintf("a%02d", 1:20), sprintf("b%02d", 1:20))
    H <- as_H(vapply(genomes, function(g) g[, 1], numeric(L)),
              vapply(genomes, function(g) g[, 2], numeric(L)), samples)
    popmap <- list(A = samples[1:20], B = samples[21:40])
    planted <- c("a17 a18", "a19 a20")
    count <- function(pw) {
      same <- substr(pw$sample1, 1, 1) == substr(pw$sample2, 1, 1)
      pair <- paste(pw$sample1, pw$sample2)
      unrelated <- same & !pair %in% planted
      c(n_unrelated = sum(unrelated), above = sum(pw$kinship[unrelated] > 0.0442),
        full_sib = pw$kinship[pair == planted[1]], half_sib = pw$kinship[pair == planted[2]])
    }
    pooled <- count(kinship_check(H, method = "beta", threshold = NULL, verbose = FALSE)$pairwise)
    within <- count(kinship_check(H, popmap, method = "beta", threshold = NULL,
                                  verbose = FALSE)$pairwise)
    c(n = pooled[["n_unrelated"]], pooled_above = pooled[["above"]],
      within_above = within[["above"]], full_sib = within[["full_sib"]],
      half_sib = within[["half_sib"]])
  }
  runs <- t(replicate(20L, one_run()))
  share <- runs[, "pooled_above"] / runs[, "n"]
  cat("pooled_beta: two populations at FST = 0.1, 20 individuals each, 5,000 loci, 20 repeats\n")
  cat(sprintf("  unrelated same-population pairs above 0.0442, pooled: %.0f%% (%.0f-%.0f%% over repeats; %d pairs)\n",
              100 * mean(share), 100 * min(share), 100 * max(share), runs[1, "n"]))
  cat(sprintf("  the same, each population its own reference: %.1f%%\n",
              100 * mean(runs[, "within_above"] / runs[, "n"])))
  cat(sprintf("  planted pairs, within populations: full sibs %.3f, half sibs %.3f (true 0.25, 0.125)\n",
              mean(runs[, "full_sib"]), mean(runs[, "half_sib"])))
}
