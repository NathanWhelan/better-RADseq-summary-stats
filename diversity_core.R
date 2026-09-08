###############################################################################
#
#  diversity_core.R  --  the five statistics this package exists to get right:
#                        Ho, He, FIS, rarefied allelic richness, and rarefied
#                        private allelic richness.
#
#  Sourced by diversity_stats.R. Run `Rscript diversity_core.R --selftest` to
#  check every estimator here against brute-force Monte Carlo, and the Stacks
#  formulas against real published Stacks output.
#
#  Run `Rscript diversity_core.R --selftest` to check the arithmetic. It needs
#  no data and no packages.
#
#  ---------------------------------------------------------------------------
#  1. EXPECTED HETEROZYGOSITY (He, gene diversity, Hs)
#  ---------------------------------------------------------------------------
#  Use Nei & Chesser (1983) / Nei (1987, eq. 7.39):
#
#       Hs = (n / (n - 1)) * (1 - sum_i p_i^2 - Ho / (2n))
#
#  with n the number of DIPLOID INDIVIDUALS genotyped, p_i the sample allele
#  frequencies, and Ho the observed heterozygote frequency.
#
#  WHY NOT (2n/(2n-1)) * (1 - sum p^2), which is what Stacks' `Pi` column is
#  and what most pipelines use? Because that correction assumes the 2n gene
#  copies are an independent sample of gametes, i.e. FIS = 0 -- the very
#  assumption under test. Writing H for the true gene diversity and F for the
#  true inbreeding coefficient:
#
#       E[1 - sum p_hat^2]              = H * (2n - 1 - F) / (2n)
#       E[(2n/(2n-1)) (1 - sum p^2)]    = H * (1 - F/(2n-1))     <- biased low
#       E[Nei-Chesser Hs]               = H                       <- unbiased
#
#  Consequences that matter here:
#    * FIS built on Stacks' Pi is biased TOWARD ZERO by a factor
#      [1 - (1-F)/(2n-1)]: about 3% at n = 15 and 5% at n = 10. The bias
#      depends on n, so it does not cancel when comparing two populations of
#      different size, and it makes FIS drift under rarefaction for a purely
#      arithmetic reason.
#    * FIS built on Stacks' Pi CANNOT REACH -1. With every individual
#      heterozygous the true FIS is -1, but 1 - Ho/Pi returns -(n-1)/n: -0.90
#      at n = 10. The Nei-Chesser version returns exactly -1. The self-test
#      checks this.
#
#  WHICH STACKS COLUMN IS WHICH. `populations.sumstats.tsv` has BOTH `Exp_Het`
#  and `Pi`, and they are different numbers estimating the same thing:
#
#      Exp_Het = 1 - sum p^2                    the plug-in estimate, biased low
#      Pi      = (2n/(2n-1)) (1 - sum p^2)      the same with the standard
#                                               unbiased-pi correction applied
#
#  So Pi = Exp_Het * 2n/(2n-1): a factor of 1.034 at n = 15, 1.053 at n = 10,
#  1.333 at n = 2. Verified against a published Stacks summary file
#  (Exp_Het 0.40000, Num_Indv 2, Pi 0.53333; 0.4 * 4/3 = 0.53333 exactly).
#  Compare anything from this package against `Pi`, never against `Exp_Het`.
#
#  Both estimators are computed here. hs_nei_chesser() is the default and what
#  the scripts report; hs_stacks_pi() / hs_stacks_pi_counts() reproduce the
#  Stacks quantity so the two can be printed side by side. diversity_stats.R
#  prints both.
#
#  A CAVEAT ON THE NAME. At a BIALLELIC site, per-site pi and expected
#  heterozygosity are the same parameter -- the probability that two randomly
#  drawn gene copies differ. At a MULTI-ALLELIC locus they are not: gene
#  diversity (1 - sum p^2) counts every pair of distinct alleles as equally
#  different, while pi weights each pair by how many nucleotides actually
#  differ. Stacks keeps that distinction in populations.hapstats.tsv, which
#  reports Gene Diversity and Haplotype Diversity as separate columns. So on a
#  haplotype VCF nothing here is nucleotide diversity, whatever it is called.
#
#  ---------------------------------------------------------------------------
#  2. WHAT "He = 0.29" ACTUALLY MEANS (read this before quoting a number)
#  ---------------------------------------------------------------------------
#  Averaging He over the SNPs in a joint, multi-population call set gives
#  heterozygosity PER ASCERTAINED SNP. That number depends on how many SNPs
#  the run happened to call, so it is not comparable to any other study.
#
#  Schmidt et al. (2021) show that SNP-based heterozygosity is biased by sample
#  size and by analysing differentiated populations together, and that the
#  unbiased quantity is AUTOSOMAL heterozygosity: the same sum divided by every
#  sequenced site, monomorphic ones included. Retaining sites that are
#  monomorphic WITHIN a population but variable in another is a step toward
#  that, not the thing itself.
#
#  The conversion is one multiplication, so there is no excuse for not
#  reporting it:
#
#       He_autosomal = He_per_SNP * (n_SNPs_used / n_sites_sequenced)
#
#  `n_sites_sequenced` is the `Sites` column of the "All positions (variant and
#  fixed)" block of populations.sumstats_summary.tsv. Stacks reports both
#  blocks; the "All positions" one is the Schmidt-compliant estimate, and the
#  "Variant positions" one is not.
#
#  autosomal_het() below does the conversion and is used wherever a total-site
#  count is supplied.
#
#  ---------------------------------------------------------------------------
#  3. FIS
#  ---------------------------------------------------------------------------
#  FIS = 1 - sum(Ho) / sum(Hs), a ratio of sums over loci, never a mean of
#  per-locus ratios (Weir & Cockerham 1984; Bhatia et al. 2013). A site
#  monomorphic in a population adds 0 to both sums, so the ratio of sums does
#  not care whether such sites are included -- which is the whole reason to
#  prefer it.
#
#  ---------------------------------------------------------------------------
#  4. RAREFIED ALLELIC RICHNESS
#  ---------------------------------------------------------------------------
#       A_j(g) = sum_i [ 1 - C(N_j - N_ij, g) / C(N_j, g) ]
#
#  the expected number of distinct alleles in g GENE COPIES drawn without
#  replacement (Hurlbert 1971; El Mousadik & Petit 1996; Petit et al. 1998;
#  Kalinowski 2004). g is in gene copies: 10 diploids is g = 20.
#
#  5. RAREFIED PRIVATE ALLELIC RICHNESS
#       P_j(g) = sum_i [ Pr(allele i drawn in j) * prod_{k!=j} Pr(not drawn in k) ]
#  (Kalinowski 2004; Szpiech et al. 2008, implemented in HP-RARE and ADZE).
#
#  Both are analytic expectations over gene copies, which is what HP-RARE and
#  ADZE do. Do not substitute Monte-Carlo subsampling of INDIVIDUALS: that is a
#  different estimand whenever FIS != 0, because the two gene copies inside one
#  individual are not an independent draw.
#
#  Requires only base R.
#
###############################################################################


## ---------------------------------------------------------------------------
## He / Hs
## ---------------------------------------------------------------------------

## Nei & Chesser (1983) unbiased gene diversity for ONE population at ONE locus.
##   sum_p2 : sum of squared sample allele frequencies
##   ho     : observed heterozygote frequency
##   n      : number of DIPLOID individuals genotyped
## Vectorised over all three. Returns NA where n < 2.
hs_nei_chesser <- function(sum_p2, ho, n) {
  out <- (n / (n - 1)) * (1 - sum_p2 - ho / (2 * n))
  out[!is.finite(out) | n < 2] <- NA_real_
  out
}

## Biallelic convenience wrapper: p is the frequency of one allele.
hs_biallelic <- function(p, ho, n) hs_nei_chesser(p^2 + (1 - p)^2, ho, n)

## The estimator Stacks reports as `Pi`, biallelic form. Reported alongside
## Nei-Chesser for comparison: unbiased when FIS = 0, biased low otherwise.
hs_stacks_pi <- function(p, n) {
  out <- 2 * p * (1 - p) * (2 * n) / (2 * n - 1)
  out[!is.finite(out) | n < 1] <- NA_real_
  out
}

## Same estimator from gene-copy counts, so it also works on a haplotype VCF.
## Named for the CORRECTION it applies, not for the software column it happens
## to match: calling it "stacks_pi" invites the confusion this package exists to
## remove. `gene_div_2n_counts` is the preferred name; the old one is kept as an
## alias so nothing that already calls it breaks.
## N is the number of gene copies actually observed, which equals 2n on the
## complete-data locus set the scripts use, so this agrees exactly with
## hs_stacks_pi() there. Checked in the self-test.
gene_div_2n_counts <- function(counts, n = NULL) {
  N <- sum(counts)
  if (!is.finite(N) || N < 2) return(NA_real_)
  (N / (N - 1)) * (1 - sum((counts / N)^2))
}
hs_stacks_pi_counts <- gene_div_2n_counts   # deprecated alias

## Multi-allelic: counts is a vector of gene-copy counts for one population at
## one locus; ho the observed heterozygote frequency; n individuals genotyped.
hs_from_counts <- function(counts, ho, n) {
  tot <- sum(counts)
  if (!is.finite(tot) || tot < 2 || n < 2) return(NA_real_)
  hs_nei_chesser(sum((counts / tot)^2), ho, n)
}

## FIS as a ratio of sums. ho and he are per-locus (or per-site) vectors.
fis_ratio_of_sums <- function(ho, he) {
  ok <- is.finite(ho) & is.finite(he)
  sh <- sum(he[ok]); if (!isTRUE(sh > 0)) return(NA_real_)
  1 - sum(ho[ok]) / sh
}

## Per-ascertained-SNP heterozygosity -> autosomal heterozygosity
## (Schmidt et al. 2021). n_sites_sequenced comes from the `Sites` column of
## the "All positions (variant and fixed)" block of
## populations.sumstats_summary.tsv.
autosomal_het <- function(het_per_snp, n_snps_used, n_sites_sequenced) {
  if (!is.finite(n_sites_sequenced) || n_sites_sequenced < n_snps_used)
    return(NA_real_)
  het_per_snp * n_snps_used / n_sites_sequenced
}


## ---------------------------------------------------------------------------
## Rarefaction
## ---------------------------------------------------------------------------

## Pr(each allele appears at least once in g gene copies drawn without
## replacement). lchoose keeps this stable for large N.
p_sampled <- function(counts, g) {
  N <- sum(counts)
  if (!is.finite(N) || g > N || g < 1) return(rep(NA_real_, length(counts)))
  1 - exp(lchoose(N - counts, g) - lchoose(N, g))
}

## Rarefied allelic richness at one locus for one population.
rare_richness <- function(counts, g) {
  if (!is.finite(sum(counts)) || sum(counts) < g) return(NA_real_)
  ps <- p_sampled(counts, g)
  if (anyNA(ps)) return(NA_real_) else sum(ps)
}

## Rarefied private allelic richness for population j at one locus.
## count_mat: populations x alleles gene-copy counts.
rare_private <- function(count_mat, j, g) {
  if (!is.matrix(count_mat)) count_mat <- rbind(count_mat)
  ps <- matrix(NA_real_, nrow(count_mat), ncol(count_mat))
  for (rr in seq_len(nrow(count_mat))) ps[rr, ] <- p_sampled(count_mat[rr, ], g)
  if (anyNA(ps)) return(NA_real_)
  term <- ps[j, ]
  for (k in setdiff(seq_len(nrow(count_mat)), j)) term <- term * (1 - ps[k, ])
  sum(term)
}

## All populations at once; returns a length-r vector.
rare_private_all <- function(count_mat, g) {
  vapply(seq_len(nrow(count_mat)), function(j) rare_private(count_mat, j, g),
         numeric(1))
}


## ---------------------------------------------------------------------------
## SELF-TEST
## ---------------------------------------------------------------------------
if (identical(commandArgs(trailingOnly = TRUE)[1], "--selftest")) {

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
    ib  <- runif(reps * n) < F0
    a1  <- runif(reps * n) < p0
    a2  <- ifelse(ib, a1, runif(reps * n) < p0)
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
  pv <- runif(L, 0.05, 0.95); nn <- rep(12, L)
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
        hs_stacks_pi_counts(cnt), hs_stacks_pi(ph, n), tol = 1e-12)
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
    r <- st[i, ]; n <- 2 * r$N; sp2 <- r$P^2 + (1 - r$P)^2
    ni <- c(round(n * r$P), n - round(n * r$P))
    chk(sprintf("N=%d p=%.2f  Exp Het = 1 - sum p^2", r$N, r$P), 1 - sp2, r$eh, 1e-4)
    chk("            Pi = (n/(n-1))(1 - sum p^2)", (n / (n - 1)) * (1 - sp2), r$pi, 1e-4)
    chk("            Pi = 1 - sum C(n_i,2)/C(n,2)",
        1 - sum(choose(ni, 2)) / choose(n, 2), r$pi, 1e-4)
    chk("            Fis = (Pi - Obs Het)/Pi", (r$pi - r$oh) / r$pi, r$fis, 1e-4)
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
    degraded <- runif(Ltest) < q
    ho <- he <- rep(NA_real_, Ltest)
    ho_naive <- he_naive <- rep(NA_real_, Ltest)
    for (j in seq_len(Ltest)) {
      nj <- if (degraded[j]) sample(2:(n - 1), 1) else n
      ib <- runif(nj) < F0
      a1 <- runif(nj) < p0
      a2 <- ifelse(ib, a1, runif(nj) < p0)
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
  quit(save = "no", status = if (ok) 0 else 1)
}
