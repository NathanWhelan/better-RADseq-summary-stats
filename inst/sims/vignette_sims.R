###############################################################################
#
#  inst/sims/vignette_sims.R -- the simulations quoted in the vignettes.
#
#  Each part reproduces one table or number in vignette("rationale") or
#  vignette("reviewer-faq"). Run all of them, or name the parts wanted:
#
#    Rscript inst/sims/vignette_sims.R                  # everything
#    Rscript inst/sims/vignette_sims.R file_choice power
#
#  Parts:
#    file_choice   FIS precision from the haplotype VCF, all SNPs, and one
#                  SNP per RAD locus
#    snp_density   He from SNP and haplotype records as SNPs per locus vary,
#                  and the between-population ratio
#    denominator   pooled vs within-population denominator for He, and
#                  thinning to one SNP per locus
#    boot_modes    coverage of 95% intervals from boot = "loci",
#                  "individuals" and "both"
#    power         power of Welch's t on individual heterozygosity
#
#  Needs the package installed (or devtools::load_all()). Every part sets its
#  own seed. Run time: a few minutes per part, boot_modes longest.
#
###############################################################################

suppressMessages(library(RADdiversity))
args <- commandArgs(trailingOnly = TRUE)
parts <- if (length(args)) args else c("file_choice", "snp_density", "denominator",
                                       "boot_modes", "power")

## ---------------------------------------------------------------------------
## Simulation engine: RAD loci with haplotypes, two populations
## ---------------------------------------------------------------------------

## One dataset. Each RAD locus carries 1..max_snps SNPs and a pool of up to 8
## haplotypes with random frequencies. Population B either shares A's
## frequencies or is founded by `founders` individuals drawn from A (a
## bottleneck: fewer haplotypes, lower diversity). An individual's two
## haplotypes are identical by descent at a locus with probability F (drawn
## once per individual per locus, as inbreeding acts on whole haplotypes).
## Returns the SNP, haplotype and one-SNP-per-locus data objects and a popmap.
sim_rad <- function(n = c(A = 15L, B = 10L), n_tags = 1200L, max_snps = 4L, F = 0.1,
                    founders = NULL, snp_weights = NULL) {
  ids <- c(sprintf("A%02d", seq_len(n[[1]])), sprintf("B%02d", seq_len(n[[2]])))
  pop <- rep(c("A", "B"), n)
  n_ind <- length(ids)
  F_ind <- rep_len(F, n_ind)
  if (is.null(snp_weights)) snp_weights <- rev(seq_len(max_snps))
  snp_a1 <- snp_a2 <- hap_a1 <- hap_a2 <- vector("list", n_tags)
  n_hap_alleles <- integer(n_tags)
  for (t in seq_len(n_tags)) {
    k <- sample.int(max_snps, 1L, prob = snp_weights)
    m <- min(2L^k, 8L)
    repeat {
      haps <- unique(matrix(stats::rbinom(m * k, 1L, 0.5), ncol = k))
      if (nrow(haps) >= 2L && all(colSums(haps) %in% seq_len(nrow(haps) - 1L))) break
    }
    m <- nrow(haps)
    f_A <- stats::rgamma(m, 0.6)
    f_A <- f_A / sum(f_A)
    f_B <- if (is.null(founders)) f_A else {
      copies <- tabulate(sample.int(m, 2L * founders, replace = TRUE, prob = f_A), m)
      copies / sum(copies)
    }
    freq <- rbind(A = f_A, B = f_B)
    h1 <- h2 <- integer(n_ind)
    for (i in seq_len(n_ind)) {
      h1[i] <- sample.int(m, 1L, prob = freq[pop[i], ])
      h2[i] <- if (stats::runif(1) < F_ind[i]) h1[i] else sample.int(m, 1L, prob = freq[pop[i], ])
    }
    snp_a1[[t]] <- haps[h1, , drop = FALSE] + 1L      # individuals x SNPs, allele 1/2
    snp_a2[[t]] <- haps[h2, , drop = FALSE] + 1L
    hap_a1[[t]] <- h1
    hap_a2[[t]] <- h2
    n_hap_alleles[t] <- m
  }
  k_per_tag <- vapply(snp_a1, ncol, integer(1))
  as_H <- function(A1, A2, locus_raw, n_alleles, single_base) {
    dimnames(A1) <- dimnames(A2) <- list(NULL, ids)
    storage.mode(A1) <- storage.mode(A2) <- "integer"
    alleles <- lapply(n_alleles, function(a)
      if (single_base) c("A", "C")[seq_len(a)] else paste0("H", seq_len(a), "AC"))
    structure(list(A1 = A1, A2 = A2, locus = paste0("r", seq_len(nrow(A1))),
                   locus_raw = locus_raw, alleles = alleles, n_alleles = n_alleles,
                   samples = ids, n_records_read = nrow(A1)), class = "raddiv_vcf")
  }
  snps <- as_H(t(do.call(cbind, snp_a1)), t(do.call(cbind, snp_a2)),
               rep(paste0("L", seq_len(n_tags)), k_per_tag), rep(2L, sum(k_per_tag)), TRUE)
  haps <- as_H(do.call(rbind, hap_a1), do.call(rbind, hap_a2), paste0("L", seq_len(n_tags)),
               n_hap_alleles, FALSE)
  first <- cumsum(c(1L, k_per_tag[-n_tags]))
  one_snp <- as_H(snps$A1[first, , drop = FALSE], snps$A2[first, , drop = FALSE],
                  snps$locus_raw[first], rep(2L, n_tags), TRUE)
  list(snps = snps, haps = haps, one_snp = one_snp,
       pops = split(ids, factor(pop, levels = c("A", "B"))))
}

run_div <- function(H, pops, g = 16, nboot = 0, ...)
  suppressWarnings(diversity_stats(H, pops, g = g, nboot = nboot, verbose = FALSE, ...))

say <- function(...) cat(sprintf(...), "\n", sep = "")

## ---------------------------------------------------------------------------
## file_choice: FIS precision by file
## ---------------------------------------------------------------------------
if ("file_choice" %in% parts) {
  set.seed(101)
  say("\n== file_choice: known FIS 0.10, n = 15 and 10, 1,200 RAD loci, 20 replicates ==")
  reps <- 20L
  rows <- replicate(reps, {
    d <- sim_rad(F = 0.1, snp_weights = c(0.45, 0.30, 0.17, 0.08))
    vapply(c("haps", "snps", "one_snp"), function(f) {
      r <- run_div(d[[f]], d$pops, nboot = 1000)$per_population
      c(fis = r$Fis[1], width = r$Fis_hi[1] - r$Fis_lo[1], records = nrow(d[[f]]$A1))
    }, numeric(3))
  })
  tab <- apply(rows, c(1, 2), mean)
  print(round(rbind(FIS = tab["fis", ], bias = tab["fis", ] - 0.1, CI_width = tab["width", ],
                    records = tab["records", ]), 4))
}

## ---------------------------------------------------------------------------
## snp_density: He from SNP and haplotype records as SNP density varies
## ---------------------------------------------------------------------------
if ("snp_density" %in% parts) {
  set.seed(202)
  say("\n== snp_density: B founded by 4 individuals of A; n = 15 and 10; 1,200 loci; 5 replicates ==")
  out <- t(vapply(c(2L, 4L, 6L, 8L), function(max_snps) {
    rowMeans(replicate(5L, {
      d <- sim_rad(F = 0, founders = 4L, max_snps = max_snps)
      s <- run_div(d$snps, d$pops)$per_population
      h <- run_div(d$haps, d$pops)$per_population
      c(mean_snps = nrow(d$snps$A1) / 1200, He_snp_A = s$He[1], He_snp_B = s$He[2],
        He_hap_A = h$He[1], He_hap_B = h$He[2], Ho_snp_A = s$Ho[1], Ho_hap_A = h$Ho[1])
    }))
  }, numeric(7)))
  out <- cbind(max_snps = c(2, 4, 6, 8), out,
               ratio_snp = out[, "He_snp_A"] / out[, "He_snp_B"],
               ratio_hap = out[, "He_hap_A"] / out[, "He_hap_B"])
  print(round(out, 4))
  say("Range across SNP densities: SNP He A %.0f%%, haplotype He A %.0f%%",
      100 * (max(out[, "He_snp_A"]) / min(out[, "He_snp_A"]) - 1),
      100 * (max(out[, "He_hap_A"]) / min(out[, "He_hap_A"]) - 1))
}

## ---------------------------------------------------------------------------
## denominator: pooled vs within-population He; thinning
## ---------------------------------------------------------------------------
if ("denominator" %in% parts) {
  set.seed(303)
  say("\n== denominator: B founded by 3 individuals of A; n = 15 and 10; 1,200 loci; 5 replicates ==")
  res <- t(replicate(5L, {
    d <- sim_rad(F = 0, founders = 3L)
    pooled <- run_div(d$snps, d$pops)$per_population
    within <- pooled$He / (pooled$pct_poly / 100)
    thin <- run_div(filter_thin_one_snp(d$snps, method = "random", verbose = FALSE), d$pops)$per_population
    c(pooled_A = pooled$He[1], pooled_B = pooled$He[2], within_A = within[1], within_B = within[2],
      pct_poly_B = pooled$pct_poly[2], thin_A = thin$He[1], thin_B = thin$He[2])
  }))
  res <- cbind(res, ratio_pooled = res[, "pooled_A"] / res[, "pooled_B"],
               ratio_within = res[, "within_A"] / res[, "within_B"],
               ratio_thin = res[, "thin_A"] / res[, "thin_B"])
  print(round(res, 4))
  print(round(colMeans(res), 4))
}

## ---------------------------------------------------------------------------
## boot_modes: coverage of bootstrap intervals
## ---------------------------------------------------------------------------
if ("boot_modes" %in% parts) {
  set.seed(404)
  say("\n== boot_modes: biallelic loci with fixed frequencies; F = 0.1 for everyone; 150 replicates ==")
  L <- 1500L
  p <- stats::runif(L, 0.1, 0.9)
  truth <- c(He = mean(2 * p * (1 - p)), Fis = 0.1)
  one <- function(n, mode) {
    ids <- c(sprintf("A%02d", seq_len(n)), sprintf("B%02d", seq_len(n)))
    a1 <- matrix(1L + (stats::runif(L * 2 * n) < p), L)
    ibd <- matrix(stats::runif(L * 2 * n) < 0.1, L)
    a2 <- ifelse(ibd, a1, 1L + (stats::runif(L * 2 * n) < p))
    storage.mode(a2) <- "integer"
    dimnames(a1) <- dimnames(a2) <- list(NULL, ids)
    H <- structure(list(A1 = a1, A2 = a2, locus = paste0("r", seq_len(L)),
                        locus_raw = paste0("r", seq_len(L)), alleles = rep(list(c("A", "C")), L),
                        n_alleles = rep(2L, L), samples = ids), class = "raddiv_vcf")
    pops <- list(A = ids[seq_len(n)], B = ids[n + seq_len(n)])
    r <- run_div(H, pops, g = 2 * n, nboot = 200, boot = mode)$per_population
    c(He = truth[["He"]] >= r$He_lo[1] && truth[["He"]] <= r$He_hi[1],
      Fis = truth[["Fis"]] >= r$Fis_lo[1] && truth[["Fis"]] <= r$Fis_hi[1])
  }
  for (n in c(10L, 15L)) for (mode in c("loci", "individuals", "both")) {
    cover <- rowMeans(replicate(150L, one(n, mode)))
    say("n = %2d  boot = %-11s  He coverage %5.1f%%   FIS coverage %5.1f%%", n, mode,
        100 * cover[["He"]], 100 * cover[["Fis"]])
  }
}

## ---------------------------------------------------------------------------
## power: Welch's t on individual heterozygosity
## ---------------------------------------------------------------------------
if ("power" %in% parts) {
  set.seed(505)
  say("\n== power: n = 15 and 10, 1,500 loci, SD of F among individuals 0.10; 1,000 replicates ==")
  L <- 1500L
  individual_het <- function(F_values, p)
    vapply(F_values, function(f) mean(stats::runif(L) < 2 * p * (1 - p) * (1 - f)), numeric(1))
  draw_F <- function(n, mean_F) pmin(pmax(stats::rnorm(n, mean_F, 0.10), 0), 0.9)
  for (d_F in c(0, 0.05, 0.08, 0.10)) {
    rejected <- replicate(1000L, {
      p <- stats::runif(L, 0.05, 0.95)
      a <- individual_het(draw_F(15, 0.10), p)
      b <- individual_het(draw_F(10, 0.10 + d_F), p)
      stats::t.test(a, b)$p.value < 0.05
    })
    say("mean F differs by %.2f (about %.1f SD of F): Welch rejects %.0f%% of the time",
        d_F, d_F / 0.10, 100 * mean(rejected))
  }
}
