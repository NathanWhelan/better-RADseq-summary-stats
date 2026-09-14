###############################################################################
#
#  inst/sims/uncertainty_sources.R -- which standard error should be reported
#  for each statistic?
#
#  THE QUESTION. A population's He (or Ho, FIS, Ar, privAr) is uncertain for
#  two reasons: the RAD loci typed are a sample of the genome, and the
#  individuals caught are a sample of the population. diversity_stats() offers
#    _se           jackknife over RAD loci (individuals held fixed)
#    _lo, _hi      bootstrap over RAD loci (individuals held fixed)
#    _se_ind       jackknife over individuals (loci held fixed)
#    _se_combined  sqrt(_se^2 + _se_ind^2), both together
#  and, for comparison only, bootstraps over individuals (boot =
#  "individuals") and over both (boot = "both").
#
#  HOW IT IS ANSWERED. Simulated populations with a KNOWN true value. Each
#  simulated study is repeated many times, and we count how often the 95%
#  interval (estimate +/- 1.96 SE, or the bootstrap's _lo to _hi) contains the
#  truth. A good method does so about 95% of the time. We also report the
#  average SE divided by the true SE (the SD of the estimates across
#  repeats): 1 is right, above 1 too wide, below 1 too narrow.
#
#  TWO KINDS OF TRUTH.
#    "genome"      new loci AND new individuals in every repeat. The truth is
#                  the population's genome-wide value. This is what a reported
#                  He means, and what the package's advice is based on.
#    "these loci"  the same loci in every repeat, new individuals only. The
#                  truth is the population's value at those loci.
#  An earlier version of this check (coverage_se.R) used only "these loci",
#  which is why it recommended the individual SE alone; against the genome
#  truth that undercovers He.
#
#  True values: He = mean over loci of 1 - sum(p^2), Ho = He (1 - E[F]),
#  FIS = E[F]. Ar and privAr have no simple formula, so their "truth" is the
#  mean estimate over all repeats; their coverage therefore checks only the
#  width of the interval, not bias.
#
#  PARTS (name them on the command line, or run all):
#    main         SD of F among individuals 0, 0.1, 0.25; both truths; every
#                 method including the individual and two-way bootstraps
#    robust       genome truth under other conditions: many rare variants,
#                 8 and 30 individuals, 3,000 loci, 20% missing genotypes,
#                 4-allele haplotype loci, no inbreeding variance, and a
#                 sample made of full-sib families
#    boundary     how often the individual SE is too large when records have
#                 only 2 genotyped individuals (calibrates the warning in
#                 diversity_stats(); see .jackknife_boundary())
#
#  Run from the package source:
#    Rscript inst/sims/uncertainty_sources.R [reps] [part ...]
#  e.g. Rscript inst/sims/uncertainty_sources.R 300 robust
#  Default 250 repeats, all parts, on up to 4 cores (unix). "main" is the
#  slowest (bootstraps over individuals): about 15 minutes at 250 repeats.
#  Needs the package installed, or devtools::load_all().
#
###############################################################################

if (!requireNamespace("RADdiversity", quietly = TRUE)) stop("Install RADdiversity first.")
suppressMessages(library(RADdiversity))
args <- commandArgs(trailingOnly = TRUE)
reps <- if (length(args) && grepl("^[0-9]+$", args[1])) as.integer(args[1]) else 250L
parts <- setdiff(args, as.character(reps))
if (!length(parts)) parts <- c("main", "robust", "boundary")
cores <- if (.Platform$OS.type == "unix") min(4L, parallel::detectCores()) else 1L
options(width = 200)

## ---------------------------------------------------------------------------
## Simulating data
## ---------------------------------------------------------------------------

## Individual inbreeding coefficients: normal around 0.1, truncated to [0, 0.9].
draw_F <- function(m, sdF) pmin(pmax(stats::rnorm(m, 0.1, sdF), 0), 0.9)
mean_F <- function(sdF) { set.seed(99); mean(draw_F(1e6, sdF)) }

## Allele frequencies for L loci (rows) and k alleles (columns).
##   uniform  biallelic, p from 0.1 to 0.9
##   rare     biallelic, U-shaped: many rare variants
##   haps     4 alleles with uneven frequencies (like RAD haplotypes)
allele_freqs <- function(L, kind) switch(kind,
  uniform = { p <- stats::runif(L, 0.1, 0.9); cbind(p, 1 - p) },
  rare    = { p <- pmin(pmax(stats::rbeta(L, 0.15, 0.15), 0.005), 0.995); cbind(p, 1 - p) },
  haps    = { x <- matrix(stats::rgamma(L * 4, 0.6), L); x / rowSums(x) })

## Genome-wide true He for a kind of frequencies, from 200,000 loci.
true_He <- function(kind) { set.seed(1); f <- allele_freqs(2e5, kind); mean(1 - rowSums(f^2)) }

## Genotypes of n individuals: loci x individuals allele matrices a1, a2.
## Each individual's two alleles are identical by descent at a locus with
## probability F. With `families` > 0 the sample is that many full-sib
## families of equal size instead (no inbreeding).
genotypes <- function(freq, n, sdF, families = 0) {
  L <- nrow(freq)
  cum <- t(apply(freq, 1, cumsum))
  draw <- function(m) {
    u <- matrix(stats::runif(L * m), L)
    a <- matrix(1L, L, m)
    for (j in seq_len(ncol(cum) - 1)) a <- a + (u > cum[, j])
    a
  }
  if (families > 0) {
    a1 <- a2 <- NULL
    for (f in seq_len(families)) {
      mother <- draw(2); father <- draw(2)
      for (s in seq_len(n / families)) {
        a1 <- cbind(a1, mother[cbind(seq_len(L), 1 + (stats::runif(L) < 0.5))])
        a2 <- cbind(a2, father[cbind(seq_len(L), 1 + (stats::runif(L) < 0.5))])
      }
    }
    return(list(a1 = a1, a2 = a2))
  }
  F_ind <- draw_F(n, sdF)
  a1 <- draw(n); a2 <- draw(n)
  ibd <- matrix(stats::runif(L * n) < rep(F_ind, each = L), L)
  a2[ibd] <- a1[ibd]
  list(a1 = a1, a2 = a2)
}

## Two populations drawn from the same frequencies, as a raddiv_vcf object.
## Each locus is its own RAD locus. `missing`: share of genotypes set missing.
make_data <- function(freq, n, sdF, families = 0, missing = 0) {
  A <- genotypes(freq, n, sdF, families)
  B <- genotypes(freq, n, sdF, families)
  a1 <- cbind(A$a1, B$a1); a2 <- cbind(A$a2, B$a2)
  gone <- matrix(stats::runif(length(a1)) < missing, nrow(a1))
  a1[gone] <- NA; a2[gone] <- NA
  samples <- c(paste0("a", seq_len(n)), paste0("b", seq_len(n)))
  dimnames(a1) <- dimnames(a2) <- list(NULL, samples)
  L <- nrow(freq); k <- ncol(freq)
  list(H = structure(list(A1 = a1, A2 = a2, locus = paste0("l", seq_len(L)),
                          locus_raw = paste0("l", seq_len(L)),
                          alleles = rep(list(c("AC", "CT", "GA", "TC")[seq_len(k)]), L),
                          n_alleles = rep(k, L), samples = samples), class = "raddiv_vcf"),
       pops = list(popA = samples[seq_len(n)], popB = samples[n + seq_len(n)]))
}

## ---------------------------------------------------------------------------
## An individual jackknife for Ar and privAr
## ---------------------------------------------------------------------------
## diversity_stats() no longer reports one (this simulation is why). It is
## kept here, written the way the package's jackknife for Ho/He/FIS works, so
## the evidence for that decision can be reproduced. Population A only.
jack_richness <- function(d, g, min_n = 2L) {
  H <- d$H; pops <- d$pops
  rows <- seq_len(nrow(H$A1))
  counts <- RADdiversity:::.pop_counts(H, pops, rows)
  n_typed <- RADdiversity:::.typed_by_pop(H, pops, rows)
  n_het <- RADdiversity:::.het_by_pop(H, pops, rows)
  p_drawn <- lapply(counts, RADdiversity:::.p_sampled_mat, g = g)
  absent <- RADdiversity:::.absent_elsewhere(p_drawn, 1)
  use <- n_typed[, 1] >= min_n
  stat <- function(C, n, h) {
    rs <- RADdiversity:::.one_pop_record_stats(C, n, h, absent, min_n, g, use = use)
    RADdiversity:::.diversity_from_sums(colSums(RADdiversity:::.diversity_pieces(rs)),
                                        "x")[1, c("Ar_x", "Pr_x")]
  }
  C <- counts[[1]]
  values <- vapply(pops[[1]], function(id) {
    a <- H$A1[rows, id]; b <- H$A2[rows, id]
    typed <- !is.na(a); tr <- which(typed)
    Cw <- C
    Cw[cbind(tr, a[tr])] <- Cw[cbind(tr, a[tr])] - 1L
    Cw[cbind(tr, b[tr])] <- Cw[cbind(tr, b[tr])] - 1L
    stat(Cw, n_typed[, 1] - typed, n_het[, 1] - (typed & a != b))
  }, numeric(2))
  n <- ncol(values)
  sqrt((n - 1) / n * rowSums((values - rowMeans(values))^2))
}

## ---------------------------------------------------------------------------
## One simulated study, and summaries over repeats
## ---------------------------------------------------------------------------
stat_names <- c("Ho", "He", "Fis", "Ar", "privAr")

## Population A's estimates and uncertainty from one dataset.
## boot = TRUE also runs the individual and two-way bootstraps.
one_study <- function(d, g, boot, nboot = 200L) {
  run <- function(mode, nb) suppressWarnings(
    diversity_stats(d$H, d$pops, g = g, nboot = nb, boot = mode, se_individuals = mode == "loci",
                    verbose = FALSE))
  r <- run("loci", if (boot) nboot else 0L)
  pp <- r$per_population[1, ]; rr <- r$richness[1, ]
  ind_rich <- jack_richness(d, g)
  out <- list(
    est    = c(pp$Ho, pp$He, pp$Fis, rr$Ar, rr$privAr),
    se     = c(pp$Ho_se, pp$He_se, pp$Fis_se, rr$Ar_se, rr$privAr_se),
    se_ind = c(pp$Ho_se_ind, pp$He_se_ind, pp$Fis_se_ind, ind_rich))
  if (boot) {
    ci <- function(res) {
      p <- res$per_population[1, ]; q <- res$richness[1, ]
      rbind(c(p$Ho_lo, p$He_lo, p$Fis_lo, q$Ar_lo, q$privAr_lo),
            c(p$Ho_hi, p$He_hi, p$Fis_hi, q$Ar_hi, q$privAr_hi))
    }
    out$ci_loci <- ci(r)
    out$ci_individuals <- ci(run("individuals", nboot))
    out$ci_both <- ci(run("both", nboot))
  }
  out
}

## Coverage (%) and SE ratio from a list of one_study() results.
summarise <- function(studies, truth, setting) {
  get <- function(name) t(vapply(studies, `[[`, numeric(5), name))
  est <- get("est")
  true_sd <- apply(est, 2, stats::sd)
  covers_se <- function(se) round(100 * colMeans(abs(sweep(est, 2, truth)) <= 1.96 * se, na.rm = TRUE))
  covers_ci <- function(name) {
    lo <- t(vapply(studies, function(s) s[[name]][1, ], numeric(5)))
    hi <- t(vapply(studies, function(s) s[[name]][2, ], numeric(5)))
    round(100 * colMeans(sweep(lo, 2, truth) <= 0 & sweep(hi, 2, truth) >= 0, na.rm = TRUE))
  }
  se <- get("se"); se_ind <- get("se_ind"); se_comb <- sqrt(se^2 + se_ind^2)
  rows <- list(
    "coverage %, locus SE"        = covers_se(se),
    "coverage %, individual SE"   = covers_se(se_ind),
    "coverage %, combined SE"     = covers_se(se_comb),
    "SE / true SE, locus"         = round(colMeans(se) / true_sd, 2),
    "SE / true SE, individual"    = round(colMeans(se_ind) / true_sd, 2),
    "SE / true SE, combined"      = round(colMeans(se_comb) / true_sd, 2))
  if (!is.null(studies[[1]]$ci_loci))
    rows <- c(rows, list("coverage %, locus bootstrap"      = covers_ci("ci_loci"),
                         "coverage %, individual bootstrap" = covers_ci("ci_individuals"),
                         "coverage %, two-way bootstrap"    = covers_ci("ci_both")))
  tab <- do.call(rbind, rows)
  colnames(tab) <- stat_names
  data.frame(setting = setting, measure = rownames(tab), tab, row.names = NULL, check.names = FALSE)
}

## ---------------------------------------------------------------------------
## Part "main": two truths x three levels of inbreeding variance
## ---------------------------------------------------------------------------
if ("main" %in% parts) {
  grid <- expand.grid(truth = c("genome", "these loci"), sdF = c(0, 0.1, 0.25),
                      stringsAsFactors = FALSE)
  L <- 1000; n <- 15; g <- 20
  run_setting <- function(i) {
    truth_kind <- grid$truth[i]; sdF <- grid$sdF[i]
    set.seed(1000 + i)
    fixed <- allele_freqs(L, "uniform")
    studies <- lapply(seq_len(reps), function(r) {
      freq <- if (truth_kind == "genome") allele_freqs(L, "uniform") else fixed
      one_study(make_data(freq, n, sdF), g, boot = TRUE)
    })
    He_t <- if (truth_kind == "genome") true_He("uniform") else mean(1 - rowSums(fixed^2))
    Fb <- mean_F(sdF)
    est <- t(vapply(studies, `[[`, numeric(5), "est"))
    truth <- c(He_t * (1 - Fb), He_t, Fb, mean(est[, 4]), mean(est[, 5]))
    summarise(studies, truth, sprintf("%s, SD(F) = %.2f", truth_kind, sdF))
  }
  res <- do.call(rbind, parallel::mclapply(seq_len(nrow(grid)), run_setting, mc.cores = cores))
  cat(sprintf("\n== main: 2 populations of %d, %d loci, g = %d, %d repeats per setting ==\n",
              n, L, g, reps))
  print(res, row.names = FALSE)
}

## ---------------------------------------------------------------------------
## Part "robust": genome truth under other conditions
## ---------------------------------------------------------------------------
if ("robust" %in% parts) {
  scenarios <- list(
    list(name = "baseline (n 15, 1,000 loci, SD(F) 0.1)", n = 15, L = 1000, sdF = 0.10, kind = "uniform", miss = 0,   fam = 0, g = 20),
    list(name = "many rare variants",                     n = 15, L = 1000, sdF = 0.10, kind = "rare",    miss = 0,   fam = 0, g = 20),
    list(name = "8 individuals",                          n = 8,  L = 1000, sdF = 0.10, kind = "uniform", miss = 0,   fam = 0, g = 10),
    list(name = "30 individuals, 3,000 loci",             n = 30, L = 3000, sdF = 0.10, kind = "uniform", miss = 0,   fam = 0, g = 40),
    list(name = "20% missing genotypes",                  n = 15, L = 1000, sdF = 0.10, kind = "uniform", miss = 0.2, fam = 0, g = 16),
    list(name = "4-allele haplotype loci",                n = 15, L = 1000, sdF = 0.10, kind = "haps",    miss = 0,   fam = 0, g = 20),
    list(name = "no inbreeding variance (SD(F) 0)",       n = 15, L = 1000, sdF = 0,    kind = "uniform", miss = 0,   fam = 0, g = 20),
    list(name = "5 full-sib families of 3",               n = 15, L = 1000, sdF = 0,    kind = "uniform", miss = 0,   fam = 5, g = 20))
  run_scenario <- function(i) {
    s <- scenarios[[i]]
    set.seed(2000 + i)
    studies <- lapply(seq_len(reps), function(r)
      one_study(make_data(allele_freqs(s$L, s$kind), s$n, s$sdF, s$fam, s$miss), s$g, boot = FALSE))
    Fb <- if (s$fam > 0) 0 else mean_F(s$sdF)
    He_t <- true_He(s$kind)
    est <- t(vapply(studies, `[[`, numeric(5), "est"))
    truth <- c(He_t * (1 - Fb), He_t, Fb, mean(est[, 4]), mean(est[, 5]))
    tab <- summarise(studies, truth, s$name)
    ## Bias of the estimate itself, in true SEs (relatives bias He and FIS).
    bias <- round((colMeans(est) - truth) / apply(est, 2, stats::sd), 1)
    rbind(tab, data.frame(setting = s$name, measure = "bias of estimate, in true SEs",
                          t(stats::setNames(bias, stat_names)), check.names = FALSE))
  }
  res <- do.call(rbind, parallel::mclapply(seq_along(scenarios), run_scenario, mc.cores = cores))
  cat(sprintf("\n== robust: genome truth, %d repeats per scenario ==\n", reps))
  print(res, row.names = FALSE)
}

## ---------------------------------------------------------------------------
## Part "boundary": records with only 2 genotyped individuals
## ---------------------------------------------------------------------------
## Leaving out one of 2 genotyped individuals leaves He undefined, so the
## record drops out of that jackknife replicate. Same loci in every repeat
## ("these loci" truth), so the individual SE alone is the right comparison.
## The share of such records is varied through the missing-genotype rate.
if ("boundary" %in% parts) {
  L <- 1000; n <- 6; g <- 4
  missing_rates <- c(0, 0.2, 0.3, 0.4, 0.5, 0.6)
  run_rate <- function(i) {
    set.seed(3000 + i)
    fixed <- allele_freqs(L, "uniform")
    out <- t(replicate(reps, {
      d <- make_data(fixed, n, 0.1, missing = missing_rates[i])
      r <- suppressWarnings(diversity_stats(d$H, d$pops, g = g, nboot = 0, se_individuals = TRUE,
                                            verbose = FALSE))
      c(r$per_population$Ho[1], r$per_population$He[1], r$per_population$Fis[1],
        r$per_population$Ho_se_ind[1], r$per_population$He_se_ind[1],
        r$per_population$Fis_se_ind[1], r$settings$jackknife_boundary[1])
    }))
    ratio <- round(colMeans(out[, 4:6]) / apply(out[, 1:3], 2, stats::sd), 2)
    data.frame(missing = missing_rates[i], share_with_2 = round(mean(out[, 7]), 3),
               Ho = ratio[1], He = ratio[2], Fis = ratio[3])
  }
  res <- do.call(rbind, parallel::mclapply(seq_along(missing_rates), run_rate, mc.cores = cores))
  cat(sprintf("\n== boundary: individual SE / true SE, 2 populations of %d, %d loci, %d repeats ==\n",
              n, L, reps))
  print(res, row.names = FALSE)
}
