###############################################################################
#
#  inst/sims/jackknife_richness.R -- are the delete-one-individual standard
#  errors of allelic richness (Ar) and private allelic richness (privAr) the
#  right size?
#
#  Two populations of 10 individuals are sampled again and again from one
#  fixed set of 3,000 RAD loci with 4 alleles each (random frequencies, the
#  same in both populations). Genotypes are then set missing at random at a
#  given rate. For each setting, diversity_stats(se_individuals = TRUE) is run
#  on every replicate, and the average `_se_ind` of population A is divided by
#  the standard deviation of its estimate across replicates -- the true
#  standard error. A ratio of 1 is right; above 1 the reported SE is too
#  large, below 1 too small.
#
#  `boundary` is the share of population A's records with a defined Ar that
#  have fewer than g + 2 gene copies: there, deleting one individual leaves
#  fewer than g copies, so the record drops out of some jackknife replicates
#  and not others. diversity_stats() warns when this share is large (see
#  .jackknife_boundary() in R/diversity_internals.R).
#
#  Quoted in vignette("rationale"), section 4 ("Individual standard errors
#  for Ar and privAr").
#
#  Run from the package source:   Rscript inst/sims/jackknife_richness.R [reps]
#  (default 200 replicates per setting; settings run in parallel on up to 4
#  cores; about 10 minutes). Needs the package installed, or run under
#  devtools::load_all().
#
###############################################################################

if (!requireNamespace("RADdiversity", quietly = TRUE)) stop("Install RADdiversity first.")
suppressMessages(library(RADdiversity))
args <- commandArgs(trailingOnly = TRUE)
reps <- if (length(args)) as.integer(args[1]) else 200L

L <- 3000; n <- 10; k <- 4
set.seed(42)
freq <- t(replicate(L, { x <- stats::rgamma(k, 0.6); x / sum(x) }))
cum <- t(apply(freq, 1, cumsum))
samp <- c(paste0("a", seq_len(n)), paste0("b", seq_len(n)))
pops <- list(popA = samp[seq_len(n)], popB = samp[n + seq_len(n)])

draw_H <- function(miss) {
  pick <- function() {
    u <- matrix(stats::runif(L * 2 * n), L)
    1L + (u > cum[, 1]) + (u > cum[, 2]) + (u > cum[, 3])
  }
  A1 <- pick(); A2 <- pick()
  missing <- matrix(stats::runif(L * 2 * n) < miss, L)
  A1[missing] <- NA; A2[missing] <- NA
  dimnames(A1) <- dimnames(A2) <- list(NULL, samp)
  structure(list(A1 = A1, A2 = A2, locus = paste0("r", seq_len(L)),
                 locus_raw = paste0("L", seq_len(L)),
                 alleles = rep(list(c("AC", "CT", "GA", "TC")), L),
                 n_alleles = rep(4L, L), samples = samp), class = "raddiv_vcf")
}

settings <- data.frame(missing = c(0, 0, 0.01, 0.03, 0.10, 0.10, 0.10, 0.10, 0.10),
                       g       = c(10, 18, 18,   18,   10,   12,   14,   16,   18))

one_setting <- function(i) {
  miss <- settings$missing[i]; g <- settings$g[i]
  set.seed(1000 + i)
  out <- t(replicate(reps, {
    H <- draw_H(miss)
    copies <- 2 * rowSums(!is.na(H$A1[, pops$popA]))
    defined <- copies >= g
    r <- suppressWarnings(diversity_stats(H, pops, g = g, nboot = 0, se_individuals = TRUE,
                                          verbose = FALSE))
    c(Ar = r$richness$Ar[1], Ar_se_ind = r$richness$Ar_se_ind[1],
      privAr = r$richness$privAr[1], privAr_se_ind = r$richness$privAr_se_ind[1],
      He = r$per_population$He[1], He_se_ind = r$per_population$He_se_ind[1],
      Fis = r$per_population$Fis[1], Fis_se_ind = r$per_population$Fis_se_ind[1],
      boundary = mean(copies[defined] < g + 2))
  }))
  ratio <- function(s) mean(out[, paste0(s, "_se_ind")]) / stats::sd(out[, s])
  data.frame(missing = miss, g = g, boundary = round(mean(out[, "boundary"]), 3),
             Ar = round(ratio("Ar"), 2), privAr = round(ratio("privAr"), 2),
             He = round(ratio("He"), 2), Fis = round(ratio("Fis"), 2))
}

cores <- if (.Platform$OS.type == "unix") min(4L, parallel::detectCores()) else 1L
rows <- parallel::mclapply(seq_len(nrow(settings)), one_setting, mc.cores = cores)
result <- do.call(rbind, rows)
cat(sprintf("Reported individual SE / true SE, population A (%d replicates per row):\n", reps))
print(result, row.names = FALSE)
