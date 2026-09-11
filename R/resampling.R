###############################################################################
#
#  R/resampling.R -- the block jackknife and block bootstrap over RAD loci,
#  shared by diversity_stats() and differentiation_stats().
#
#  THE IDEA. Every statistic these functions report over loci is built from
#  SUMS over records: FIS = 1 - sum(Ho)/sum(Hs), FST = sum(a)/sum(a+b+c), a
#  mean over loci is a sum divided by a count, and so on. Resampling whole
#  RAD loci therefore only ever changes how many times each locus's sums are
#  added in. So each caller first sums its per-record pieces once per RAD
#  locus, giving a matrix S (one row per locus, one column per piece), plus a
#  function f() that turns a row of summed pieces into the statistics. Then:
#
#    point estimate            f(colSums(S))
#    jackknife, drop locus b   f(colSums(S) - S[b, ])        -- all b at once
#    bootstrap replicate       f(counts %*% S), counts[b] = times locus b drawn
#
#  This gives the same numbers as re-summing the resampled records every
#  replicate (up to floating-point rounding in the last digits), without
#  touching the per-record matrices again. The old record-by-record jackknife
#  cost (number of loci) x (number of records); this costs (number of loci)
#  x (number of pieces).
#
#  ONLY for resampling whole loci. When individuals are resampled, the
#  per-record values themselves change every replicate and must be
#  recomputed (diversity_stats()'s stat_from_resampled()).
#
###############################################################################

## Not exported. Delete-one-block jackknife standard error. `S` has one row
## per RAD locus; f() maps a matrix of summed pieces (one row per estimate)
## to a matrix of statistics (one row per estimate, named columns).
##   SE = sqrt( (nL - 1)/nL * sum_b (theta_(-b) - mean_b theta_(-b))^2 )
.jack_block_sums <- function(S, f) {
  nL <- nrow(S)
  jk <- f(sweep(-S, 2L, colSums(S), "+"))     # row b = estimate without locus b
  theta_bar <- colMeans(jk)
  sqrt(((nL - 1) / nL) * colSums(sweep(jk, 2L, theta_bar)^2))
}

## Not exported. Delete-one-INDIVIDUAL jackknife standard errors for each
## population's Ho, He, FIS, Ar and privAr (diversity_stats(se_individuals =
## TRUE)): the uncertainty from which individuals were sampled, which the
## loci jackknife and bootstrap leave out. For each individual of population
## p, its two allele copies and its heterozygosity are removed from p's
## per-locus counts -- only p's own statistics change -- and p's statistics
## are recomputed with the point estimate's formulas (the min_n and >= g
## gene-copy rules included). Nobody is ever duplicated, so this carries none
## of the Nei-Chesser / rarefaction bias that rules out a bootstrap over
## individuals ("Bootstrap mode" in vignette("rationale")).
##   SE = sqrt( (n - 1)/n * sum_i (theta_(-i) - mean_i theta_(-i))^2 )
## `cmats`: per-record allele-count tables (populations x alleles), for the
## records `idx`; `n_typed`: typed individuals per record (rows = idx) and
## population. Returns list(se = named SEs, full = the same statistics
## recomputed from everyone, which must equal the point estimates).
.jack_individuals <- function(H, pops, idx, cmats, n_typed, min_n, g) {
  r <- length(pops); L <- length(idx)
  kmax <- max(vapply(cmats, ncol, integer(1)))
  ## Population p's counts as one records x alleles matrix (zero-padded).
  counts_of <- function(p) matrix(vapply(cmats, function(m) {
    v <- numeric(kmax); v[seq_len(ncol(m))] <- m[p, ]; v }, numeric(kmax)),
    nrow = L, byrow = TRUE)
  cnt <- lapply(seq_len(r), counts_of)
  ## Pr(each allele appears in g copies drawn from N), NA where N < g -- the
  ## same formula as p_sampled(), for every record at once.
  psamp <- function(C, N) {
    out <- 1 - exp(lchoose(N - C, g) - lchoose(N, g))
    out[N < g, ] <- NA_real_
    out
  }
  ps_all <- lapply(seq_len(r), function(p) psamp(cnt[[p]], 2 * n_typed[, p]))
  het <- lapply(pops, function(ids)
    rowSums(H$A1[idx, ids, drop = FALSE] != H$A2[idx, ids, drop = FALSE], na.rm = TRUE))
  ## One population's Ho, He, FIS, Ar, privAr from its counts (C), typed
  ## individuals (n) and heterozygote counts per record; Q = the other
  ## populations' probability of NOT sampling each allele (for privAr).
  pop_stats <- function(C, n, hetc, Q) {
    ho <- hetc / n
    hs <- hs_nei_chesser(rowSums((C / (2 * n))^2), ho, n)
    ok <- n >= min_n & is.finite(hs) & is.finite(ho)
    ps <- psamp(C, 2 * n)
    ar <- rowSums(ps); pr <- rowSums(ps * Q)
    ar[n < min_n] <- NA_real_; pr[n < min_n] <- NA_real_
    shs <- sum(hs[ok])
    c(Ho  = if (any(ok)) sum(ho[ok]) / sum(ok) else NA_real_,
      He  = if (any(ok)) shs / sum(ok) else NA_real_,
      Fis = if (shs > 1e-12) 1 - sum(ho[ok]) / shs else NA_real_,
      Ar  = if (any(!is.na(ar))) mean(ar, na.rm = TRUE) else NA_real_,
      Pr  = if (any(!is.na(pr))) mean(pr, na.rm = TRUE) else NA_real_)
  }
  se <- full <- numeric(0)
  for (p in seq_len(r)) {
    Q <- Reduce(`*`, lapply(setdiff(seq_len(r), p), function(k) 1 - ps_all[[k]]))
    C <- cnt[[p]]; n <- n_typed[, p]; hc <- het[[p]]
    jk <- vapply(pops[[p]], function(id) {
      a <- H$A1[idx, id]; b <- H$A2[idx, id]
      ty <- which(!is.na(a))
      Cn <- C
      ## Two separate steps, so a homozygote loses both of its copies (one
      ## vectorised assignment with a repeated index would subtract once).
      Cn[cbind(ty, a[ty])] <- Cn[cbind(ty, a[ty])] - 1
      Cn[cbind(ty, b[ty])] <- Cn[cbind(ty, b[ty])] - 1
      typed <- !is.na(a)
      pop_stats(Cn, n - typed, hc - (typed & a != b), Q)
    }, numeric(5))
    nj <- ncol(jk)
    nm <- paste0(c("Ho_", "He_", "Fis_", "Ar_", "Pr_"), names(pops)[p])
    se <- c(se, stats::setNames(sqrt((nj - 1) / nj * rowSums((jk - rowMeans(jk))^2)), nm))
    full <- c(full, stats::setNames(pop_stats(C, n, hc, Q), nm))
  }
  list(se = se, full = full)
}

## Not exported. Block bootstrap: each replicate draws nL loci with
## replacement, using exactly one sample.int(nL, nL, replace = TRUE) call per
## replicate (the same random numbers as resampling the records directly),
## and f() turns that replicate's summed pieces into statistics. Replicates
## are processed in batches so the nL x batch matrix of draw counts stays
## under `max_cells` entries. Returns an nboot x (statistics) matrix.
.boot_block_sums <- function(S, nboot, f, max_cells = 5e6) {
  nL <- nrow(S)
  batch <- max(1L, min(nboot, floor(max_cells / nL)))
  out <- vector("list", ceiling(nboot / batch))
  done <- 0L; k <- 0L
  while (done < nboot) {
    m <- min(batch, nboot - done)
    counts <- vapply(seq_len(m), function(i)
      as.numeric(tabulate(sample.int(nL, nL, replace = TRUE), nbins = nL)),
      numeric(nL))
    k <- k + 1L
    out[[k]] <- f(crossprod(counts, S))       # m replicates x statistics
    done <- done + m
  }
  do.call(rbind, out)
}
