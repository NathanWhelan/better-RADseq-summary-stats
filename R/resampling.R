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
