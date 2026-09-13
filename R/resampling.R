###############################################################################
#
#  R/resampling.R -- the block jackknife and block bootstrap over RAD loci,
#  shared by diversity_stats(), differentiation_stats() and pi_allsites().
#
#  THE IDEA. Every statistic these functions report is built from SUMS over
#  records: FIS = 1 - sum(Ho)/sum(Hs), FST = sum(a)/sum(a+b+c), and a mean
#  over loci is a sum divided by a count. Resampling whole RAD loci therefore
#  only changes how many times each locus's sums are added in. So each caller
#  first sums its per-record pieces once per RAD locus, giving a matrix S
#  (one row per locus, one column per piece), plus a function f() that turns
#  summed pieces into statistics. Then:
#
#    point estimate              f(colSums(S))
#    jackknife, leave out b      f(colSums(S) - S[b, ])       -- all b at once
#    bootstrap replicate         f(counts %*% S), counts[b] = times locus b drawn
#
#  This gives the same numbers as re-summing the resampled records in every
#  replicate (up to rounding in the last few digits), and costs
#  (number of loci) x (number of pieces) instead of re-reading every record.
#
#  This shortcut works only when whole LOCI are resampled. When individuals
#  are resampled, each record's own values change and must be recomputed
#  (see .boot_individuals() and .jack_individuals() in
#  R/diversity_internals.R).
#
###############################################################################

## Not exported. Delete-one-block jackknife standard error.
##   S  one row per RAD locus, one column per summed piece
##   f  maps a matrix of summed pieces (one row per estimate) to a matrix of
##      statistics (one row per estimate, named columns)
## Returns a named vector of standard errors:
##   SE = sqrt( (nL - 1)/nL * sum_b (theta_(-b) - mean_b theta_(-b))^2 )
.jack_block_sums <- function(S, f) {
  n_loci <- nrow(S)
  ## Row b of `leave_one_out` = the totals with locus b removed.
  leave_one_out <- sweep(-S, 2L, colSums(S), "+")
  estimates <- f(leave_one_out)
  deviations <- sweep(estimates, 2L, colMeans(estimates))
  sqrt(((n_loci - 1) / n_loci) * colSums(deviations^2))
}

## Not exported. Runs `nboot` bootstrap replicates in batches.
##   n_units   how many units (RAD loci, or individuals) are resampled
##   f         a function of `counts`, an n_units x m matrix whose column i
##             says how many times each unit was drawn in replicate i; it
##             returns an m x (statistics) matrix
##   cells_per_replicate, max_cells
##             batches are sized so that one batch needs at most `max_cells`
##             matrix cells, keeping memory use bounded
## Each replicate uses exactly one sample.int(n_units, n_units, replace = TRUE)
## call, so results match drawing the replicates one at a time.
## Returns an nboot x (statistics) matrix.
.boot_in_batches <- function(n_units, nboot, f, cells_per_replicate = n_units,
                             max_cells = 5e6) {
  batch_size <- max(1L, min(nboot, floor(max_cells / max(1, cells_per_replicate))))
  results <- vector("list", ceiling(nboot / batch_size))
  done <- 0L
  batch <- 0L
  while (done < nboot) {
    m <- min(batch_size, nboot - done)
    counts <- vapply(seq_len(m), function(i)
      as.numeric(tabulate(sample.int(n_units, n_units, replace = TRUE), nbins = n_units)),
      numeric(n_units))
    counts <- matrix(counts, nrow = n_units, ncol = m)
    batch <- batch + 1L
    results[[batch]] <- f(counts)
    done <- done + m
  }
  do.call(rbind, results)
}

## Not exported. Block bootstrap over RAD loci: each replicate draws the
## loci (rows of S) with replacement and f() turns that replicate's summed
## pieces into statistics. Returns an nboot x (statistics) matrix.
.boot_block_sums <- function(S, nboot, f, max_cells = 5e6) {
  .boot_in_batches(nrow(S), nboot, function(counts) f(crossprod(counts, S)),
                   max_cells = max_cells)
}

## Not exported. 2.5% and 97.5% quantiles of each column of a bootstrap
## matrix (replicates x statistics): a statistics x 2 matrix (lo, hi). With
## no replicates, every interval is NA.
.percentile_ci <- function(boot, stat_names) {
  if (is.null(boot))
    return(matrix(NA_real_, length(stat_names), 2L, dimnames = list(stat_names, c("lo", "hi"))))
  ci <- t(apply(boot, 2L, stats::quantile, c(0.025, 0.975), na.rm = TRUE, names = FALSE))
  dimnames(ci) <- list(colnames(boot), c("lo", "hi"))
  ci
}
