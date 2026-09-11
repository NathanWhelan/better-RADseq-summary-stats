## The block-sum jackknife and bootstrap (R/resampling.R) must give the same
## answers as resampling the records themselves. A ratio of sums over 40
## records grouped into 15 blocks stands in for FIS / FST.
ratio_of_sums <- function(tot) {
  tot <- matrix(tot, ncol = 2)
  cbind(ratio = tot[, 1] / tot[, 2])
}
set.seed(1)
x <- runif(40); y <- runif(40) + 1
blk <- rep(1:15, length.out = 40)
S <- rowsum(cbind(x, y), blk)

test_that(".jack_block_sums() equals a brute-force delete-one-block jackknife", {
  drop_one <- vapply(1:15, function(b) sum(x[blk != b]) / sum(y[blk != b]), numeric(1))
  se_brute <- sqrt(14 / 15 * sum((drop_one - mean(drop_one))^2))
  expect_equal(unname(RADdiversity:::.jack_block_sums(S, ratio_of_sums)), se_brute,
               tolerance = 1e-12)
})

test_that(".boot_block_sums() reproduces resampling whole blocks with the same random numbers", {
  rows_of <- split(seq_along(x), blk)
  set.seed(2)
  brute <- replicate(50, {
    i <- unlist(rows_of[sample.int(15, 15, replace = TRUE)])
    sum(x[i]) / sum(y[i])
  })
  set.seed(2)
  ## max_cells = 30 forces 25 batches of 2 replicates, so the batching is
  ## exercised too.
  fast <- RADdiversity:::.boot_block_sums(S, 50, ratio_of_sums, max_cells = 30)
  expect_equal(dim(fast), c(50L, 1L))
  expect_equal(as.vector(fast), brute, tolerance = 1e-12)
})
