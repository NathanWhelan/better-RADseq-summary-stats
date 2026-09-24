## read_distances() and isolation_by_distance(): the two distance layouts,
## matching populations by name, the Mantel test checked against
## vegan::mantel(), the regression, and the summary.

ex <- function(name) system.file("extdata", name, package = "RADdiversity")

write_lines <- function(lines, ext = ".csv") {
  path <- tempfile(fileext = ext)
  writeLines(lines, path)
  path
}

## Three populations: a-b 15, a-c 27, b-c 12.
km3 <- matrix(c(0, 15, 27, 15, 0, 12, 27, 12, 0), 3,
              dimnames = list(c("a", "b", "c"), c("a", "b", "c")))

## n populations along a line, with FST rising with distance plus noise.
toy <- function(n, noise = 0.01, seed = 1) {
  set.seed(seed)
  pos <- cumsum(stats::runif(n, 5, 20))
  pops <- paste0("p", seq_len(n))
  km <- abs(outer(pos, pos, "-"))
  e <- matrix(stats::runif(n * n, 0, noise), n)
  fst <- km / 2000 + (e + t(e)) / 2
  diag(fst) <- 0
  dimnames(km) <- dimnames(fst) <- list(pops, pops)
  list(fst = fst, km = km)
}

## ---------------------------------------------------------------------------
## read_distances()
## ---------------------------------------------------------------------------

test_that("read_distances() reads a full, an upper-half and a lower-half square matrix", {
  full <- c(",a,b,c", "a,0,15,27", "b,15,0,12", "c,27,12,0")
  upper <- c(",a,b,c", "a,0,15,27", "b,,0,12", "c,,,0")
  lower <- c(",a,b,c", "a,0,,", "b,15,0,", "c,27,12,0")
  for (lines in list(full, upper, lower))
    expect_equal(read_distances(write_lines(lines), verbose = FALSE), km3)
  expect_message(read_distances(write_lines(upper)), "lower half was blank")
  expect_message(read_distances(write_lines(lower)), "upper half was blank")
})

test_that("read_distances() matches columns to rows by name, with or without a top-left cell", {
  shuffled <- c(",c,a,b", "b,12,15,0", "a,27,0,15", "c,0,27,12")
  m <- read_distances(write_lines(shuffled), verbose = FALSE)
  expect_equal(m[c("a", "b", "c"), c("a", "b", "c")], km3)
  no_corner <- c("a,b,c", "a,0,15,27", "b,15,0,12", "c,27,12,0")
  expect_equal(read_distances(write_lines(no_corner), verbose = FALSE), km3)
  labelled_corner <- c("pop,a,b,c", "a,0,15,27", "b,15,0,12", "c,27,12,0")
  expect_equal(read_distances(write_lines(labelled_corner), verbose = FALSE), km3)
})

test_that("read_distances() reads a list of pairs, with or without a header, in either order", {
  with_header <- c("pop1,pop2,distance", "a,b,15", "a,c,27", "b,c,12")
  expect_equal(read_distances(write_lines(with_header), verbose = FALSE), km3)
  no_header <- c("a,b,15", "c,a,27", "c,b,12")
  expect_equal(read_distances(write_lines(no_header), verbose = FALSE), km3)
  ## a pair listed twice (either order) with the same distance is fine
  expect_equal(read_distances(write_lines(c(with_header, "b,a,15")), verbose = FALSE), km3)
  ## a pair with no distance is NA, with a message
  gap <- c("pop1,pop2,distance", "a,b,15", "b,c,12")
  messages <- capture_messages(m <- read_distances(write_lines(gap)))
  expect_true(any(grepl("1 pair(s) have no distance: a - c", messages, fixed = TRUE)))
  expect_true(is.na(m["a", "c"]))
})

test_that("read_distances() reads TAB files, a byte order mark, Windows line ends and numeric names", {
  tab <- write_lines(c("\ta\tb\tc", "a\t0\t15\t27", "b\t\t0\t12", "c\t\t\t0"), ".tsv")
  expect_equal(read_distances(tab, verbose = FALSE), km3)
  ## As Excel writes a "CSV UTF-8" file on Windows
  excel <- tempfile(fileext = ".csv")
  con <- file(excel, "wb")
  writeBin(c(as.raw(c(0xef, 0xbb, 0xbf)),
             charToRaw(",a,b,c\r\na,0,15,27\r\nb,15,0,12\r\nc,27,12,0\r\n")), con)
  close(con)
  expect_equal(read_distances(excel, verbose = FALSE), km3)
  ## Populations named by numbers, in both layouts
  m <- read_distances(write_lines(c(",1,2,3", "1,0,15,27", "2,,0,12", "3,,,0")), verbose = FALSE)
  expect_identical(rownames(m), c("1", "2", "3"))
  expect_equal(unname(m), unname(km3))
  m <- read_distances(write_lines(c("1,2,15", "1,3,27", "2,3,12")), verbose = FALSE)
  expect_equal(unname(m), unname(km3))
})

test_that("read_distances() takes a matrix, a dist object or a data frame, at full precision", {
  expect_equal(read_distances(km3, verbose = FALSE), km3)
  expect_equal(read_distances(stats::as.dist(km3), verbose = FALSE), km3)
  expect_equal(read_distances(as.data.frame(km3), verbose = FALSE), km3)
  pairs <- data.frame(pop1 = c("a", "a", "b"), pop2 = c("b", "c", "c"), distance = c(15, 27, 12))
  expect_equal(read_distances(pairs, verbose = FALSE), km3)
  ## as read.csv() reads a square file without row.names = 1
  df <- utils::read.csv(write_lines(c(",a,b,c", "a,0,15,27", "b,,0,12", "c,,,0")),
                        check.names = FALSE)
  expect_equal(read_distances(df, verbose = FALSE), km3)
  thirds <- km3 + 1 / 3
  diag(thirds) <- 0
  expect_identical(read_distances(thirds, verbose = FALSE), thirds)
})

test_that("read_distances() stops on mistakes with a message that says what is wrong", {
  read <- function(lines) read_distances(write_lines(lines), verbose = FALSE)
  expect_error(read(c(",a,b,c", "a,0,15,27", "b,16,0,12", "c,27,12,0")),
               "a - b is 15 in one half of the matrix and 16 in the other")
  expect_error(read(c("a,b,15", "b,a,16", "a,c,1", "b,c,2")), "listed twice, with 15 and 16")
  expect_error(read(c(",a,b,c", "a,0,-5,27", "b,,0,12", "c,,,0")), "negative")
  expect_error(read(c(",a,b,c", "a,0,15 km,27", "b,,0,12", "c,,,0")),
               "from a to b is \"15 km\", which is not a number")
  expect_error(read(c(",a,b,c", "a,3,15,27", "b,,0,12", "c,,,0")), "from a to itself is 3")
  expect_error(read(c("pop1,pop2,distance", "a,b,", "a,c,27")), "row 2 \\(a, b\\) has a blank")
  expect_error(read(c("a,b,c,d", "1,2,3,4")), "Use one of these two layouts")
  expect_error(read_distances(tempfile()), "Distance file not found")
  expect_error(read_distances(list(1)), "must be a file path")
  err <- tryCatch(read(c("a,b,c,d", "1,2,3,4")), error = function(e) e)
  expect_null(conditionCall(err))
})

## ---------------------------------------------------------------------------
## isolation_by_distance(): the test
## ---------------------------------------------------------------------------

test_that("every ordering of n populations is listed once", {
  o <- RADdiversity:::.all_orderings(4)
  expect_equal(dim(o), c(24L, 4L))
  expect_equal(nrow(unique(o)), 24L)
  expect_true(all(apply(o, 1, function(r) setequal(r, 1:4))))
})

test_that("r and the exact p-value equal vegan::mantel()'s, Pearson and Spearman", {
  skip_if_not_installed("vegan")
  for (n in c(5, 6)) for (noise in c(0.01, 0.06)) {
    t <- toy(n, noise)
    lin <- t$fst / (1 - t$fst)
    for (method in c("pearson", "spearman")) {
      res <- isolation_by_distance(t$fst, t$km, habitat = "1D", method = method, verbose = FALSE)
      ## vegan tries every ordering here too (complete enumeration)
      v <- suppressMessages(vegan::mantel(stats::as.dist(lin), stats::as.dist(t$km),
                                          method = method, permutations = 9999))
      expect_equal(res$result$mantel_r, unname(v$statistic))
      expect_equal(res$result$p_value, v$signif)
      expect_true(res$result$exact)
      expect_equal(res$result$orderings, factorial(n))
    }
  }
})

test_that("with 8 populations the p from random orderings is close to the exact p", {
  t <- toy(8, noise = 0.1, seed = 4)            # exact p about 0.16
  res <- isolation_by_distance(t$fst, t$km, habitat = "1D", seed = 11, verbose = FALSE)
  expect_false(res$result$exact)
  expect_equal(res$result$orderings, 9999)
  ## the exact p, from all 40,320 orderings
  lin <- t$fst / (1 - t$fst)
  low <- lower.tri(lin)
  x <- t$km[low]
  r_all <- apply(RADdiversity:::.all_orderings(8), 1, function(o) stats::cor(lin[o, o][low], x))
  exact_p <- mean(r_all >= res$result$mantel_r - 1e-8)
  expect_gt(exact_p, 0.1)                        # far from the floor, so this is a real check
  expect_lt(abs(res$result$p_value - exact_p), 0.02)   # about 5 SE of a 9,999-ordering p
})

test_that("FST is used as FST/(1 - FST), D as it is, and a 2-D habitat uses the natural log", {
  t <- toy(5)
  fst <- isolation_by_distance(t$fst, t$km, habitat = "1D", verbose = FALSE)
  expect_equal(fst$pairs$FST_linearized, fst$pairs$FST / (1 - fst$pairs$FST))
  fit <- stats::lm(FST_linearized ~ distance, data = fst$pairs)
  expect_equal(c(fst$result$intercept, fst$result$slope), unname(stats::coef(fit)))
  expect_equal(fst$result$mantel_r, stats::cor(fst$pairs$FST_linearized, fst$pairs$distance))

  d <- isolation_by_distance(t$fst, t$km, habitat = "2D", stat = "D", verbose = FALSE)
  expect_equal(d$pairs$log_distance, log(d$pairs$distance))
  expect_equal(d$pairs$D, fst$pairs$FST)              # the same numbers, used as they are
  expect_null(d$pairs$FST_linearized)
  fit <- stats::lm(D ~ log_distance, data = d$pairs)
  expect_equal(c(d$result$intercept, d$result$slope), unname(stats::coef(fit)))
  expect_equal(d$result$mantel_r, stats::cor(d$pairs$D, d$pairs$log_distance))
  expect_identical(d$result$geographic, "ln(distance)")
})

test_that("populations are matched by name, never by order", {
  t <- toy(6, noise = 0.03)
  a <- isolation_by_distance(t$fst, t$km, habitat = "1D", verbose = FALSE)
  o <- c(4, 2, 6, 1, 5, 3)
  b <- isolation_by_distance(t$fst, t$km[o, o], habitat = "1D", verbose = FALSE)
  expect_equal(a$result, b$result)
  expect_equal(a$pairs, b$pairs)

  extra <- matrix(1, 7, 7, dimnames = list(c(rownames(t$km), "x"), c(rownames(t$km), "x")))
  extra[1:6, 1:6] <- t$km
  diag(extra) <- 0
  messages <- capture_messages(res <- isolation_by_distance(t$fst, extra, habitat = "1D"))
  expect_true(any(grepl("(left out): x", messages, fixed = TRUE)))
  expect_equal(res$result, a$result)
  expect_identical(res$settings$dropped_from_distances, "x")

  expect_error(isolation_by_distance(t$fst, t$km[-2, -2], habitat = "1D", verbose = FALSE),
               "in the genetic data but not in the distances: p2")
  renamed <- t$km
  rownames(renamed)[3] <- colnames(renamed)[3] <- "P 3"
  expect_error(isolation_by_distance(t$fst, renamed, habitat = "1D", verbose = FALSE),
               "p3 \\(did you mean \"P 3\"\\?\\)")
})

test_that("exclude leaves populations out, and an unknown name stops", {
  t <- toy(6, noise = 0.03)
  a <- isolation_by_distance(t$fst, t$km, habitat = "1D", exclude = "p2", verbose = FALSE)
  b <- isolation_by_distance(t$fst[-2, -2], t$km[-2, -2], habitat = "1D", verbose = FALSE)
  expect_equal(a$result, b$result)
  expect_identical(a$settings$excluded, "p2")
  expect_false(any(c(a$pairs$pop1, a$pairs$pop2) == "p2"))
  expect_error(isolation_by_distance(t$fst, t$km, habitat = "1D", exclude = "p9"),
               "not in the genetic data: p9")
  expect_error(isolation_by_distance(t$fst, t$km, habitat = "1D", exclude = 2), "`exclude` must be")
})

test_that("a missing habitat and impossible requests stop with a clear message", {
  t <- toy(5)
  expect_error(isolation_by_distance(t$fst, t$km), "habitat = \"1D\" for a river")
  expect_error(isolation_by_distance(t$fst, t$km), "habitat = \"2D\" for an area")
  expect_error(isolation_by_distance(t$fst, t$km, habitat = "3D"), "`habitat` must be one of")
  expect_error(isolation_by_distance(t$fst, t$km, habitat = "1D", stat = "beta"),
               "`stat` must be one of")
  expect_error(isolation_by_distance(t$fst, t$km, habitat = "1D", method = "kendall"),
               "`method` must be one of")
  expect_error(isolation_by_distance(t$fst, t$km, habitat = "1D", nperm = 0), "`nperm`")

  zero <- t$km
  zero["p1", "p2"] <- zero["p2", "p1"] <- 0
  expect_error(isolation_by_distance(t$fst, zero, habitat = "2D", verbose = FALSE),
               "distance 0: p1 - p2")
  expect_silent(isolation_by_distance(t$fst, zero, habitat = "1D", verbose = FALSE))

  no_fst <- t$fst
  no_fst["p1", "p3"] <- no_fst["p3", "p1"] <- NA
  expect_error(isolation_by_distance(no_fst, t$km, habitat = "1D", verbose = FALSE),
               "No FST for: p1 - p3")
  no_km <- t$km
  no_km["p1", "p4"] <- no_km["p4", "p1"] <- NA
  expect_error(isolation_by_distance(t$fst, no_km, habitat = "1D", verbose = FALSE),
               "No distance for: p1 - p4")
  one <- t$fst
  one["p1", "p2"] <- one["p2", "p1"] <- 1
  expect_error(isolation_by_distance(one, t$km, habitat = "1D", verbose = FALSE),
               "FST is 1 or more for: p1 - p2")
  expect_error(isolation_by_distance(t$fst[1:2, 1:2], t$km, habitat = "1D", verbose = FALSE),
               "at least 3 populations")
  expect_error(isolation_by_distance("fst.csv", t$km, habitat = "1D"),
               "must be the result of differentiation_stats\\(\\)")
  lopsided <- t$fst
  lopsided["p1", "p2"] <- 0.5
  expect_error(isolation_by_distance(lopsided, t$km, habitat = "1D"), "not symmetric: p1 - p2")
})

test_that("seed makes random orderings reproducible without changing the session's random numbers", {
  t <- toy(8, noise = 0.06, seed = 4)
  set.seed(5)
  expected_next <- stats::runif(1)
  set.seed(5)
  a <- isolation_by_distance(t$fst, t$km, habitat = "1D", nperm = 999, seed = 1, verbose = FALSE)
  expect_equal(stats::runif(1), expected_next)
  b <- isolation_by_distance(t$fst, t$km, habitat = "1D", nperm = 999, seed = 1, verbose = FALSE)
  expect_identical(a$result$p_value, b$result$p_value)
})

test_that("the example data show isolation by distance along the river, for FST and D", {
  dif <- differentiation_stats(ex("ibd_example.snps.vcf.gz"), ex("ibd_example_popmap.tsv"),
                               nboot = 0, beta = FALSE, verbose = FALSE)
  km <- ex("ibd_example_distances.csv")
  fst <- isolation_by_distance(dif, km, habitat = "1D", verbose = FALSE)
  ## Pinned: the example data and the exact test do not change between runs.
  expect_equal(fst$result$mantel_r, 0.8717925, tolerance = 1e-6)
  expect_equal(fst$result$p_value, 1 / 720)
  expect_equal(fst$result$slope, 0.000635, tolerance = 1e-3)
  expect_identical(fst$settings$negative_pairs, "site4 - site5")
  expect_equal(fst$pairs$FST, dif$pairwise_fst[cbind(fst$pairs$pop1, fst$pairs$pop2)])
  expect_match(fst$settings$source, "SNP VCF")

  d <- isolation_by_distance(dif, km, habitat = "1D", stat = "D", verbose = FALSE)
  expect_equal(d$pairs$D, dif$pairwise_D[cbind(d$pairs$pop1, d$pairs$pop2)])
  expect_gt(d$result$mantel_r, 0.8)
  expect_equal(d$result$p_value, 1 / 720)
})

test_that("a genetic matrix with one half blank, or names in a first column, gives the same result", {
  t <- toy(6, noise = 0.03)
  full <- isolation_by_distance(t$fst, t$km, habitat = "1D", verbose = FALSE)
  lower <- t$fst
  lower[upper.tri(lower)] <- NA                  # as many programs write pairwise FST
  diag(lower) <- NA
  a <- isolation_by_distance(lower, t$km, habitat = "1D", verbose = FALSE)
  expect_equal(a$result, full$result)
  expect_equal(a$pairs, full$pairs)
  upper <- t$fst
  upper[lower.tri(upper)] <- NA
  expect_equal(isolation_by_distance(upper, t$km, habitat = "1D", verbose = FALSE)$result,
               full$result)
  ## as read.csv() reads a square file without row.names = 1
  df <- data.frame(pop = rownames(t$fst), t$fst, check.names = FALSE)
  expect_equal(isolation_by_distance(df, t$km, habitat = "1D", verbose = FALSE)$result,
               full$result)
  ## a pair blank in both halves is still reported
  lower["p3", "p1"] <- NA
  expect_error(isolation_by_distance(lower, t$km, habitat = "1D", verbose = FALSE),
               "No FST for: p1 - p3")
})

test_that("a half-blank FST file with a blank diagonal, read with read.csv(), gives the same result", {
  ## read.csv() reads a column with no value at all (the last column of a
  ## lower half, the first of an upper half) as TRUE/FALSE, not as numbers.
  t <- toy(6, noise = 0.03)
  full <- isolation_by_distance(t$fst, t$km, habitat = "1D", verbose = FALSE)
  as_csv <- function(m) {
    cells <- ifelse(is.na(m), "", format(m, digits = 17))
    write_lines(c(paste0(",", paste(colnames(m), collapse = ",")),
                  paste0(rownames(m), ",", apply(cells, 1, paste, collapse = ","))))
  }
  lower <- t$fst
  lower[upper.tri(lower, diag = TRUE)] <- NA
  upper <- t$fst
  upper[lower.tri(upper, diag = TRUE)] <- NA
  for (half in list(lower, upper)) {
    f <- as_csv(half)
    df <- utils::read.csv(f)
    expect_true(any(vapply(df, is.logical, logical(1))))          # the case being tested
    expect_equal(isolation_by_distance(df, t$km, habitat = "1D", verbose = FALSE)$result,
                 full$result, tolerance = 1e-12)
    named <- utils::read.csv(f, row.names = 1)
    expect_equal(isolation_by_distance(named, t$km, habitat = "1D", verbose = FALSE)$result,
                 full$result, tolerance = 1e-12)
  }
})

test_that("8 or more populations use random orderings, however large nperm is", {
  t <- toy(8, noise = 0.06, seed = 4)
  res <- isolation_by_distance(t$fst, t$km, habitat = "1D", nperm = 50000, seed = 1,
                               verbose = FALSE)
  expect_false(res$result$exact)
  expect_equal(res$result$orderings, 50000)
})

test_that("the smallest possible p allows for symmetric distances", {
  ## 4 evenly spaced sites: the reverse order gives the same r, so even perfect
  ## isolation by distance cannot reach p = 1/24.
  p4 <- c("a", "b", "c", "d")
  km <- abs(outer(c(0, 10, 20, 30), c(0, 10, 20, 30), "-"))
  dimnames(km) <- list(p4, p4)
  res <- isolation_by_distance(km / 1000, km, habitat = "1D", verbose = FALSE)
  expect_equal(res$result$p_value, 2 / 24)
  expect_equal(res$settings$min_p, 2 / 24)
  out <- capture.output(print(summary(res)))
  expect_true(any(grepl("^  look +With these 4 populations and distances no ordering can give p",
                        out)))
  ## uneven spacing: the floor is 1/24 again
  km2 <- abs(outer(c(0, 10, 25, 31), c(0, 10, 25, 31), "-"))
  dimnames(km2) <- list(p4, p4)
  expect_equal(isolation_by_distance(km2 / 1000, km2, habitat = "1D", verbose = FALSE)$settings$min_p,
               1 / 24)
})

## ---------------------------------------------------------------------------
## print(), summary() and plot()
## ---------------------------------------------------------------------------

test_that("the short summary fits 80 columns, says the slope has no SE, and $tables$results is what it prints", {
  t <- toy(6, noise = 0.03)
  res <- isolation_by_distance(t$fst, t$km, habitat = "2D", verbose = FALSE)
  for (details in c(FALSE, TRUE)) {
    out <- capture.output(print(summary(res, details = details)))
    expect_lte(max(nchar(out)), 80)
    expect_true(any(grepl("slope (no SE)", out, fixed = TRUE)))
    expect_true(any(grepl("^CHECKS", out)))
    expect_true(any(grepl("Cite Mantel (1967) and Rousset (1997)", out, fixed = TRUE)))
    expect_true(any(grepl("4 x pi x density x sigma^2", out, fixed = TRUE)))
  }
  out <- capture.output(print(summary(res)))
  expect_lt(length(out), 45)
  tab <- summary(res)$tables$results
  expect_named(tab, c("genetic", "geographic", "populations", "pairs", "Mantel r", "p",
                      "orderings", "slope (no SE)"))
  for (cell in unlist(tab[1, ], use.names = FALSE))
    expect_true(any(grepl(cell, out, fixed = TRUE)))
  details <- capture.output(print(summary(res, details = TRUE)))
  expect_true(any(grepl("EVERY PAIR", details)))
  expect_true(any(grepl("Guillot & Rousset 2013", details)))

  out_d <- capture.output(print(summary(
    isolation_by_distance(t$fst, t$km, habitat = "1D", stat = "D", verbose = FALSE))))
  expect_true(any(grepl("only the FST slope estimates", out_d)))
  expect_true(any(grepl("Cite Mantel (1967) and Jost (2008)", out_d, fixed = TRUE)))
})

test_that("with 3 populations the summary says the test cannot detect isolation by distance", {
  t <- toy(3)
  out <- capture.output(print(summary(isolation_by_distance(t$fst, t$km, habitat = "1D",
                                                            verbose = FALSE))))
  expect_true(any(grepl("^  look +With these 3 populations and distances", out)))
  expect_match(gsub("\\s+", " ", paste(out, collapse = " ")), "no ordering can give p below 0.17")
})

test_that("print() gives the result and plot() draws it", {
  t <- toy(5)
  res <- isolation_by_distance(t$fst, t$km, habitat = "1D", verbose = FALSE)
  out <- capture.output(print(res))
  expect_match(out[1], "Isolation by distance: 5 populations, 10 pairs, 1-D habitat")
  expect_true(any(grepl("Mantel r = ", out)))
  pdf(NULL)
  on.exit(dev.off())
  expect_identical(plot(res, xlab = "River km", main = ""), res)
})
