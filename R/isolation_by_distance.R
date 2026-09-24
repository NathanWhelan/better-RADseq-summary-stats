###############################################################################
#
#  R/isolation_by_distance.R -- do populations farther apart differ more?
#
#  When individuals mostly mate near where they were born, populations far
#  apart exchange fewer genes than neighbours, so genetic differentiation
#  grows with geographic distance: isolation by distance (Wright 1943). This
#  file tests for that pattern with the pairwise FST or Jost's D from
#  differentiation_stats() and a table of distances between populations.
#
#    read_distances()         reads the distances -- a square matrix (one half
#                             may be blank) or a list of pairs -- and returns a
#                             symmetric matrix named by population
#    isolation_by_distance()  the Mantel test, and the regression slope
#
#  (Not to be confused with "IBD" = identity by descent, as in PLINK's PI_HAT;
#  see R/kinship.R. The words are written out in full here for that reason.)
#
#  WHAT IS COMPARED (Rousset 1997). Pairwise FST is used as FST/(1 - FST) and
#  compared with the distance in a 1-D habitat (a river, a coastline) or with
#  the natural log of the distance in a 2-D habitat. Under isolation by
#  distance these relationships are close to straight lines, with slope
#  1/(4 x density x sigma^2) in 1-D and 1/(4 x pi x density x sigma^2) in 2-D,
#  where sigma^2 is the mean squared parent-offspring distance along one axis.
#
#  WHY JOST'S D IS USED AS IT IS. Write Q0 for the chance that two genes from
#  one population are identical and Qr for two genes from populations at
#  distance r. Then FST/(1 - FST) = (Q0 - Qr)/(1 - Q0) and, for a pair of
#  populations, D = (Q0 - Qr)/Q0. Rousset's straight line is a property of
#  Q0 - Qr. Both statistics divide it by a number that is the same for every
#  pair when the populations are equally diverse, so D = FST/(1 - FST) x
#  Hs/(1 - Hs): D already has the straight-line shape. D/(1 - D) =
#  (Q0 - Qr)/Qr would divide by a number that falls with distance and bend
#  the line upward, so it is not offered. D's slope is the FST slope scaled by
#  Hs/(1 - Hs); it is reported, but only the FST slope estimates dispersal.
#  (Checked with exact identity probabilities in stepping-stone models; see
#  vignette("rationale"), "Isolation by distance".)
#
#  THE TEST (Mantel 1967). r is the correlation between the genetic and the
#  geographic values over all pairs of populations. The pairs share
#  populations, so r cannot be tested as an ordinary correlation. Instead the
#  population labels of the genetic table are shuffled and r recomputed, which
#  gives the values r takes when genetics has nothing to do with geography.
#  The p-value is one-sided, because isolation by distance predicts r > 0.
#  With 7 or fewer populations every ordering is used (at most 5,040), so the
#  test is exact; with more, `nperm` random orderings are used and
#  p = (b + 1)/(nperm + 1), where b counts those giving an r at least as large
#  as the data's (Phipson & Smyth 2010). Guillot & Rousset (2013) show that
#  this simple Mantel test is valid for isolation by distance, and that the
#  partial Mantel test is not; no partial test is offered.
#
#  SOURCES. Everything here is written from the publications above; no code
#  is taken from another package. The package tests compare the results with
#  vegan::mantel() (Oksanen et al.), which gives the same r and, with 7 or
#  fewer populations, the same p-value.
#
#  REFERENCES
#    Wright, S. (1943) Isolation by distance. Genetics 28: 114-138.
#    Mantel, N. (1967) The detection of disease clustering and a generalized
#      regression approach. Cancer Research 27: 209-220.
#    Rousset, F. (1997) Genetic differentiation and estimation of gene flow
#      from F-statistics under isolation by distance. Genetics 145: 1219-1228.
#    Phipson, B. & Smyth, G.K. (2010) Permutation P-values should never be
#      zero. Statistical Applications in Genetics and Molecular Biology 9: 39.
#    Guillot, G. & Rousset, F. (2013) Dismantling the Mantel tests. Methods in
#      Ecology and Evolution 4: 336-344.
#    Jost, L. (2008) GST and its relatives do not measure differentiation.
#      Molecular Ecology 17: 4015-4026.
#
###############################################################################

## Not exported. With this many populations or fewer, every ordering is tried
## (7! = 5,040), so the Mantel test is exact.
.ibd_max_exact <- 7L

## ---------------------------------------------------------------------------
## Reading the distances
## ---------------------------------------------------------------------------

#' Read distances between populations
#'
#' Reads the geographic distances for [isolation_by_distance()] and returns
#' them as a symmetric matrix named by population. Print the result to check
#' that your file was read the way you meant.
#'
#' Two layouts are accepted, as a CSV file (comma-separated) or a
#' tab-separated file. Population names must be spelled exactly as in the
#' popmap. Their order does not matter: rows and columns are matched by name.
#'
#' **A square matrix**, with population names across the first row and down
#' the first column. The top-left cell may be blank or missing. One half of the
#' matrix may be left blank; it is filled from the other half. If both halves
#' are filled, they must agree.
#'
#' ```
#' ,site1,site2,site3
#' site1,0,15,27
#' site2,,0,12
#' site3,,,0
#' ```
#'
#' **A list of pairs**, one row per pair of populations: two names and the
#' distance. The header row is optional, and the order of the two names does
#' not matter.
#'
#' ```
#' pop1,pop2,distance
#' site1,site2,15
#' site1,site3,27
#' site2,site3,12
#' ```
#'
#' Use any unit (km, m, river km); the slope from [isolation_by_distance()] is
#' per unit of your distances. Use the distance the organisms actually travel:
#' along the river for stream animals, for example, not a straight line.
#'
#' @param distances A file path (CSV or tab-separated), a data frame, a
#'   matrix with population names as row and column names, or a `dist`
#'   object with labels.
#' @param verbose Print a line saying what was read. Default `TRUE`.
#' @return A symmetric numeric matrix with population names as row and column
#'   names and 0 on the diagonal. A pair with no distance in the input is
#'   `NA` (and a message says so); [isolation_by_distance()] stops if it needs
#'   that pair.
#' @seealso [isolation_by_distance()]
#' @examples
#' km <- system.file("extdata", "ibd_example_distances.csv", package = "RADdiversity")
#' read_distances(km)
#'
#' # The same distances as a list of pairs
#' pairs <- data.frame(pop1 = c("a", "a", "b"), pop2 = c("b", "c", "c"),
#'                     distance = c(15, 27, 12))
#' read_distances(pairs)
#' @export
read_distances <- function(distances, verbose = TRUE) {
  .check_flag(verbose, "verbose")
  if (inherits(distances, "dist")) {
    if (is.null(attr(distances, "Labels")))
      stop("The dist object has no population names (labels). Give it labels, ",
           "or give the distances as a file or a matrix with names.", call. = FALSE)
    distances <- as.matrix(distances)
  }
  if (is.character(distances)) {
    .check_string(distances, "distances")
    if (!file.exists(distances))
      stop("Distance file not found: ", distances, "\n  Check the path and try again.", call. = FALSE)
    cells <- .distance_file_cells(distances)
    source <- basename(distances)
  } else if (is.matrix(distances) || is.data.frame(distances)) {
    cells <- .distance_table_cells(distances)
    source <- if (is.matrix(distances)) "the matrix" else "the data frame"
  } else {
    stop("`distances` must be a file path, a data frame, a matrix or a dist object (got an ",
         "object of class ", paste(class(distances), collapse = "/"), ").", call. = FALSE)
  }
  parsed <- .parse_distance_cells(cells, source)
  m <- parsed$matrix
  n <- nrow(m)
  missing_pairs <- .pair_labels(m, is.na(m) & upper.tri(m))
  .inform(verbose, sprintf("Read distances between %d populations (%d pairs) from %s in %s%s.",
                           n, n * (n - 1L) / 2L, parsed$layout, source, parsed$note))
  if (length(missing_pairs))
    .inform(verbose, sprintf("  %d pair(s) have no distance: %s", length(missing_pairs),
                             .shorten(missing_pairs)))
  m
}

## Not exported. The cells of a distance file as a character matrix (blank
## cells are NA), with the file's first row as the first row. The separator
## is a TAB if the first line has one, else a comma, else a semicolon. A byte
## order mark (written by Excel) and Windows line endings are removed.
.distance_file_cells <- function(path) {
  con <- file(path, encoding = "UTF-8-BOM")
  lines <- tryCatch(readLines(con, warn = FALSE), finally = close(con))
  lines <- sub("\r$", "", lines)
  lines <- lines[nzchar(trimws(lines))]
  if (!length(lines)) stop("The distance file is empty: ", path, call. = FALSE)
  sep <- if (grepl("\t", lines[1], fixed = TRUE)) "\t"
         else if (grepl(",", lines[1], fixed = TRUE)) ","
         else if (grepl(";", lines[1], fixed = TRUE)) ";"
         else stop("Could not find the columns of ", path, ": separate them with commas ",
                   "(a CSV file) or TABs. See ?read_distances.", call. = FALSE)
  ## Read with as many columns as the longest line, so that a header row
  ## without a top-left cell cannot shift the others.
  n_cols <- max(utils::count.fields(textConnection(lines), sep = sep, quote = "\"",
                                    comment.char = "", blank.lines.skip = TRUE), na.rm = TRUE)
  cells <- utils::read.table(text = lines, sep = sep, header = FALSE, quote = "\"",
                             colClasses = "character", na.strings = c("", "NA"),
                             strip.white = TRUE, fill = TRUE, comment.char = "",
                             col.names = paste0("V", seq_len(n_cols)))
  as.matrix(cells)
}

## Not exported. The cells of a data frame or matrix of distances, in the same
## form as .distance_file_cells(): the column names become the first row, and
## meaningful row names become the first column. Numbers are written with 17
## significant digits, so nothing is lost on the way back.
.distance_table_cells <- function(x) {
  as_text <- function(v) if (is.numeric(v)) ifelse(is.na(v), NA_character_, sprintf("%.17g", v))
                         else trimws(as.character(v))
  if (is.matrix(x)) {
    if (is.null(rownames(x)) || is.null(colnames(x)))
      stop("A distance matrix needs population names as its row and column names.", call. = FALSE)
    x <- as.data.frame(x, stringsAsFactors = FALSE, check.names = FALSE)
  }
  body <- vapply(x, as_text, character(nrow(x)))
  if (!is.matrix(body)) body <- matrix(body, nrow = nrow(x))
  has_row_names <- !identical(rownames(x), as.character(seq_len(nrow(x))))
  if (has_row_names) rbind(c(NA_character_, names(x)), cbind(rownames(x), body))
  else rbind(names(x), body)
}

## Not exported. Turns the cells into the symmetric distance matrix. Tries the
## square layout first, then the list of pairs; stops with both layouts shown
## when neither fits. Returns list(matrix, layout, note).
.parse_distance_cells <- function(cells, source) {
  cells[] <- trimws(cells)
  cells[!is.na(cells) & !nzchar(cells)] <- NA_character_
  ## Drop empty rows and columns (Excel often adds them).
  cells <- cells[rowSums(!is.na(cells)) > 0L, colSums(!is.na(cells)) > 0L, drop = FALSE]
  square <- if (nrow(cells) >= 3L) .distance_square_layout(cells) else NULL
  if (!is.null(square)) return(.distances_from_square(square$values, square$names, source))
  if (ncol(cells) == 3L && nrow(cells) >= 1L) return(.distances_from_pairs(cells, source))
  stop("Could not read ", source, " as distances. Use one of these two layouts ",
       "(names as in the popmap):\n",
       "  a square matrix           or   a list of pairs\n",
       "  ,site1,site2,site3             pop1,pop2,distance\n",
       "  site1,0,15,27                  site1,site2,15\n",
       "  site2,,0,12                    site1,site3,27\n",
       "  site3,,,0                      site2,site3,12\n",
       "See ?read_distances.", call. = FALSE)
}

## Not exported. If the cells are a square matrix -- the names down the first
## column are the same set as the names across the first row, with or without
## a top-left cell -- returns list(names, values) with the values' columns in
## the order of `names`. Otherwise NULL.
.distance_square_layout <- function(cells) {
  row_names <- cells[-1L, 1L]
  m <- length(row_names)
  if (anyNA(row_names) || anyDuplicated(row_names) || ncol(cells) < m + 1L) return(NULL)
  header <- cells[1L, ]
  ## With a top-left cell (blank or a label), the names are in columns 2..m+1;
  ## without one, the first row is one cell shorter and they are in 1..m.
  col_names <- if (!is.na(header[1L]) && header[1L] %in% row_names && ncol(cells) == m + 1L &&
                   is.na(header[m + 1L])) header[seq_len(m)]
               else header[1L + seq_len(m)]
  if (anyNA(col_names) || anyDuplicated(col_names) || !setequal(col_names, row_names)) return(NULL)
  if (ncol(cells) > m + 1L) return(NULL)
  values <- cells[-1L, 1L + seq_len(m), drop = FALSE]
  list(names = row_names, values = values[, match(row_names, col_names), drop = FALSE])
}

## Not exported. The symmetric matrix from the values of a square layout (a
## character or numeric matrix whose rows and columns follow `pops`).
.distances_from_square <- function(values, pops, source) {
  num <- suppressWarnings(matrix(as.numeric(values), nrow(values)))
  bad <- which(!is.na(values) & is.na(num), arr.ind = TRUE)
  if (nrow(bad))
    stop(sprintf("In %s, the distance from %s to %s is \"%s\", which is not a number.",
                 source, pops[bad[1, 1]], pops[bad[1, 2]], values[bad[1, 1], bad[1, 2]]),
         "\n  Use a point for decimals (12.5, not 12,5) and no units or thousands separators.",
         call. = FALSE)
  self <- diag(num)
  if (any(!is.na(self) & self != 0)) {
    k <- which(!is.na(self) & self != 0)[1L]
    stop(sprintf("In %s, the distance from %s to itself is %s; it must be 0 or blank.",
                 source, pops[k], format(self[k])),
         "\n  A non-zero value on the diagonal usually means the rows and columns do not line up.",
         call. = FALSE)
  }
  mirror <- t(num)
  both <- !is.na(num) & !is.na(mirror) & upper.tri(num)
  differ <- both & abs(num - mirror) > 1e-8 * pmax(1, abs(num))
  if (any(differ)) {
    k <- which(differ, arr.ind = TRUE)[1L, ]
    stop(sprintf("In %s, %s - %s is %s in one half of the matrix and %s in the other.",
                 source, pops[k[1]], pops[k[2]], format(num[k[1], k[2]]), format(mirror[k[1], k[2]])),
         "\n  Fill in one half only, or make the two halves agree.", call. = FALSE)
  }
  filled <- ifelse(is.na(num), mirror, num)
  diag(filled) <- 0
  dimnames(filled) <- list(pops, pops)
  .check_distance_values(filled, source)
  note <- if (all(is.na(num[lower.tri(num)])) && any(!is.na(num[upper.tri(num)])))
            "; the lower half was blank and was filled from the upper half"
          else if (all(is.na(num[upper.tri(num)])) && any(!is.na(num[lower.tri(num)])))
            "; the upper half was blank and was filled from the lower half"
          else ""
  list(matrix = filled, layout = "a square matrix", note = note)
}

## Not exported. The symmetric matrix from a list of pairs: three columns,
## two population names and a distance. A first row whose third cell holds
## text that is not a number is taken as the header.
.distances_from_pairs <- function(cells, source) {
  first_row <- 1L
  if (!is.na(cells[1L, 3L]) && is.na(suppressWarnings(as.numeric(cells[1L, 3L])))) {
    cells <- cells[-1L, , drop = FALSE]
    first_row <- 2L
  }
  if (!nrow(cells)) stop(source, " has a header but no pairs.", call. = FALSE)
  row_label <- function(k) sprintf("row %d (%s, %s)", k + first_row - 1L, cells[k, 1L], cells[k, 2L])
  no_name <- which(is.na(cells[, 1L]) | is.na(cells[, 2L]))
  if (length(no_name))
    stop(sprintf("In %s, %s is missing a population name.", source, row_label(no_name[1L])),
         call. = FALSE)
  d <- suppressWarnings(as.numeric(cells[, 3L]))
  bad <- which(is.na(d))
  if (length(bad))
    stop(sprintf("In %s, %s has %s as its distance, which is not a number.", source,
                 row_label(bad[1L]), if (is.na(cells[bad[1L], 3L])) "a blank"
                                     else paste0("\"", cells[bad[1L], 3L], "\"")),
         "\n  Use a point for decimals (12.5, not 12,5) and no units or thousands separators.",
         call. = FALSE)
  pops <- unique(c(rbind(cells[, 1L], cells[, 2L])))
  m <- matrix(NA_real_, length(pops), length(pops), dimnames = list(pops, pops))
  diag(m) <- 0
  for (k in seq_len(nrow(cells))) {
    i <- cells[k, 1L]
    j <- cells[k, 2L]
    if (i == j) {
      if (d[k] != 0)
        stop(sprintf("In %s, %s gives a distance from a population to itself; it must be 0.",
                     source, row_label(k)), call. = FALSE)
      next
    }
    if (!is.na(m[i, j]) && abs(m[i, j] - d[k]) > 1e-8 * max(1, abs(d[k])))
      stop(sprintf("In %s, the pair %s - %s is listed twice, with %s and %s.", source, i, j,
                   format(m[i, j]), format(d[k])), call. = FALSE)
    m[i, j] <- m[j, i] <- d[k]
  }
  .check_distance_values(m, source)
  list(matrix = m, layout = "a list of pairs", note = "")
}

## Not exported. Stops on a negative distance.
.check_distance_values <- function(m, source) {
  neg <- which(!is.na(m) & m < 0 & upper.tri(m), arr.ind = TRUE)
  if (nrow(neg))
    stop(sprintf("In %s, the distance between %s and %s is negative (%s). Distances must be 0 or more.",
                 source, rownames(m)[neg[1, 1]], colnames(m)[neg[1, 2]],
                 format(m[neg[1, 1], neg[1, 2]])), call. = FALSE)
  invisible(NULL)
}

## Not exported. "a - b" labels for the cells of a square matrix flagged in
## `which_cells` (a logical matrix), row population first.
.pair_labels <- function(m, which_cells) {
  k <- which(which_cells, arr.ind = TRUE)
  if (!nrow(k)) return(character(0))
  k <- k[order(k[, 1], k[, 2]), , drop = FALSE]
  paste(rownames(m)[k[, 1]], "-", colnames(m)[k[, 2]])
}

## Not exported. The first few items of a list, then "and N more".
.shorten <- function(x, n = 5L) {
  if (length(x) <= n) return(paste(x, collapse = ", "))
  paste0(paste(x[seq_len(n)], collapse = ", "), sprintf(", and %d more", length(x) - n))
}

## ---------------------------------------------------------------------------
## The test
## ---------------------------------------------------------------------------

#' Isolation by distance: do populations farther apart differ more?
#'
#' Tests whether genetic differentiation between populations grows with the
#' geographic distance between them, with a Mantel test (Mantel 1967), and
#' gives the slope of that relationship (Rousset 1997). It uses the pairwise
#' FST or Jost's D from [differentiation_stats()] and a table of distances
#' read by [read_distances()].
#'
#' @details
#' **What is compared.** With `stat = "FST"` (the default), each pair's FST
#' is used as FST/(1 - FST), as Rousset (1997) showed. In a 1-D habitat it is
#' compared with the distance, in a 2-D habitat with the natural log of the
#' distance. Under isolation by distance both relationships are close to
#' straight lines, and the slope estimates dispersal: 1/slope estimates
#' 4 x density x sigma^2 in 1-D and 4 x pi x density x sigma^2 in 2-D, where
#' density is individuals per unit of length or area and sigma^2 is the mean
#' squared distance between parent and offspring along one axis. These hold
#' best for pairs farther apart than sigma, in populations near equilibrium.
#'
#' **Jost's D** (`stat = "D"`) is used as it is. For a pair of populations
#' with similar diversity, D equals FST/(1 - FST) times Hs/(1 - Hs), so D
#' already has the straight-line shape; D/(1 - D) would bend it and is not
#' offered. The slope of D is reported to describe the pattern, but only the
#' FST slope estimates dispersal (see `vignette("rationale")`).
#'
#' **The test.** `r` is the correlation between the genetic and the geographic
#' values over all pairs of populations. Pairs share populations, so `r`
#' cannot be tested as an ordinary correlation. Instead, the populations are
#' shuffled among the geographic positions and `r` recomputed each time; `p`
#' is the share of orderings giving an `r` at least as large as the data's
#' (one-sided: isolation by distance predicts a positive `r`). With 7 or fewer
#' populations every ordering is tried (at most 5,040), so `p` is exact and
#' the same on every run. With more, `nperm` random orderings are used and
#' p = (b + 1)/(nperm + 1), where b counts those at least as large (Phipson &
#' Smyth 2010). With few populations the test has little power: with 3
#' populations `p` can never be below 1/6.
#'
#' The simple Mantel test is valid for isolation by distance (Guillot &
#' Rousset 2013). Partial Mantel tests, which try to control for a third
#' distance matrix, are not valid for spatially structured data and are not
#' offered.
#'
#' **The slope** is an ordinary least-squares line through all pairs, a point
#' estimate with no standard error here. Genepop (Rousset 2008) gives a
#' bootstrap interval for it.
#'
#' **Negative values.** A pairwise FST or D estimate can be slightly below 0
#' when two populations barely differ. Such values are kept as they are:
#' setting them to 0 would bias the slope.
#'
#' @param genetic The result of [differentiation_stats()], or a square matrix
#'   of pairwise values with population names as row and column names (for
#'   example FST from another program). One half of the matrix may be blank
#'   (`NA`), as many programs write it; it is filled from the other half. A
#'   data frame with the names in its first column also works.
#' @param distances The distances between populations: a file path, data
#'   frame, matrix or `dist` object, in either layout that [read_distances()]
#'   accepts. Population names must match those in `genetic`; their order does
#'   not matter.
#' @param habitat Required. `"1D"` for populations along a line -- a river,
#'   stream, coastline or other narrow strip -- where genetic values are
#'   compared with the distance. `"2D"` for populations spread over an area,
#'   where they are compared with the natural log of the distance (every
#'   distance must then be above 0).
#' @param stat `"FST"` (the default), used as FST/(1 - FST), or `"D"` (Jost's
#'   D), used as it is.
#' @param method `"pearson"` (the default) for the correlation of the values,
#'   or `"spearman"` for the correlation of their ranks, which is less swayed
#'   by a few extreme pairs. The slope always uses the values.
#' @param nperm Random orderings of the populations when there are 8 or more.
#'   Default `9999`. With 7 or fewer every ordering is used and `nperm` is
#'   ignored.
#' @param exclude Names of populations to leave out, for example one far from
#'   the others, to see how much it drives the result.
#' @param seed `NULL` (the default: use the session's random numbers; call
#'   `set.seed()` first for reproducible results), or a whole number, used
#'   for this call only.
#' @param verbose Print progress messages. Default `TRUE`.
#' @return An object of class `raddiv_ibd`, a list of:
#'   \describe{
#'     \item{result}{One row: what was compared (`genetic`, `geographic`,
#'       `habitat`), `populations`, `pairs`, `method`, `mantel_r`, `p_value`,
#'       `orderings` (how many were used), `exact` (`TRUE` when every ordering
#'       was used), `slope` and `intercept`.}
#'     \item{pairs}{One row per pair of populations: `pop1`, `pop2`,
#'       `distance`, `log_distance` (2-D only), and `FST` and `FST_linearized`
#'       (= FST/(1 - FST)), or `D`.}
#'     \item{settings}{The settings of this run, the populations used and
#'       those left out.}
#'   }
#'   `summary()` shows the result, what to report and a list of checks;
#'   `summary(x)$tables$results` is that table as a data frame, for a paper.
#'   `plot()` draws the pairs and the fitted line.
#' @references
#' Wright, S. (1943) Isolation by distance. *Genetics* 28:114-138.
#' \doi{10.1093/genetics/28.2.114}
#'
#' Mantel, N. (1967) The detection of disease clustering and a generalized
#' regression approach. *Cancer Research* 27:209-220.
#'
#' Rousset, F. (1997) Genetic differentiation and estimation of gene flow from
#' F-statistics under isolation by distance. *Genetics* 145:1219-1228.
#' \doi{10.1093/genetics/145.4.1219}
#'
#' Rousset, F. (2008) genepop'007: a complete re-implementation of the genepop
#' software for Windows and Linux. *Molecular Ecology Resources* 8:103-106.
#' \doi{10.1111/j.1471-8286.2007.01931.x}
#'
#' Phipson, B. & Smyth, G.K. (2010) Permutation P-values should never be zero:
#' calculating exact P-values when permutations are randomly drawn.
#' *Statistical Applications in Genetics and Molecular Biology* 9:39.
#' \doi{10.2202/1544-6115.1585}
#'
#' Guillot, G. & Rousset, F. (2013) Dismantling the Mantel tests. *Methods in
#' Ecology and Evolution* 4:336-344. \doi{10.1111/2041-210x.12018}
#'
#' Jost, L. (2008) GST and its relatives do not measure differentiation.
#' *Molecular Ecology* 17:4015-4026. \doi{10.1111/j.1365-294X.2008.03887.x}
#' @seealso [read_distances()] for the distance file, [differentiation_stats()]
#'   for FST and D.
#' @examples
#' # A small made-up example: 4 populations along a river.
#' fst <- matrix(c(0,    0.02, 0.05, 0.08,
#'                 0.02, 0,    0.03, 0.06,
#'                 0.05, 0.03, 0,    0.02,
#'                 0.08, 0.06, 0.02, 0), 4,
#'               dimnames = list(c("a", "b", "c", "d"), c("a", "b", "c", "d")))
#' km <- data.frame(pop1 = c("a", "a", "a", "b", "b", "c"),
#'                  pop2 = c("b", "c", "d", "c", "d", "d"),
#'                  distance = c(10, 25, 40, 15, 30, 15))
#' isolation_by_distance(fst, km, habitat = "1D", verbose = FALSE)
#'
#' # From VCF to test, with the example data: six sites along a river.
#' vcf    <- system.file("extdata", "ibd_example.snps.vcf.gz", package = "RADdiversity")
#' popmap <- system.file("extdata", "ibd_example_popmap.tsv", package = "RADdiversity")
#' km     <- system.file("extdata", "ibd_example_distances.csv", package = "RADdiversity")
#' dif <- differentiation_stats(vcf, popmap, nboot = 0, beta = FALSE, verbose = FALSE)
#' ibd <- isolation_by_distance(dif, km, habitat = "1D")
#' summary(ibd)
#' plot(ibd, xlab = "River km")
#' isolation_by_distance(dif, km, habitat = "1D", stat = "D", verbose = FALSE)
#' @export
isolation_by_distance <- function(genetic, distances, habitat, stat = "FST", method = "pearson",
                                  nperm = 9999L, exclude = NULL, seed = NULL, verbose = TRUE) {
  ## ---- 1. Check the arguments ----------------------------------------------
  if (missing(habitat))
    stop("Say what kind of habitat the populations live in:\n",
         "  habitat = \"1D\" for a river, stream, coastline or other narrow strip\n",
         "           (genetic values are compared with the distance);\n",
         "  habitat = \"2D\" for an area, such as a forest, grassland or lake bottom\n",
         "           (they are compared with the log of the distance; Rousset 1997).",
         call. = FALSE)
  if (missing(genetic) || missing(distances))
    stop("Give `genetic` (the result of differentiation_stats()) and `distances` ",
         "(see ?read_distances).", call. = FALSE)
  .check_choice(habitat, "habitat", c("1D", "2D"))
  .check_choice(stat, "stat", c("FST", "D"))
  .check_choice(method, "method", c("pearson", "spearman"))
  nperm <- .check_count(nperm, "nperm", min = 1)
  .check_flag(verbose, "verbose")
  .check_seed(seed)
  if (!is.null(exclude) && !(is.character(exclude) && !anyNA(exclude)))
    stop("`exclude` must be population names, e.g. exclude = c(\"site3\") (got: ",
         .show_value(exclude), ").", call. = FALSE)
  if (!is.null(seed)) {
    restore_rng <- .save_rng_state()
    on.exit(restore_rng(), add = TRUE)
    set.seed(seed)
  }

  ## ---- 2. The genetic and the geographic tables ------------------------------
  gen <- .ibd_genetic_matrix(genetic, stat)
  geo <- read_distances(distances, verbose = verbose && is.character(distances))

  ## ---- 3. Match the populations by name ---------------------------------------
  all_pops <- rownames(gen$matrix)
  unknown <- setdiff(exclude, all_pops)
  if (length(unknown))
    stop("`exclude` names population(s) not in the genetic data: ", paste(unknown, collapse = ", "),
         .near_match_hint(unknown, all_pops),
         "\n  The populations are: ", paste(all_pops, collapse = ", "), call. = FALSE)
  pops <- setdiff(all_pops, exclude)
  if (length(pops) < 3L)
    stop("Isolation by distance needs at least 3 populations (3 pairs); ",
         if (length(exclude)) "after `exclude`, " else "", "there ",
         if (length(pops) == 1L) "is 1." else paste0("are ", length(pops), "."), call. = FALSE)
  no_distance <- setdiff(pops, rownames(geo))
  if (length(no_distance))
    stop("These populations are in the genetic data but not in the distances: ",
         paste(no_distance, collapse = ", "),
         .near_match_hint(no_distance, setdiff(rownames(geo), pops)),
         "\n  Names must be spelled exactly as in the popmap. Add the population(s) to the ",
         "distances, or leave them out with exclude = c(",
         paste0("\"", no_distance, "\"", collapse = ", "), ").", call. = FALSE)
  dropped <- setdiff(rownames(geo), all_pops)
  if (length(dropped))
    .inform(verbose, "  In the distances but not in the genetic data (left out): ",
            paste(dropped, collapse = ", "))

  ## ---- 4. One value per pair ---------------------------------------------------
  idx <- utils::combn(length(pops), 2L)
  pair_pop1 <- pops[idx[1L, ]]
  pair_pop2 <- pops[idx[2L, ]]
  g <- gen$matrix[cbind(pair_pop1, pair_pop2)]
  d <- geo[cbind(pair_pop1, pair_pop2)]
  no_value <- is.na(g)
  if (any(no_value))
    stop(sprintf("No %s for: %s. ", stat, .shorten(paste(pair_pop1, "-", pair_pop2)[no_value])),
         if (gen$from_result) "differentiation_stats() found no record typed in both populations of these pairs. "
         else "", "Leave one population of each pair out with `exclude`.", call. = FALSE)
  if (any(is.na(d)))
    stop("No distance for: ", .shorten(paste(pair_pop1, "-", pair_pop2)[is.na(d)]),
         ". Add these pairs to the distances.", call. = FALSE)
  if (stat == "FST" && any(g >= 1))
    stop("FST is 1 or more for: ", .shorten(paste(pair_pop1, "-", pair_pop2)[g >= 1]),
         ", so FST/(1 - FST) is not defined.", call. = FALSE)
  if (habitat == "2D" && any(d <= 0))
    stop("With habitat = \"2D\" the log of the distance is used, and these pairs are at ",
         "distance 0: ", .shorten(paste(pair_pop1, "-", pair_pop2)[d <= 0]),
         ". Merge populations sampled at the same place, or leave one out with `exclude`.",
         call. = FALSE)
  y <- if (stat == "FST") g / (1 - g) else g
  x <- if (habitat == "2D") log(d) else d

  ## ---- 5. Mantel test and slope ----------------------------------------------
  n_pops <- length(pops)
  exact <- n_pops <= .ibd_max_exact
  .inform(verbose, sprintf("Isolation by distance: %d populations, %d pairs; Mantel test with %s ...",
                           n_pops, length(y),
                           if (exact) sprintf("all %s orderings (exact)", .big(factorial(n_pops)))
                           else sprintf("%s random orderings", .big(nperm))))
  mantel <- .mantel_test(y, x, idx, n_pops, method, exact, nperm)
  x_centred <- x - mean(x)
  slope <- sum(x_centred * (y - mean(y))) / sum(x_centred^2)
  intercept <- mean(y) - slope * mean(x)

  ## ---- 6. The result ----------------------------------------------------------
  pairs <- data.frame(pop1 = pair_pop1, pop2 = pair_pop2, distance = d, stringsAsFactors = FALSE)
  if (habitat == "2D") pairs$log_distance <- x
  if (stat == "FST") {
    pairs$FST <- g
    pairs$FST_linearized <- y
  } else {
    pairs$D <- g
  }
  result <- data.frame(genetic = if (stat == "FST") "FST/(1-FST)" else "D",
                       geographic = if (habitat == "2D") "ln(distance)" else "distance",
                       habitat = habitat, populations = n_pops, pairs = length(y), method = method,
                       mantel_r = mantel$r, p_value = mantel$p, orderings = mantel$orderings,
                       exact = exact, slope = slope, intercept = intercept,
                       stringsAsFactors = FALSE)
  settings <- list(habitat = habitat, stat = stat, method = method, nperm = nperm, exact = exact,
                   orderings = mantel$orderings, min_p = mantel$min_p, seed = seed,
                   populations = pops, excluded = exclude, dropped_from_distances = dropped,
                   negative_pairs = paste(pair_pop1, "-", pair_pop2)[g < 0],
                   source = gen$source)
  structure(list(result = result, pairs = pairs, settings = settings), class = "raddiv_ibd")
}

## Not exported. The genetic values as a square matrix named by population,
## from a differentiation_stats() result or from a matrix the user gives.
## Returns list(matrix, from_result, source).
.ibd_genetic_matrix <- function(genetic, stat) {
  if (inherits(genetic, "raddiv_differentiation")) {
    m <- if (stat == "FST") genetic$pairwise_fst else genetic$pairwise_D
    return(list(matrix = m, from_result = TRUE,
                source = sprintf("differentiation_stats() on a %s VCF",
                                 if (isTRUE(genetic$settings$is_haplotype)) "haplotype" else "SNP")))
  }
  if (is.data.frame(genetic)) {
    ## A column with no value at all (the last column of a lower half with a
    ## blank diagonal, or the first of an upper half) is read by read.csv() as
    ## TRUE/FALSE rather than as numbers. It holds only blanks, so it is
    ## numeric NA.
    blank <- vapply(genetic, function(v) is.logical(v) && all(is.na(v)), logical(1))
    genetic[blank] <- lapply(genetic[blank], as.numeric)
    ## As read.csv() reads a square file without row.names = 1: the names are
    ## in the first column.
    first <- genetic[[1L]]
    if (ncol(genetic) > 1L && (is.character(first) || is.factor(first)) &&
        all(vapply(genetic[-1L], is.numeric, logical(1)))) {
      genetic <- as.matrix(genetic[-1L])
      rownames(genetic) <- trimws(as.character(first))
    } else {
      genetic <- as.matrix(genetic)
    }
  }
  if (!(is.matrix(genetic) && is.numeric(genetic)))
    stop("`genetic` must be the result of differentiation_stats(), or a square matrix of ",
         "pairwise values with population names (got an object of class ",
         paste(class(genetic), collapse = "/"), ").", call. = FALSE)
  if (nrow(genetic) != ncol(genetic) || is.null(rownames(genetic)) || is.null(colnames(genetic)) ||
      !identical(sort(rownames(genetic)), sort(colnames(genetic))) || anyDuplicated(rownames(genetic)))
    stop("`genetic` must be a square matrix with the same population names as its row and ",
         "column names.", call. = FALSE)
  genetic <- genetic[, rownames(genetic), drop = FALSE]
  off <- row(genetic) != col(genetic)
  asym <- off & !is.na(genetic) & !is.na(t(genetic)) & abs(genetic - t(genetic)) > 1e-8
  if (any(asym))
    stop("`genetic` is not symmetric: ", .pair_labels(genetic, asym & upper.tri(genetic))[1L],
         " differs between its two halves.", call. = FALSE)
  ## Many programs fill in only one half of the matrix: take each blank cell
  ## from its mirror image.
  filled <- ifelse(is.na(genetic), t(genetic), genetic)
  dimnames(filled) <- dimnames(genetic)
  list(matrix = filled, from_result = FALSE, source = "a matrix given by the user")
}

## Not exported. For names that were not found, a hint naming any candidate
## that differs only in capitals, spaces, dots, dashes or underscores ("" if
## there is none).
.near_match_hint <- function(missing_names, candidates) {
  key <- function(s) gsub("[[:space:]._-]+", "", tolower(s))
  hits <- vapply(missing_names, function(nm) {
    same <- candidates[key(candidates) == key(nm)]
    if (length(same)) sprintf("%s (did you mean \"%s\"?)", nm, same[1L]) else NA_character_
  }, character(1))
  hits <- hits[!is.na(hits)]
  if (length(hits)) paste0("\n  Close matches: ", paste(hits, collapse = "; ")) else ""
}

## Not exported. The one-sided Mantel test. `y` and `x` are the genetic and
## geographic values of the pairs, in the order of `idx` (utils::combn() of
## the population numbers). An ordering gives geographic position k the
## genetic values of population ordering[k]; r is recomputed each time.
## Every ordering pairs the same set of genetic values with the distances in a
## different arrangement, so their mean and spread do not change and r is a
## sum of products divided by a fixed number. With `exact`, every ordering of
## the n_pops populations is used (the data's own ordering among them), and
## p = (orderings with r at least the data's) / n_pops!; `min_p` is then the
## smallest p these genetic values could get with these distances. Otherwise `nperm`
## random orderings, and p = (b + 1)/(nperm + 1) (Phipson & Smyth 2010).
## "At least the data's" allows for rounding error, so an ordering that gives
## the same r by a different sum is counted. Spearman: the same, on ranks.
.mantel_test <- function(y, x, idx, n_pops, method, exact, nperm) {
  if (method == "spearman") {
    y <- rank(y)
    x <- rank(x)
  }
  ## The genetic values as a symmetric matrix, so an ordering can be applied to
  ## its rows and columns together.
  Y <- matrix(0, n_pops, n_pops)
  Y[t(idx)] <- y
  Y[t(idx[2:1, , drop = FALSE])] <- y
  x_centred <- x - mean(x)
  scale <- sqrt(sum(x_centred^2) * sum((y - mean(y))^2))
  r_of <- function(orderings) {                      # one row per ordering
    cells <- (orderings[, idx[2L, ], drop = FALSE] - 1L) * n_pops + orderings[, idx[1L, ], drop = FALSE]
    moved <- Y[as.vector(cells)]
    drop(matrix(moved, nrow(orderings)) %*% x_centred) / scale
  }
  r <- r_of(matrix(seq_len(n_pops), 1L))
  if (!is.finite(r))
    stop("The Mantel r cannot be computed: every pair has the same ",
         if (stats::sd(x) == 0) "distance." else "genetic value.", call. = FALSE)
  tolerance <- sqrt(.Machine$double.eps)
  if (exact) {
    all_orderings <- .all_orderings(n_pops)
    r_all <- r_of(all_orderings)
    orderings <- nrow(all_orderings)
    p <- sum(r_all >= r - tolerance) / orderings
    ## The smallest p these data allow: the share of orderings that give the
    ## largest r. It is 1/n! unless several orderings give that r, as when the
    ## distances are symmetric (evenly spaced sites give the same r in reverse).
    min_p <- sum(r_all >= max(r_all) - tolerance) / orderings
  } else {
    at_least <- 0
    done <- 0L
    while (done < nperm) {                           # in blocks, to bound memory
      block <- min(1000L, nperm - done)
      random <- t(vapply(seq_len(block), function(i) sample.int(n_pops), integer(n_pops)))
      at_least <- at_least + sum(r_of(random) >= r - tolerance)
      done <- done + block
    }
    orderings <- nperm
    p <- (at_least + 1) / (nperm + 1)
    min_p <- 1 / (nperm + 1)
  }
  list(r = r, p = p, orderings = orderings, min_p = min_p)
}

## Not exported. Every ordering of 1..n, one per row (n! rows).
.all_orderings <- function(n) {
  if (n == 1L) return(matrix(1L, 1L, 1L))
  shorter <- .all_orderings(n - 1L)
  do.call(rbind, lapply(seq_len(n), function(first) {
    rest <- setdiff(seq_len(n), first)
    cbind(first, matrix(rest[as.vector(shorter)], nrow(shorter)), deparse.level = 0)
  }))
}

## ---------------------------------------------------------------------------
## print(), summary() and plot()
## ---------------------------------------------------------------------------

## Not exported. How the orderings are described, in the RESULTS table
## ("all 720 (exact)", "9,999 random") and in a sentence ("all 720 orderings
## (an exact test)", "9,999 random orderings").
.ibd_orderings_text <- function(st, sentence = FALSE) {
  if (sentence)
    return(if (st$exact) sprintf("all %s orderings (an exact test)", .big(st$orderings))
           else sprintf("%s random orderings", .big(st$orderings)))
  if (st$exact) sprintf("all %s (exact)", .big(st$orderings)) else sprintf("%s random", .big(st$orderings))
}

## Not exported. "1-D" or "2-D", for printing.
.ibd_habitat_text <- function(habitat) if (habitat == "2D") "2-D" else "1-D"

## Not exported. A slope for printing: 3 significant digits, never in
## scientific notation (0.000512, 12.3).
.ibd_format_slope <- function(x) ifelse(is.finite(x), formatC(signif(x, 3), digits = 3, format = "fg"), "NA")

#' @rdname isolation_by_distance
#' @param x,object A `raddiv_ibd` object, as returned by
#'   `isolation_by_distance()`.
#' @param ... For `plot()`: graphical settings passed to [graphics::plot()],
#'   such as `xlab = "River km"` or `main = ""`. Otherwise ignored.
#' @export
print.raddiv_ibd <- function(x, ...) {
  st <- x$settings
  res <- x$result
  cat(sprintf("Isolation by distance: %d populations, %d pairs, %s habitat\n",
              res$populations, res$pairs, .ibd_habitat_text(st$habitat)))
  cat(sprintf("  %s against %s; one-sided Mantel test (%s),\n  %s\n",
              res$genetic, res$geographic, if (st$method == "pearson") "Pearson" else "Spearman",
              .ibd_orderings_text(st, sentence = TRUE)))
  cat(sprintf("  Mantel r = %.3f, p = %s; slope = %s (no SE)\n",
              res$mantel_r, .format_p(res$p_value), .ibd_format_slope(res$slope)))
  cat("\n$pairs\n")
  if (nrow(x$pairs) <= 10) {
    .print_table(.round_table(x$pairs))
  } else {
    cat(sprintf("  (first 10 of %d pairs; the full table is x$pairs)\n", nrow(x$pairs)))
    .print_table(.round_table(utils::head(x$pairs, 10)))
  }
  cat("\nsummary() shows what to report and what to check; plot() draws the pairs.\n")
  invisible(x)
}

#' @rdname isolation_by_distance
#' @param details For `summary()`: `FALSE` (the default) prints a short view --
#'   the result, what to report and a list of checks. `TRUE` adds every pair
#'   and the notes on the method. Either way, the object returned holds the
#'   result table in `$tables$results`.
#' @export
summary.raddiv_ibd <- function(object, details = FALSE, ...) {
  .check_flag(details, "details")
  structure(list(result = object, details = details,
                 tables = list(results = .ibd_results(object))),
            class = "summary.raddiv_ibd")
}

#' @export
print.summary.raddiv_ibd <- function(x, ...) {
  .ibd_brief(x$result, details = isTRUE(x$details))
  invisible(x)
}

## Not exported. The RESULTS table: one row, the Mantel r and p, and the slope,
## which is a point estimate and has no SE (the header says so).
.ibd_results <- function(res) {
  r <- res$result
  tab <- data.frame(genetic = r$genetic, geographic = r$geographic, populations = r$populations,
                    pairs = r$pairs, stringsAsFactors = FALSE, check.names = FALSE)
  tab[["Mantel r"]] <- sprintf("%.3f", r$mantel_r)
  tab[["p"]] <- .format_p(r$p_value)
  tab[["orderings"]] <- .ibd_orderings_text(res$settings)
  tab[["slope (no SE)"]] <- .ibd_format_slope(r$slope)
  tab
}

## Not exported. The legend under the RESULTS table.
.ibd_legend <- function(st) {
  per <- if (st$habitat == "2D") "unit of ln(distance)" else "unit of distance"
  slope <- if (st$stat == "FST") {
    sprintf(paste("how fast FST/(1-FST) rises per %s; a point estimate with no SE.",
                  "1/slope estimates %s x density x sigma^2 (Rousset 1997; sigma^2 = mean",
                  "squared parent-offspring distance along one axis)."),
            per, if (st$habitat == "2D") "4 x pi" else "4")
  } else {
    sprintf("how fast D rises per %s; a point estimate with no SE. It describes the pattern; only the FST slope estimates dispersal.",
            per)
  }
  .legend_entries(c(
    "Mantel r" = paste(if (st$method == "pearson") "correlation" else "rank (Spearman) correlation",
                       "of the genetic and geographic values over all pairs (Mantel 1967)."),
    p = paste("one-sided: the share of orderings of the populations that give an r at least",
              "this large", if (st$exact) "(every ordering was tried, so p is exact)." else "(Phipson & Smyth 2010)."),
    slope = slope))
}

## Not exported. What to report.
.ibd_report_advice <- function(st) {
  if (st$stat == "FST")
    c("Mantel r and its one-sided p with the number of orderings, the habitat",
      sprintf("(FST/(1-FST) against %s), and the slope (no SE).",
              if (st$habitat == "2D") "ln(distance)" else "distance"),
      "Cite Mantel (1967) and Rousset (1997).")
  else
    c("Mantel r and its one-sided p with the number of orderings, and that D was",
      sprintf("compared with %s as it is. Cite Mantel (1967) and Jost (2008).",
              if (st$habitat == "2D") "ln(distance)" else "distance"))
}

## Not exported. The checks, as .check() lists.
.ibd_checks <- function(res) {
  st <- res$settings
  n <- length(st$populations)
  checks <- list(.check("ok", sprintf("All %d populations matched to the distances by name.", n)))
  if (length(st$excluded))
    checks[[length(checks) + 1L]] <- .check("info", paste("Left out with `exclude`:",
                                                          paste(st$excluded, collapse = ", ")))
  if (length(st$dropped_from_distances))
    checks[[length(checks) + 1L]] <- .check("info", paste(
      "In the distances but not in the genetic data (left out):",
      paste(st$dropped_from_distances, collapse = ", ")))
  if (st$min_p > 0.05) {
    checks[[length(checks) + 1L]] <- .check("look", sprintf(
      "With these %d populations and distances no ordering can give p below %s, so the test cannot detect isolation by distance. Sample more populations.",
      n, .format_p(st$min_p)))
  } else if (st$exact) {
    checks[[length(checks) + 1L]] <- .check("info", sprintf(
      "Every ordering of the %d populations was tried (%s: an exact test), so the smallest possible p is %s.",
      n, .big(st$orderings), .format_p(st$min_p)))
  }
  neg <- st$negative_pairs
  if (length(neg)) {
    what <- if (length(neg) == 1L) sprintf("1 pair has %s below 0", st$stat)
            else sprintf("%d pairs have %s below 0", length(neg), st$stat)
    checks[[length(checks) + 1L]] <- .check("info", sprintf(
      "%s (%s). This happens by chance when two populations barely differ. Such values are kept as they are, because setting them to 0 would bias the slope.",
      what, .shorten(neg, 3L)))
  }
  checks
}

## Not exported. The short summary, and with `details` the full one: every
## pair and the notes on the method.
.ibd_brief <- function(res, details = FALSE) {
  st <- res$settings
  cat(sprintf("ISOLATION BY DISTANCE: %d populations, %d pairs, %s habitat\n",
              res$result$populations, res$result$pairs, .ibd_habitat_text(st$habitat)))
  cat(sprintf("  from %s\n", st$source))
  .section("RESULTS")
  .print_compact(.ibd_results(res))
  .legend(.ibd_legend(st))
  .section("WHAT TO REPORT")
  .legend(.ibd_report_advice(st))
  .section("CHECKS")
  .check_lines(.ibd_checks(res))
  if (!details) {
    .section("Not shown here: every pair and the notes on the method.")
    .legend("summary(x, details = TRUE)   lists every pair and explains the method",
            "x$pairs, plot(x)              every pair, as a table and as a plot")
    return(invisible(NULL))
  }
  .section("EVERY PAIR")
  .print_table(.round_table(res$pairs))
  .section("NOTES")
  notes <- list(.report_text$ibd_mantel,
                if (st$habitat == "2D") .report_text$ibd_2d else .report_text$ibd_1d,
                if (st$stat == "D") .report_text$ibd_d else .report_text$ibd_slope,
                .report_text$ibd_partial)
  for (note in notes) {
    .legend(note)
    cat("\n")
  }
  invisible(NULL)
}

#' @rdname isolation_by_distance
#' @export
plot.raddiv_ibd <- function(x, ...) {
  st <- x$settings
  res <- x$result
  pr <- x$pairs
  defaults <- list(
    x = if (st$habitat == "2D") pr$log_distance else pr$distance,
    y = if (st$stat == "FST") pr$FST_linearized else pr$D,
    xlab = if (st$habitat == "2D") "ln(distance)" else "distance",
    ylab = if (st$stat == "FST") expression(F[ST] / (1 - F[ST])) else "Jost's D",
    main = sprintf("Isolation by distance (%s habitat)", .ibd_habitat_text(st$habitat)),
    pch = 16, col = "steelblue")
  do.call(graphics::plot, utils::modifyList(defaults, list(...)))
  graphics::abline(res$intercept, res$slope, col = "firebrick", lwd = 2)
  graphics::legend("topleft", bty = "n",
                   legend = c(sprintf("Mantel r = %.3f", res$mantel_r),
                              sprintf("p = %s, %s", .format_p(res$p_value),
                                      .ibd_orderings_text(st, sentence = TRUE))))
  invisible(x)
}
