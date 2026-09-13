###############################################################################
#
#  R/utils.R -- small helpers shared by every file in the package.
#
#  Nothing here computes a statistic. These functions check arguments, print
#  progress messages, handle the random-number seed, round and write result
#  tables, and hold the two numerical tolerances used throughout. Keeping
#  them in one place means every exported function checks its arguments,
#  talks to the user and writes files in exactly the same way.
#
###############################################################################

## ---------------------------------------------------------------------------
## Numerical tolerances
## ---------------------------------------------------------------------------

## Not exported. Tolerance for comparing a computed fraction with a threshold
## the user typed. A call rate of 9 out of 10 individuals is stored as
## 0.8999999999999999 in floating point, so `0.9 >= 0.9` can come out FALSE;
## subtracting this tolerance from the threshold keeps a value sitting
## exactly ON the threshold on the "passes" side.
.threshold_tol <- 1e-9

## Not exported. Below this, a denominator is treated as zero (the statistic
## is undefined and reported as NA). A delete-one jackknife subtracts one
## locus's sums from the total, which can turn an exact 0 into 1e-17.
.zero_tol <- 1e-12

## ---------------------------------------------------------------------------
## Argument checks
##
## Each check stops with a message that names the argument, says what it
## should be, and shows what was actually passed. They check only the TYPE
## and LENGTH of an argument; limits that need an explanation (for example
## "g must be at most twice the smallest population") are checked by the
## function itself, with a message that explains the limit.
## ---------------------------------------------------------------------------

## Not exported. Shows a bad argument value in an error message: at most 5
## elements, so a mistaken whole-matrix argument does not flood the console.
.show_value <- function(x) {
  if (is.null(x)) return("NULL")
  shown <- paste(utils::head(format(x), 5), collapse = ", ")
  if (length(x) > 5) shown <- paste0(shown, ", ...")
  if (length(x) != 1L) shown <- paste0(shown, " (length ", length(x), ")")
  shown
}

## Not exported. `x` must be one whole number >= `min` (and <= `max`).
## Returns it as an integer, so `g = 20` and `g = 20L` behave identically.
.check_count <- function(x, name, min = 0, max = Inf) {
  ok <- is.numeric(x) && length(x) == 1L && !is.na(x) && is.finite(x) &&
    x == round(x) && x >= min && x <= max
  if (!ok)
    stop("`", name, "` must be a single whole number",
         if (is.finite(max)) paste0(" from ", min, " to ", max)
         else paste0(" >= ", min),
         " (got: ", .show_value(x), ").", call. = FALSE)
  as.integer(x)
}

## Not exported. `x` must be one number from `min` to `max` (inclusive).
## `Inf` is allowed only when `max = Inf`.
.check_number <- function(x, name, min = -Inf, max = Inf) {
  ok <- is.numeric(x) && length(x) == 1L && !is.na(x) && x >= min && x <= max
  if (!ok)
    stop("`", name, "` must be a single number",
         if (is.finite(min) && is.finite(max)) paste0(" from ", min, " to ", max)
         else if (is.finite(min)) paste0(" >= ", min)
         else if (is.finite(max)) paste0(" <= ", max),
         " (got: ", .show_value(x), ").", call. = FALSE)
  x
}

## Not exported. `x` must be exactly TRUE or FALSE.
.check_flag <- function(x, name) {
  if (!(is.logical(x) && length(x) == 1L && !is.na(x)))
    stop("`", name, "` must be TRUE or FALSE (got: ", .show_value(x), ").",
         call. = FALSE)
  x
}

## Not exported. `x` must be exactly one of `choices`. Deliberately NOT
## match.arg(): match.arg() accepts abbreviations, so a typo such as
## `boot = "individual"` would silently select a different mode instead of
## stopping.
.check_choice <- function(x, name, choices) {
  if (!(is.character(x) && length(x) == 1L && !is.na(x) && x %in% choices))
    stop("`", name, "` must be one of: ",
         paste0("\"", choices, "\"", collapse = ", "),
         " (got: ", .show_value(x), ").", call. = FALSE)
  x
}

## Not exported. `x` must be one non-empty character string (e.g. a path).
.check_string <- function(x, name) {
  if (!(is.character(x) && length(x) == 1L && !is.na(x) && nzchar(x)))
    stop("`", name, "` must be a single character string (got: ",
         .show_value(x), ").", call. = FALSE)
  x
}

## Not exported. `seed` must be NULL (use the session's random numbers) or
## one whole number.
.check_seed <- function(seed) {
  if (!is.null(seed)) .check_count(seed, "seed", min = -.Machine$integer.max)
  seed
}

## ---------------------------------------------------------------------------
## Messages
## ---------------------------------------------------------------------------

## Not exported. Every progress message in the package goes through this, so
## `verbose = FALSE` silences all of them. message() (not cat()) is used so
## that suppressMessages() works too.
.inform <- function(verbose, ...) {
  if (isTRUE(verbose)) message(...)
  invisible(NULL)
}

## Not exported. A whole number with thousands separators: 12345 -> "12,345".
.big <- function(x) format(x, big.mark = ",", scientific = FALSE, trim = TRUE)

## ---------------------------------------------------------------------------
## Random numbers
## ---------------------------------------------------------------------------

## Not exported. Captures the caller's current random-number state so it can
## be put back with on.exit(). Used only when a function is given an explicit
## `seed`:
##
##   if (!is.null(seed)) {
##     restore_rng <- .save_rng_state()
##     on.exit(restore_rng(), add = TRUE)
##     set.seed(seed)
##   }
##
## so that passing `seed` makes a result reproducible WITHOUT changing the
## random numbers the user's own code gets afterwards. With `seed = NULL` (the
## default everywhere) nothing is saved or restored: the function simply uses
## the session's random numbers, and `set.seed()` before the call makes the
## result reproducible, as in base R.
.save_rng_state <- function() {
  if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
    old <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    function() assign(".Random.seed", old, envir = .GlobalEnv)
  } else {
    ## No random number had been drawn yet this session: put back that same
    ## "no state" rather than leaving a .Random.seed that was not there.
    function() {
      if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
        rm(".Random.seed", envir = .GlobalEnv)
    }
  }
}

## ---------------------------------------------------------------------------
## Optional packages
## ---------------------------------------------------------------------------

## Not exported. TRUE when hierfstat is installed. A wrapper rather than a
## direct requireNamespace() call so that tests can pretend hierfstat is
## missing (testthat::with_mocked_bindings() can replace a function defined
## in this package, but not one from base R).
.hierfstat_available <- function() requireNamespace("hierfstat", quietly = TRUE)

## ---------------------------------------------------------------------------
## Result tables: rounding and writing
## ---------------------------------------------------------------------------

## Not exported. Rounds the numeric (double) columns of a data frame for
## display or for a TSV file. Result objects themselves always keep full
## precision; only print() and the file writer call this.
##   digits       decimal places for every double column not listed below
##   round_cols   named vector: decimal places for particular columns
##   signif_cols  named vector: significant digits for particular columns
## Integer, character and logical columns are never changed.
.round_table <- function(df, digits = 4, round_cols = NULL, signif_cols = NULL) {
  if (is.null(df)) return(NULL)
  for (col in names(df)) {
    x <- df[[col]]
    if (!is.double(x)) next
    df[[col]] <- if (col %in% names(signif_cols)) signif(x, signif_cols[[col]])
                 else if (col %in% names(round_cols)) round(x, round_cols[[col]])
                 else round(x, digits)
  }
  df
}

## Not exported. Writes each non-NULL table in `tables` (a named list of
## data frames, already rounded) to `outdir` as a tab-separated file.
## `file_names` is a named character vector with the same names as `tables`.
## Creates `outdir` if needed and returns the paths written.
.write_tables <- function(tables, outdir, file_names) {
  if (!dir.exists(outdir) && !dir.create(outdir, recursive = TRUE, showWarnings = FALSE))
    stop("Could not create the output directory: ", outdir, call. = FALSE)
  written <- character(0)
  for (nm in names(tables)) {
    if (is.null(tables[[nm]])) next
    path <- file.path(outdir, file_names[[nm]])
    utils::write.table(tables[[nm]], path, sep = "\t", quote = FALSE, row.names = FALSE)
    written <- c(written, path)
  }
  written
}

## Not exported. Prints a data frame without row numbers. Callers round it
## first (with .round_table() or a function built on it).
.print_table <- function(df) {
  print(df, row.names = FALSE)
  invisible(NULL)
}
