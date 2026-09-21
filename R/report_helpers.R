###############################################################################
#
#  R/report_helpers.R -- shared building blocks for the summary() reports
#
#  summary() of a diversity_stats(), differentiation_stats(), het_between_pops()
#  or pi_allsites() result prints a SHORT view by default and the full report
#  with `details = TRUE`. The four short views follow the same rules, so the
#  pieces they share are here:
#
#    .est_se()         "estimate (SE)" text for one column of a table
#    .print_compact()  print a small table; if it is wider than the screen it
#                      wraps with the row label repeated on every block
#    .stat_blocks()    split a wide table into one small labelled table per
#                      statistic (Ho, Ho_se, Ho_lo, ... together)
#    .print_blocks()   print those blocks
#    .legend()         indented explanation lines under a table
#    .check_lines()    the CHECKS list: `ok`, `info` or `look`, one sentence each
#    .section()        a heading line
#
#  Nothing here is exported. Checked by tests/testthat/test-summary-views.R and
#  the summary tests of each function.
#
###############################################################################

## Not exported. "estimate (SE)" text for one column: `x` are the estimates and
## `se` their standard errors. One rule everywhere (screen, README recipe,
## manuscript table): the SE gets 2 significant digits and the estimate the
## same number of decimals, fixed within the column, so the numbers line up
## (0.3134 (0.0099); 0.00530 (0.00017)). The decimals come from the smallest
## usable SE in the column, rounded to 2 significant digits first (so 0.0000996
## counts as 0.00010 and gets 5 decimals, not 6); 3 when there is no usable SE.
## Missing values print as NA. The README's est_se() is the same rule.
.est_se <- function(x, se) {
  decimals <- .se_decimals(se)
  fmt <- function(v) ifelse(is.finite(v), sprintf("%.*f", decimals, v), "NA")
  paste0(fmt(x), " (", fmt(se), ")")
}

## Not exported. The number of decimals .est_se() uses for a column of SEs:
## enough for the smallest usable SE to show 2 significant digits (3 if none).
.se_decimals <- function(se) {
  usable <- is.finite(se) & se > 0
  if (any(usable)) as.integer(max(0, 1 - floor(log10(signif(min(se[usable]), 2))))) else 3L
}

## Not exported. A p-value for the short views: 2 significant digits, NA as "NA".
.format_p <- function(p) ifelse(is.finite(p), formatC(p, digits = 2, format = "g"), "NA")

## Not exported. A heading line for the short views: a blank line, then TEXT.
.section <- function(...) cat("\n", ..., "\n", sep = "")

## Not exported. Prints `df` (first column = the row label, e.g. population or
## pair) as it is, right-aligned. If it is wider than `width`, the other columns
## are packed into as many blocks as needed and the label column is repeated on
## every block, so no number is ever left without its row label.
.print_compact <- function(df, width = getOption("width")) {
  if (ncol(df) < 2L) {                 # a label and nothing else: show it, do not vanish
    print(df, row.names = FALSE)
    return(invisible(NULL))
  }
  ## Column widths as print() will show them. Numbers are formatted per column
  ## with a common number of digits (9.25e-05 prints as 0.0000925), so the width
  ## must come from format(), not from as.character(): a wrong guess lets print()
  ## wrap a block itself, and that wrapped part has no row label.
  shown <- format(df, na.encode = FALSE)
  cell_width <- function(j)
    max(nchar(names(df)[j], type = "width"), nchar(as.character(shown[[j]]), type = "width"),
        na.rm = TRUE) + 1L
  widths <- vapply(seq_along(df), cell_width, integer(1))
  start <- 2L
  while (start <= ncol(df)) {
    stop_at <- start
    used <- widths[1L] + widths[start]
    ## strictly less than the width: R wraps a matrix line that is exactly `width` wide
    while (stop_at < ncol(df) && used + widths[stop_at + 1L] < width) {
      stop_at <- stop_at + 1L
      used <- used + widths[stop_at]
    }
    print(df[, c(1L, start:stop_at), drop = FALSE], row.names = FALSE, right = TRUE)
    start <- stop_at + 1L
  }
  invisible(NULL)
}

## Not exported. Splits a wide table into one small table per statistic. Each
## block holds the label column (the first column of `df`) and every column of
## that statistic: for "Ho", the columns `Ho`, `Ho_se`, `Ho_lo`, ... A column
## belongs to a statistic when its name is the statistic or starts with it and
## a "_"; when several requested statistics match, the longest one gets it, so
## "Ar" never takes "privAr" columns, and "Ho" does not take "Ho_autosomal"
## columns when both are asked for. `always` names columns (such as the count
## `n`) that go into every block, right after the label. A statistic with no
## column in `df` gets no block, rather than a block holding only the label.
## Returns a named list of data frames.
.stat_blocks <- function(df, stats, always = character(0)) {
  label <- names(df)[1L]
  always <- intersect(always, names(df))
  owner <- vapply(names(df), function(col) {
    hits <- stats[vapply(stats, function(s) grepl(paste0("^", s, "(_|$)"), col), logical(1))]
    if (length(hits)) hits[which.max(nchar(hits))] else NA_character_
  }, character(1))
  stats <- stats[vapply(stats, function(s) any(owner == s, na.rm = TRUE), logical(1))]
  blocks <- lapply(stats, function(s)
    df[, c(label, always, setdiff(names(df)[!is.na(owner) & owner == s], always)), drop = FALSE])
  stats::setNames(blocks, stats)
}

## Not exported. Prints the blocks from .stat_blocks(), each followed by its
## note (a character vector, or NULL). `notes` is a named list keyed like
## `blocks`.
.print_blocks <- function(blocks, notes = list()) {
  for (nm in names(blocks)) {
    .print_compact(blocks[[nm]])      # wraps with the label repeated, if ever needed
    .legend(notes[[nm]])
  }
  invisible(NULL)
}

## Not exported. Explanation lines under a table, indented two spaces.
.legend <- function(...) {
  lines <- c(...)
  if (length(lines)) cat(paste0("  ", lines, "\n"), sep = "")
  invisible(NULL)
}

## Not exported. One check for the CHECKS list. `tag` is "ok" (nothing to do),
## "info" (a fact worth knowing) or "look" (act on it); `text` is one plain
## sentence saying what was found and, for "look", what to do.
.check <- function(tag, text) list(tag = tag, text = text)

## Not exported. Prints a list of .check() results as
##   ok    text
##   look  text, wrapped so the second line lines up under the first
.check_lines <- function(checks, width = 78L) {
  for (ch in checks) {
    if (is.null(ch)) next
    lines <- strwrap(ch$text, width = width - 8L)
    cat(paste0("  ", formatC(ch$tag, width = -6L), lines[1L], "\n"), sep = "")
    if (length(lines) > 1L) cat(paste0("        ", lines[-1L], "\n"), sep = "")
  }
  invisible(NULL)
}
