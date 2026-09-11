###############################################################################
#
#  R/read_sumstats_summary.R -- parser for Stacks' populations.sumstats_summary.tsv
#
#  WHAT THIS IS FOR. This package's own Ho, He and FIS estimators
#  deliberately differ from Stacks' (see vignette("rationale"), "Formulas": Nei-Chesser
#  He instead of Stacks' `Pi`, which assumes FIS = 0; FIS as a ratio of sums
#  instead of Stacks' mean of per-locus ratios). NOTHING in this package
#  treats sumstats_summary.tsv as a source of Ho, He, pi or FIS -- diversity_
#  stats() computes all of those itself. The one column in this file that is
#  a plain tally rather than an estimator is `Sites`, in the "All positions
#  (variant and fixed)" block: how many nucleotide positions Stacks retained
#  for that population. That count is genuinely useful -- it is the
#  sequenced-site denominator diversity_stats()'s `sites` argument needs for
#  the autosomal Ho/He (nucleotide diversity) conversion, and it differs
#  slightly per population for the same reason a locus can be missing from
#  one population and not another. read_sumstats_summary() below is a
#  general parser for the whole file (both blocks, every column) so that
#  number can be pulled out by name instead of retyped by hand -- see
#  diversity_stats()'s `sites` argument for the intended use. Nothing else in
#  this package calls it.
#
#  FILE SHAPE. Two blocks, always in this order:
#    # Variant positions
#    # Pop ID<TAB>Private<TAB>Num_Indv<TAB>Var<TAB>StdErr<TAB>P<TAB>Var<TAB>StdErr<TAB>...
#    <one data row per population>
#    # All positions (variant and fixed)
#    # Pop ID<TAB>Private<TAB>Sites<TAB>Variant_Sites<TAB>Polymorphic_Sites<TAB>%Polymorphic_Loci<TAB>Num_Indv<TAB>Var<TAB>StdErr<TAB>...
#    <one data row per population>
#  `Var`/`StdErr` repeat after every one of 8 statistics (Num_Indv, P,
#  Obs_Het, Obs_Hom, Exp_Het, Exp_Hom, Pi, Fis), and `%Polymorphic_Loci`
#  isn't a syntactic name -- this header cannot be read with
#  read.delim()/column names. It is matched POSITIONALLY against the known
#  token sequence below (the same idiom read_stacks_vcf() uses for the
#  GT-first-FORMAT check, R/vcf_io.R), and disambiguated into `<stat>`,
#  `<stat>_var`, `<stat>_se`. A header that doesn't match -- e.g. a Stacks
#  version that changed this file's columns -- is a loud, found-vs-expected
#  error rather than a silent misalignment.
#
###############################################################################

## Not exported. The 8 statistics that repeat as (name, Var, StdErr) triples
## in both blocks, in the fixed order Stacks writes them, and their
## R-friendly equivalents.
.sumstats_triple_stats <- c("Num_Indv", "P", "Obs_Het", "Obs_Hom",
                            "Exp_Het", "Exp_Hom", "Pi", "Fis")
.sumstats_triple_names <- c("num_indv", "p", "obs_het", "obs_hom",
                            "exp_het", "exp_hom", "pi", "fis")

## Not exported. The full expected header token sequence for one block:
## `lead_cols` first (e.g. "Pop ID", "Private"), then the 8 (stat, Var,
## StdErr) triples.
.expected_sumstats_header <- function(lead_cols) {
  triples <- as.vector(rbind(.sumstats_triple_stats, "Var", "StdErr"))
  c(lead_cols, triples)
}

## Not exported. Turns one lead-column header token into an R-friendly name:
## "%Polymorphic_Loci" -> "pct_polymorphic_loci", "Variant_Sites" ->
## "variant_sites".
.sumstats_lead_name <- function(tok) {
  tok <- sub("^%", "pct_", tok)
  tok <- tolower(gsub("[^A-Za-z0-9]+", "_", tok))
  sub("^_+", "", sub("_+$", "", tok))
}

## Not exported. Parses one block into a data frame, one row per population.
## `header`/`rows` are already tab-split; `lead_cols` are the block's own
## non-repeating leading columns (always starts "Pop ID", "Private").
.parse_sumstats_block <- function(block_name, header, rows, lead_cols) {
  expected <- .expected_sumstats_header(lead_cols)
  if (!identical(header, expected))
    stop(sprintf(
      "populations.sumstats_summary.tsv: the \"%s\" block's header does not ",
      block_name),
      "match the expected Stacks column layout -- this parser reads columns ",
      "POSITIONALLY (Var/StdErr repeat and cannot be matched by name), so a ",
      "different layout cannot be read safely.\n  Expected: ",
      paste(expected, collapse = ", "), "\n  Found:    ",
      paste(header, collapse = ", "),
      "\n  This usually means a different Stacks version changed this ",
      "file's columns; check the Stacks changelog for ",
      "`populations.sumstats_summary.tsv`.")

  n_lead <- length(lead_cols)
  n_col  <- length(expected)
  bad_len <- which(lengths(rows) != n_col)
  if (length(bad_len))
    stop(sprintf(
      "populations.sumstats_summary.tsv: the \"%s\" block's row %d has %d ",
      block_name, bad_len[1], lengths(rows)[bad_len[1]]),
      "tab-separated field(s), expected ", n_col, " to match its header. ",
      "The file is truncated or malformed.")

  m <- do.call(rbind, rows)
  out <- data.frame(population = m[, 1L], stringsAsFactors = FALSE)
  ## "Private" and any block-specific extra lead columns (Sites,
  ## Variant_Sites, Polymorphic_Sites, %Polymorphic_Loci) are single,
  ## non-repeating columns, carried straight across.
  if (n_lead > 1L) {
    for (k in 2:n_lead)
      out[[.sumstats_lead_name(lead_cols[k])]] <-
        suppressWarnings(as.numeric(m[, k]))
  }
  ## The 8 repeating triples: value, Var, StdErr, in fixed order.
  for (i in seq_along(.sumstats_triple_names)) {
    base <- n_lead + (i - 1L) * 3L
    nm <- .sumstats_triple_names[i]
    out[[nm]]                 <- suppressWarnings(as.numeric(m[, base + 1L]))
    out[[paste0(nm, "_var")]] <- suppressWarnings(as.numeric(m[, base + 2L]))
    out[[paste0(nm, "_se")]]  <- suppressWarnings(as.numeric(m[, base + 3L]))
  }
  out
}

#' Read a Stacks `populations.sumstats_summary.tsv` file
#'
#' Parses both blocks of Stacks' `populations.sumstats_summary.tsv`
#' ("Variant positions" and "All positions (variant and fixed)") into one
#' data frame per block, one row per population.
#'
#' This package's own Ho, He and FIS estimators deliberately differ from
#' Stacks' -- see `vignette("rationale")` ("Formulas"): Nei & Chesser's He instead of
#' Stacks' `Pi` (which assumes FIS = 0), and FIS as a ratio of sums instead
#' of Stacks' mean of per-locus ratios. **Nothing in this package uses the
#' values this function returns as an estimate of Ho, He, pi or FIS** --
#' [diversity_stats()] computes all of those itself. This function exists so
#' the one column here that is a plain tally rather than an estimator --
#' `sites` (how many nucleotide positions Stacks retained for that
#' population, in the "All positions" block) -- can be read by name and
#' handed to [diversity_stats()]'s `sites` argument instead of retyped by
#' hand. `sites` differs slightly per population for the same reason a locus
#' can be missing from one population and not another.
#'
#' @param path Path to `populations.sumstats_summary.tsv`.
#' @return A list with elements `variant_positions` and `all_positions`,
#'   each a data frame with one row per population and columns:
#'   `population`, `private`, (in `all_positions` only:) `sites`,
#'   `variant_sites`, `polymorphic_sites`, `pct_polymorphic_loci`, then for
#'   each of `num_indv`, `p`, `obs_het`, `obs_hom`, `exp_het`, `exp_hom`,
#'   `pi`, `fis`: the value itself plus `<name>_var` and `<name>_se`.
#' @examples
#' triple_hdr <- paste(
#'   "Num_Indv", "Var", "StdErr", "P", "Var", "StdErr", "Obs_Het", "Var",
#'   "StdErr", "Obs_Hom", "Var", "StdErr", "Exp_Het", "Var", "StdErr",
#'   "Exp_Hom", "Var", "StdErr", "Pi", "Var", "StdErr", "Fis", "Var",
#'   "StdErr", sep = "\t")
#' lines <- c(
#'   "# Variant positions",
#'   paste0("# Pop ID\tPrivate\t", triple_hdr),
#'   paste("popA", 0, 10, 0, 0, 0.9, 0, 0, 0.2, 0, 0, 0.8, 0, 0, 0.18, 0, 0,
#'         0.82, 0, 0, 0.19, 0, 0, -0.05, 0, 0, sep = "\t"),
#'   "# All positions (variant and fixed)",
#'   paste0("# Pop ID\tPrivate\tSites\tVariant_Sites\tPolymorphic_Sites\t",
#'          "%Polymorphic_Loci\t", triple_hdr),
#'   paste("popA", 0, 100000, 10, 4, 40, 10, 0, 0, 0.9, 0, 0, 0.00002, 0, 0,
#'         0.99998, 0, 0, 0.000018, 0, 0, 0.999982, 0, 0, 0.000019, 0, 0,
#'         -0.05, 0, 0, sep = "\t")
#' )
#' f <- tempfile(fileext = ".tsv")
#' writeLines(lines, f)
#' read_sumstats_summary(f)$all_positions$sites  # 1e5
#' @export
read_sumstats_summary <- function(path) {
  if (!file.exists(path))
    stop("File not found: ", path, "\n  Check the path and try again.")
  lines <- readLines(path)
  lines <- lines[nzchar(lines)]

  ## A block title line is "#"-prefixed with no tab in it; a header line is
  ## also "#"-prefixed but DOES have tabs, which is what tells the two apart.
  title_i <- grep("^#[^\t]*$", lines)
  if (length(title_i) != 2L)
    stop("populations.sumstats_summary.tsv: expected exactly 2 block title ",
         "lines (\"# Variant positions\" and \"# All positions (variant ",
         "and fixed)\"), found ", length(title_i),
         ". Is this really a Stacks sumstats_summary file?")
  titles <- sub("^#\\s*", "", lines[title_i])
  if (!grepl("^Variant positions", titles[1]) ||
      !grepl("^All positions", titles[2]))
    stop("populations.sumstats_summary.tsv: expected the blocks in the ",
         "order \"Variant positions\" then \"All positions (variant and ",
         "fixed)\", found: \"", titles[1], "\" then \"", titles[2], "\".")

  block_bounds <- c(title_i, length(lines) + 1L)
  parse_one <- function(b, lead_cols) {
    start <- block_bounds[b]; end <- block_bounds[b + 1L] - 1L
    header <- strsplit(sub("^#\\s*", "", lines[start + 1L]), "\t", fixed = TRUE)[[1]]
    body   <- lines[(start + 2L):end]
    rows   <- strsplit(body, "\t", fixed = TRUE)
    .parse_sumstats_block(titles[b], header, rows, lead_cols)
  }

  list(
    variant_positions = parse_one(1L, c("Pop ID", "Private")),
    all_positions     = parse_one(2L, c("Pop ID", "Private", "Sites",
                                        "Variant_Sites", "Polymorphic_Sites",
                                        "%Polymorphic_Loci"))
  )
}
