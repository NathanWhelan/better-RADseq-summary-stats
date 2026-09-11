###############################################################################
#
#  R/cli.R -- the command-line entry points.
#
#  The scripts in inst/scripts/ (installed to
#  system.file("scripts", package = "RADdiversity")) are one-line shims that
#  call .cli_main(). All argument parsing lives here, inside the package, so
#  every script parses flags the same way and the parsing is tested with
#  everything else (tests/testthat/test-cli.R):
#
#    Rscript diversity_stats.R  <vcf> <popmap.tsv> --g=N [--nboot=N] ...
#    Rscript het_between_pops.R <vcf> <popmap.tsv> [--min-call=X] [--outdir=DIR]
#    Rscript het_between_pops.R --selftest
#    Rscript diversity_core.R   --selftest
#
#  Nothing here calls quit(): .cli_main() returns the exit status and the
#  shim quits with it, so these functions stay safe to call from an
#  interactive R session.
#
###############################################################################

.cli_usage <- list(
  diversity_stats = c(
    "Usage: Rscript diversity_stats.R <vcf> <popmap.tsv> --g=N [options]",
    "",
    "  Works on populations.snps.vcf or populations.haps.vcf (optionally .gz).",
    "  Run it on both: the SNP VCF gives Ho, He, pct_poly and the per-site",
    "  values; the haplotype VCF gives Fis, Ar and privAr.",
    "",
    "  --g=N            REQUIRED. Rarefaction size in GENE COPIES (10 diploids",
    "                   = 20), at most twice the smallest population.",
    "  --nboot=N        bootstrap replicates (default 10000; 0 = none)",
    "  --boot=MODE      loci | individuals | both (default loci; the other two",
    "                   are for comparison only -- see ?diversity_stats)",
    "  --sites=N|FILE   sequenced sites, for per-site Ho/He: one number, or the",
    "                   path to populations.sumstats_summary.tsv (each",
    "                   population's own Sites is then read from it)",
    "  --min-n=N        minimum typed individuals per population per locus",
    "                   (default 2)",
    "  --outdir=DIR     where to write the TSV files (default: current directory)",
    "  --complete-case  use a record only if every individual of every",
    "                   population is genotyped there (Schmidt et al. 2021)"),
  het_between_pops = c(
    "Usage: Rscript het_between_pops.R <vcf> <popmap.tsv> [--min-call=X] [--outdir=DIR]",
    "       Rscript het_between_pops.R <vcf> <popmap.tsv> [min_call] [outdir]",
    "       Rscript het_between_pops.R --selftest",
    "",
    "  Tests whether heterozygosity differs between populations, using the",
    "  INDIVIDUAL as the unit of replication. Needs a VCF: per-individual",
    "  heterozygosity cannot be recovered from populations.sumstats.tsv.",
    "",
    "  --min-call=X   minimum fraction of individuals genotyped at a locus,",
    "                 within each population, for the locus to be used",
    "                 (default 0.9)",
    "  --outdir=DIR   where to write the TSV files (default: current directory)"),
  diversity_core = c(
    "Usage: Rscript diversity_core.R --selftest",
    "",
    "  Checks every estimator against brute-force Monte Carlo and against",
    "  published Stacks output.")
)

## Not exported. Splits command-line arguments into --name=value flags (only
## names in `flag_names`), bare --switches (only names in `switches`) and
## positional arguments. Anything else is an error that names the fix.
.cli_parse <- function(args, flag_names = character(), switches = character()) {
  flags <- list(); on <- character(); positional <- character()
  for (a in args) {
    nm <- sub("^--", "", a)
    m <- regmatches(a, regexec("^--([a-z-]+)=(.*)$", a))[[1]]
    if (length(m) == 3L) {
      if (!(m[2] %in% flag_names))
        stop("Unrecognized flag: --", m[2], "=... . Recognized: ",
             paste0("--", c(paste0(flag_names, "="), switches), collapse = ", "),
             call. = FALSE)
      flags[[m[2]]] <- m[3]
    } else if (startsWith(a, "--") && nm %in% switches) {
      on <- c(on, nm)
    } else if (startsWith(a, "--") && nm %in% flag_names) {
      ## e.g. `--g 20`, a common slip for `--g=20`: name the fix.
      stop("--", nm, " requires a value: --", nm, "=VALUE", call. = FALSE)
    } else if (startsWith(a, "--")) {
      stop("Unrecognized flag: ", a, call. = FALSE)
    } else {
      positional <- c(positional, a)
    }
  }
  list(flags = flags, switches = on, positional = positional)
}

## Not exported. `--sites=` is one number or a path to
## populations.sumstats_summary.tsv. Command-line values always arrive as
## text, and diversity_stats() reads any text as a file path, so a number
## must be converted here first (without this, --sites=123456 failed with
## "file not found: 123456").
.cli_sites <- function(x) {
  if (is.null(x)) return(0)
  num <- suppressWarnings(as.numeric(x))
  if (is.na(num)) x else num
}

## Not exported. Runs one command-line command and returns its exit status
## (0 = success) instead of calling quit(). An error becomes a one-line
## "Error: ..." message and status 1.
.cli_main <- function(cmd, args) {
  usage <- function() cat(.cli_usage[[cmd]], sep = "\n")
  run <- switch(cmd,
    diversity_stats = function() {
      p <- .cli_parse(args, c("g", "nboot", "boot", "sites", "min-n", "outdir"),
                      "complete-case")
      if (length(p$positional) != 2L) { usage(); return(1L) }
      f <- p$flags
      if (is.null(f$g)) {
        usage()
        stop("Missing required --g=N (rarefaction size in gene copies).", call. = FALSE)
      }
      print(diversity_stats(p$positional[1], p$positional[2], g = f$g,
                            nboot  = if (is.null(f$nboot)) 10000L else f$nboot,
                            boot   = if (is.null(f$boot)) "loci" else f$boot,
                            sites  = .cli_sites(f$sites),
                            min_n  = if (is.null(f[["min-n"]])) 2L else f[["min-n"]],
                            complete_case = "complete-case" %in% p$switches,
                            outdir = if (is.null(f$outdir)) "." else f$outdir))
      0L
    },
    het_between_pops = function() {
      p <- .cli_parse(args, c("min-call", "outdir"), "selftest")
      if ("selftest" %in% p$switches) { het_between_pops_selftest(); return(0L) }
      pos <- p$positional
      if (length(pos) < 2L || length(pos) > 4L) { usage(); return(1L) }
      ## Flags win; the older positional form (<vcf> <popmap> [min_call]
      ## [outdir]) still works.
      min_call <- p$flags[["min-call"]]
      if (is.null(min_call)) min_call <- if (length(pos) >= 3L) pos[3] else 0.9
      outdir <- p$flags$outdir
      if (is.null(outdir)) outdir <- if (length(pos) >= 4L) pos[4] else "."
      print(het_between_pops(pos[1], pos[2], min_call = min_call, outdir = outdir))
      0L
    },
    diversity_core = function() {
      p <- .cli_parse(args, switches = "selftest")
      if (!("selftest" %in% p$switches)) { usage(); return(1L) }
      if (isTRUE(diversity_core_selftest())) 0L else 1L
    },
    stop("Unknown command: ", cmd, call. = FALSE))
  tryCatch(run(), error = function(e) {
    message("Error: ", conditionMessage(e))
    1L
  })
}
