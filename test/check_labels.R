#!/usr/bin/env Rscript
## Guard against the bug the terminology pass exposed: the script's prose
## naming a column that the script does not print. Runs both file types and
## checks that every capitalised identifier mentioned in the output text also
## appears as a printed column header.

files <- c("sim.allsnps.vcf", "sim.haps.vcf")
## Pre-existing bug (now historical): this path once had to be plain
## "diversity_stats.R", because check_labels.R and diversity_stats.R both
## lived in the repo root and every invocation ran from there. As of the
## test/ reorganization, check_labels.R lives in test/ while
## diversity_stats.R stays in the repo root one directory up, and
## run_tests.sh always cd's into test/ before invoking this script -- so
## "../diversity_stats.R" is now the correct path, not the bug it used to
## be. If this script or diversity_stats.R ever moves again, re-check this
## the same way the comment above once had to: a silently-vacuous PASS (0
## columns, 0 mentions checked) is what a wrong path looks like here, not
## an error.
SCRIPT <- "../diversity_stats.R"
args  <- c("popmap.tsv", "--nboot=0", "--sites=1000000", "--g=20")
ok <- TRUE

for (f in files) {
  out <- system(paste(paste("Rscript", SCRIPT), f, paste(args, collapse = " "),
                      "2>/dev/null"), intern = TRUE)
  txt <- paste(out, collapse = "\n")

  ## column headers = tokens on lines that begin a printed data.frame
  hdr <- grep("^ *population ", out, value = TRUE)
  cols <- unique(unlist(strsplit(trimws(hdr), " +")))

  ## identifiers referred to in prose: Ho_x, He_x, Fis_x, pct_x style tokens
  ment <- unique(unlist(regmatches(txt,
    gregexpr("\\b(?:Ho|He|Fis|Ar|priv|pct)[A-Za-z0-9_]*\\b", txt))))
  ment <- setdiff(ment, c("Ho", "He", "Fis", "Ar", "priv", "pct"))
  ## words that are prose, not identifiers
  ment <- ment[grepl("_", ment)]
  ## A mention that explicitly refers to the OTHER run ("re-run on ... to get
  ## He_autosomal", "the conversion was SKIPPED") is not a promise of output.
  crosstxt <- paste(grep("SKIPPED|Re-run|Run populations|TAKE FROM|DO NOT take",
                         out, value = TRUE), collapse = "\n")
  cross <- unique(unlist(regmatches(crosstxt,
    gregexpr("\\b(?:Ho|He|Fis|Ar|priv|pct)[A-Za-z0-9_]*\\b", crosstxt))))
  ment <- setdiff(ment, cross)

  missing <- setdiff(ment, cols)
  cat(sprintf("%-18s columns: %2d   mentioned: %2d   ", f, length(cols), length(ment)))
  if (length(missing)) {
    cat("MISSING FROM OUTPUT: ", paste(missing, collapse = ", "), "\n", sep = "")
    ok <- FALSE
  } else cat("all mentioned identifiers are printed  OK\n")

  ## and the reverse: a column nobody explains
  if ("pi" %in% cols) { cat("  a column is literally named 'pi'\n"); ok <- FALSE }
}

## the script must never claim to print a column called pi
for (f in files) {
  out <- system(paste(paste("Rscript", SCRIPT), f, paste(args, collapse = " "),
                      "2>/dev/null"), intern = TRUE)
  bad <- grep("TAKE FROM THIS RUN.*\\bpi\\b|Ho, He, pi", out, value = TRUE)
  if (length(bad)) { cat("  prose promises a 'pi' column: ", bad[1], "\n"); ok <- FALSE }
}

cat(if (ok) "\nPASS: no column named in prose is missing from the output\n"
    else "\nFAIL\n")
quit(status = if (ok) 0 else 1)
