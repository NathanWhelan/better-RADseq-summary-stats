###############################################################################
#
#  inst/sims/make_example_data.R -- builds the example dataset in inst/extdata.
#
#  SIMULATED data, shaped like Stacks 2 `populations` output, for the
#  vignettes: large enough to show realistic output, small enough to ship.
#
#    example.snps.vcf.gz          one record per SNP (FORMAT GT:DP:AD)
#    example.haps.vcf.gz          one record per RAD locus (FORMAT GT; ID "."
#                                 and POS 0, as in a Stacks de novo run)
#    example.allsites.vcf.gz      every sequenced site, variable or not (GT)
#    example.sumstats_summary.tsv the populations.sumstats_summary.tsv layout
#    example_popmap.tsv           sample <TAB> population
#
#  What is built in, so the vignettes have something to find:
#    * two populations, "north" (15 individuals) and "south" (10), with south
#      carrying fewer haplotypes (it went through a bottleneck);
#    * 1,000 RAD loci of 80 bp with 0-4 SNPs each (loci with no SNP appear
#      only in the all-sites VCF);
#    * every individual inbred to its own degree: F drawn around 0.1 with SD
#      0.15, so individuals differ in inbreeding (g2 > 0). Identity by descent
#      is drawn once per individual per RAD locus, as real inbreeding acts on
#      whole haplotypes;
#    * one pair of full siblings, north_01 and north_02;
#    * missing data clustered by locus, plus genotypes backed by fewer than 4
#      reads set to missing, as a genotype caller would.
#
#  In the sumstats file only `Sites` and `Variant_Sites` matter to this
#  package; the other columns follow Stacks' definitions but are computed
#  simply (per-site averages over the simulated genotypes).
#
#  Run from the package source directory:  Rscript inst/sims/make_example_data.R
#
###############################################################################

set.seed(20260913)
out_dir <- file.path("inst", "extdata")
n_north <- 15L
n_south <- 10L
n_tags <- 1000L
tag_len <- 80L
samples <- c(sprintf("north_%02d", seq_len(n_north)), sprintf("south_%02d", seq_len(n_south)))
pop <- rep(c("north", "south"), c(n_north, n_south))
n <- length(samples)
F_ind <- pmin(pmax(stats::rnorm(n, 0.1, 0.15), 0), 0.8)
bases <- c("A", "C", "G", "T")

snp_rows <- hap_rows <- all_rows <- list()
stats_rows <- list()          # per SNP record: counts for the sumstats file
typed_tags <- matrix(FALSE, n_tags, 2, dimnames = list(NULL, c("north", "south")))

vcf_header <- function(format_note) c(
  "##fileformat=VCFv4.2",
  "##source=\"RADdiversity example data (simulated; inst/sims/make_example_data.R)\"",
  format_note,
  paste(c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", samples),
        collapse = "\t"))

for (t in seq_len(n_tags)) {
  k <- sample(0:4, 1L, prob = c(0.25, 0.35, 0.22, 0.12, 0.06))
  positions <- sort(sample(10:(tag_len - 5L), k))
  ref_base <- sample(bases, tag_len, replace = TRUE)

  ## Haplotypes (rows) at the k SNPs (0 = REF base, 1 = ALT base).
  if (k > 0L) {
    repeat {
      n_hap <- min(2L^k, sample(2:6, 1L))
      haps <- unique(matrix(stats::rbinom(n_hap * k, 1L, 0.5), ncol = k))
      if (nrow(haps) >= 2L && all(colSums(haps) > 0 & colSums(haps) < nrow(haps))) break
    }
    alt_base <- vapply(positions, function(p) sample(setdiff(bases, ref_base[p]), 1L), "")
  } else {
    haps <- matrix(0L, 1L, 0L)
  }
  n_hap <- nrow(haps)
  f_north <- stats::rgamma(n_hap, 0.8)
  f_north <- f_north / sum(f_north)
  f_south <- stats::rgamma(n_hap, 6 * f_north + 0.05)       # bottleneck: drift
  f_south <- f_south / sum(f_south)

  draw <- function(freq) sample.int(n_hap, 1L, prob = freq)
  h1 <- h2 <- integer(n)
  for (i in seq_len(n)) {
    freq <- if (pop[i] == "north") f_north else f_south
    h1[i] <- draw(freq)
    h2[i] <- if (stats::runif(1) < F_ind[i]) h1[i] else draw(freq)
  }
  ## Full siblings north_01 and north_02: one haplotype from each parent.
  mother <- c(draw(f_north), draw(f_north))
  father <- c(draw(f_north), draw(f_north))
  for (i in 1:2) {
    h1[i] <- mother[sample.int(2L, 1L)]
    h2[i] <- father[sample.int(2L, 1L)]
  }

  ## Missing data: clustered by locus, and calls on fewer than 4 reads.
  depth <- stats::rnbinom(n, mu = 22, size = 3) + 1L
  missing <- stats::runif(n) < stats::rbeta(1L, 0.6, 9) | depth < 4L
  typed_tags[t, "north"] <- any(!missing[pop == "north"])
  typed_tags[t, "south"] <- any(!missing[pop == "south"])

  ## SNP records: only SNPs variable among the typed genotypes, as Stacks writes.
  variable <- logical(k)
  snp_gt <- matrix("./.", k, n)
  for (s in seq_len(k)) {
    a <- haps[h1, s]
    b <- haps[h2, s]
    variable[s] <- length(unique(c(a[!missing], b[!missing]))) == 2L
    if (!variable[s]) next
    alt_reads <- ifelse(a + b == 0L, 0L, ifelse(a + b == 2L, depth, stats::rbinom(n, depth, 0.5)))
    cell <- paste0(pmin(a, b), "/", pmax(a, b), ":", depth, ":", depth - alt_reads, ",", alt_reads)
    cell[missing] <- "./."
    snp_gt[s, ] <- cell
    snp_rows[[length(snp_rows) + 1L]] <- paste(
      t, positions[s], paste0(t, ":", positions[s], ":+"), ref_base[positions[s]], alt_base[s],
      ".", "PASS", ".", "GT:DP:AD", paste(cell, collapse = "\t"), sep = "\t")
    stats_rows[[length(stats_rows) + 1L]] <- list(a = a, b = b, missing = missing)
  }

  ## Haplotype record: the variable SNPs' bases, one allele per distinct string.
  if (any(variable)) {
    v <- which(variable)
    hap_string <- apply(haps[, v, drop = FALSE], 1L, function(r)
      paste(ifelse(r == 1L, alt_base[v], ref_base[positions[v]]), collapse = ""))
    typed_strings <- c(hap_string[h1[!missing]], hap_string[h2[!missing]])
    alleles <- names(sort(table(typed_strings), decreasing = TRUE))
    code1 <- match(hap_string[h1], alleles) - 1L
    code2 <- match(hap_string[h2], alleles) - 1L
    cell <- paste0(pmin(code1, code2), "/", pmax(code1, code2))
    cell[missing] <- "./."
    ## As Stacks 2 writes a de novo haplotype VCF: CHROM is the locus, POS 0,
    ## ID "." (export_formats.cc, VcfHapsExport::write_batch()).
    hap_rows[[length(hap_rows) + 1L]] <- paste(
      t, 0L, ".", alleles[1], if (length(alleles) > 1L) paste(alleles[-1], collapse = ",") else ".",
      ".", "PASS", ".", "GT", paste(cell, collapse = "\t"), sep = "\t")
  }

  ## All-sites records: every base of the locus.
  fixed_gt <- ifelse(missing, "./.", "0/0")
  site_lines <- character(tag_len)
  for (p in seq_len(tag_len)) {
    s <- match(p, positions)
    if (!is.na(s) && variable[s]) {
      gt <- sub(":.*$", "", snp_gt[s, ])
      site_lines[p] <- paste(t, p, ".", ref_base[p], alt_base[s], ".", "PASS", ".", "GT",
                             paste(gt, collapse = "\t"), sep = "\t")
    } else {
      site_lines[p] <- paste(t, p, ".", ref_base[p], ".", ".", "PASS", ".", "GT",
                             paste(fixed_gt, collapse = "\t"), sep = "\t")
    }
  }
  all_rows[[t]] <- site_lines
}

write_gz <- function(lines, file) {
  con <- gzfile(file.path(out_dir, file), "w")
  on.exit(close(con))
  writeLines(lines, con)
}
write_gz(c(vcf_header("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">"),
           unlist(snp_rows)), "example.snps.vcf.gz")
write_gz(c(vcf_header("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">"),
           unlist(hap_rows)), "example.haps.vcf.gz")
write_gz(c(vcf_header("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">"),
           unlist(all_rows)), "example.allsites.vcf.gz")
writeLines(paste(samples, pop, sep = "\t"), file.path(out_dir, "example_popmap.tsv"))

## populations.sumstats_summary.tsv layout. Per SNP record and population:
## N typed, P (commoner allele's frequency), Obs_Het, Exp_Het = 1 - sum p^2,
## Pi = 2N/(2N - 1) Exp_Het, Fis = (Pi - Obs_Het)/Pi (0 where Pi = 0).
site_stats <- function(r, members) {
  typed <- members & !r$missing
  N <- sum(typed)
  if (N == 0L) return(NULL)
  p_alt <- (sum(r$a[typed]) + sum(r$b[typed])) / (2 * N)
  P <- max(p_alt, 1 - p_alt)
  ho <- mean(r$a[typed] != r$b[typed])
  he <- 1 - p_alt^2 - (1 - p_alt)^2
  pi <- if (N > 0) he * 2 * N / (2 * N - 1) else 0
  c(N = N, P = P, Obs_Het = ho, Obs_Hom = 1 - ho, Exp_Het = he, Exp_Hom = 1 - he, Pi = pi,
    Fis = if (pi > 0) (pi - ho) / pi else 0, poly = as.numeric(he > 0),
    alt_here = sum(r$a[typed]) + sum(r$b[typed]), ref_here = 2 * N - sum(r$a[typed]) - sum(r$b[typed]))
}
fmt <- function(x) formatC(x, format = "f", digits = 5)
block_rows <- function(all_positions) {
  vapply(c("north", "south"), function(p) {
    members <- pop == p
    other <- pop != p
    per_site <- do.call(rbind, lapply(stats_rows, site_stats, members = members))
    other_site <- lapply(stats_rows, site_stats, members = other)
    ## Private alleles: an allele present here and absent from the other population.
    private <- sum(vapply(seq_along(stats_rows), function(j) {
      here <- site_stats(stats_rows[[j]], members)
      there <- other_site[[j]]
      if (is.null(here) || is.null(there)) return(0)
      (here[["alt_here"]] > 0 && there[["alt_here"]] == 0) +
        (here[["ref_here"]] > 0 && there[["ref_here"]] == 0)
    }, numeric(1)))
    columns <- c("N", "P", "Obs_Het", "Obs_Hom", "Exp_Het", "Exp_Hom", "Pi", "Fis")
    sites <- sum(typed_tags[, p]) * tag_len
    variant <- nrow(per_site)
    values <- if (all_positions) {
      ## Fixed sites add P = 1, Obs_Hom = Exp_Hom = 1, and 0 for the rest.
      fixed <- c(N = NA, P = 1, Obs_Het = 0, Obs_Hom = 1, Exp_Het = 0, Exp_Hom = 1, Pi = 0, Fis = 0)
      vapply(columns, function(cl) {
        if (cl == "N") return(c(mean(per_site[, "N"]), stats::var(per_site[, "N"])))
        x <- c(per_site[, cl], rep(fixed[[cl]], sites - variant))
        c(mean(x), stats::var(x))
      }, numeric(2))
    } else {
      vapply(columns, function(cl) c(mean(per_site[, cl]), stats::var(per_site[, cl])), numeric(2))
    }
    n_for_se <- if (all_positions) sites else variant
    triples <- as.vector(rbind(fmt(values[1, ]), fmt(values[2, ]),
                               fmt(sqrt(values[2, ]) / sqrt(n_for_se))))
    lead <- if (all_positions)
      c(p, private, sites, variant, sum(per_site[, "poly"]),
        fmt(100 * sum(per_site[, "poly"]) / variant))
    else c(p, private)
    paste(c(lead, triples), collapse = "\t")
  }, character(1))
}
triple_header <- paste(as.vector(rbind(c("Num_Indv", "P", "Obs_Het", "Obs_Hom", "Exp_Het",
                                         "Exp_Hom", "Pi", "Fis"), "Var", "StdErr")),
                       collapse = "\t")
writeLines(c("# Variant positions",
             paste0("# Pop ID\tPrivate\t", triple_header),
             block_rows(FALSE),
             "# All positions (variant and fixed)",
             paste0("# Pop ID\tPrivate\tSites\tVariant_Sites\tPolymorphic_Sites\t%Polymorphic_Loci\t",
                    triple_header),
             block_rows(TRUE)),
           file.path(out_dir, "example.sumstats_summary.tsv"))

cat(sprintf("Wrote %d SNP records, %d haplotype records and %d sites for %d individuals to %s\n",
            length(snp_rows), length(hap_rows), length(unlist(all_rows)), n, out_dir))
