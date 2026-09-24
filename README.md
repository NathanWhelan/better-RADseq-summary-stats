# RADdiversity

> **Please note:** This was mostly written by Claude AI with considerable
> human input and checks. Both Sonnet 5 and Opus 5 models were used. The
> README and Vignettes have been human checked and edited. If something seems
> off, please post an issue. If you notice something that's wrong, or even
> seems potentially wrong, I would love to hear feedback. --N. Whelan (Sept 8,
> 2026).

RADdiversity calculates population-genetic statistics from RAD-seq data. It
was built for output from Stacks 2, but VCFs from other assembly pipelines work
too. Genotypes are assumed to be diploid.

## Why use this package

Routine RAD-seq summaries often get three things wrong:

1. **F<sub>IS</sub>.** Stacks averages per-site F<sub>IS</sub> ratios, and its
   expected heterozygosity assumes F<sub>IS</sub> = 0. Here, H<sub>e</sub> is
   Nei & Chesser's (1983) estimator, which is unbiased at any F<sub>IS</sub>,
   and F<sub>IS</sub> = 1 − ΣH<sub>o</sub>/ΣH<sub>e</sub> over loci (a ratio of
   sums).
2. **One package to generate useful statistics with clear methods.** Population genetics has a history of the same, or similar, terms being conflated (e.g., F<sub>IS</sub> is calculated differently by STACKS than other approaches. Sometimes, it can be unclear what method or formula is being used to calculate a metric. This package prioritizes calculating the most robust metrics for making conclusions from RADseq data. This package includes methods that have previously not previously been available in a single package.

4. **Uncertainty.** SNPs on one RAD tag are linked, so every standard error
   resamples whole RAD loci, never single SNPs. The standard errors the
   package recommends also count which individuals were sampled.

## Install

```r
install.packages("remotes")
remotes::install_github("NathanWhelan/better-RADseq-summary-stats", build_vignettes = TRUE)
```

`build_vignettes = TRUE` installs the two guides, `vignette("workflow")` and
`vignette("rationale")`. It needs the knitr and rmarkdown packages and pandoc
(RStudio includes pandoc). R ≥ 4.0. Optional:
hierfstat (Weir & Goudet's beta) and adegenet (`as_genind()`, `as_genlight()`).

## Quick start

The package ships a simulated example: two populations of 15 and 10
individuals, 1,000 RAD loci, 741 of them variable.

```r
library(RADdiversity)
haps     <- system.file("extdata", "example.haps.vcf.gz", package = "RADdiversity")
snps     <- system.file("extdata", "example.snps.vcf.gz", package = "RADdiversity")
popmap   <- system.file("extdata", "example_popmap.tsv", package = "RADdiversity")
sumstats <- system.file("extdata", "example.sumstats_summary.tsv", package = "RADdiversity")

set.seed(2024)                                                       # reproducible intervals
div_haps <- diversity_stats(haps, popmap, g = 16)                    # Fis, Ar, rarefied private alleles
div_snps <- diversity_stats(snps, popmap, g = 16, sites = sumstats)  # Ho, He, per-site values
summary(div_haps)                     # each estimate (SE), what to take, checks
diversity_table(div_snps, div_haps)   # one table for your paper, with a caption

het <- het_between_pops(snps, popmap)       # do populations differ in heterozygosity?
summary(het)
dif <- differentiation_stats(haps, popmap)  # FST and Jost's D
summary(dif)
```

On a large dataset the slow parts are the bootstrap (`nboot`, default 10,000)
and the `beta` step of `differentiation_stats()`, which calls hierfstat. In testing, this was not a huge issue. If you are concerned about slow analyses, for a
first look, use `nboot = 1000` and `beta = FALSE`. Run again with the defaults
for the results you report. `verbose = FALSE` silences progress messages.

### Using your own files

Replace the four `system.file()` lines with paths to your files.

| variable | your file | what it is |
|---|---|---|
| `snps` | `populations.snps.vcf` | Stacks output: one record per SNP. `.vcf.gz` works too. |
| `haps` | `populations.haps.vcf` | Stacks output: one record per RAD locus, with haplotype alleles. |
| `popmap` | `popmap.tsv` | The popmap you gave Stacks. No header; each line is the sample name, a TAB, then the population. Names must match the VCF. |
| `sumstats` | `populations.sumstats_summary.tsv` | Stacks output. Optional: only for per-site values (`sites =`). |

**Choosing `g`.** `g` is the rarefaction size in **gene copies**: 10 diploids
are 20. Use g = 2 × r × n, rounded down, where n is the number of individuals
in your smallest population and r is the `-r` value you gave Stacks. For
example, `-r 0.8` and a smallest population of 10 give g = 16. Then check
`div_haps$richness$privAr_n`, the number of loci the rarefied private alleles
rest on.
If it is far below the number of loci, try a smaller `g` (workflow vignette,
section 5).

Every analysis function takes a file path, or the object from
`read_stacks_vcf()`. Read the file yourself to look at it or filter it first:

```r
H    <- read_stacks_vcf(snps)            # a VCF -> R object
pops <- read_popmap(popmap, H$samples)   # checked against the VCF's sample names
H    <- filter_call_rate(H, min_call = 0.8, popmap = pops)
div_snps <- diversity_stats(H, pops, g = 16, sites = sumstats)
```

A VCF from another pipeline (ipyrad, dDocent) can be used as `snps`; skip
`haps`. If its RAD loci are not recognized, read it with
`read_stacks_vcf(snps, locus_from = "CHROM")` (see `?read_stacks_vcf`).

## Reading and saving results

* Printing a result shows its main tables, rounded.
* `summary(x)` is short: each estimate is printed as `estimate (SE)` with the one
  standard error the package recommends, what to take from this run, and
  checks marked `ok`, `info` or `look`.
* `summary(x, details = TRUE)` is the full report, with every uncertainty
  column and a note on how to read each one.
* Each table can be taken out by name at full precision:
  `div_haps$per_population`, `div_haps$richness`, `dif$pairwise`, ...
* `outdir = "results"` writes the tables as TSV files. The file names end in
  the VCF's type (`haps` or `snps`), so both runs can share one folder. For an
  object from `read_stacks_vcf()`, give `stem = "snps"` yourself.

Two names that may not be obvious:

* **Ar**, allelic richness: the number of different alleles a population
  shows at a locus in `g` gene copies, averaged over loci. Every population is
  cut to the same `g`, so a bigger sample does not look richer.
* **rarefied private alleles**: the expected number of alleles found only in
  this population when every population is cut to `g` gene copies at each
  locus, added up over the loci typed at `g` copies in every population
  (`priv_total` in the results; `privAr` is the same per locus). It is an
  expected value, so it need not be a whole number. It grows with the number
  of loci, so compare it only between populations analyzed together.

**From the shell.** Scripts installed with the package run the same analyses:

```bash
RD=$(Rscript -e 'cat(system.file("scripts", package = "RADdiversity"))')
Rscript $RD/diversity_stats.R out/populations.haps.vcf popmap.tsv --g=16 --outdir=results
Rscript $RD/het_between_pops.R out/populations.snps.vcf popmap.tsv --outdir=results
```

Run a script with no arguments to see its options.

## A table for your manuscript

```r
table1 <- diversity_table(div_snps, div_haps)
table1                                            # the table, then a caption
write.csv(table1, "table1_diversity.csv", row.names = FALSE)
table2 <- summary(dif)$tables$results             # FST (SE), D (SE) and beta for each pair
```

`table1` has one row per population: `n`, `Ho (SE)`, `He (SE)`, the per-site
values (if you gave `sites`), `Fis (SE)`, `Ar (SE)` and
`rarefied private alleles (SE)`.
It takes each column from the right VCF, and shows the one standard error the
package's simulations found most accurate. Printing it also gives a caption:
which VCF each column came from, the value of `g`, and which standard error
each column is. Add the filters and thresholds you used.

Everything else (the other SEs and intervals, `pct_poly`, counts such as
`privAr_n`) stays in the results and the `outdir` files, for a supplement.
`summary(het)$tables` and `summary(pi_res)$tables` hold the tables for the
other analyses. The workflow vignette ("A table for your manuscript") says
more.

## Isolation by distance

Do populations farther apart differ more? `isolation_by_distance()` takes the
pairwise F<sub>ST</sub> (or Jost's D) from `differentiation_stats()` and a
file of distances. The example has six sites along a river:

```r
ibd_vcf    <- system.file("extdata", "ibd_example.snps.vcf.gz", package = "RADdiversity")
ibd_popmap <- system.file("extdata", "ibd_example_popmap.tsv", package = "RADdiversity")
river_km   <- system.file("extdata", "ibd_example_distances.csv", package = "RADdiversity")

dif_river <- differentiation_stats(ibd_vcf, ibd_popmap)
read_distances(river_km)                      # check the file was read as you meant
ibd <- isolation_by_distance(dif_river, river_km, habitat = "1D")
summary(ibd)                                  # Mantel r, p, slope, what to report
plot(ibd, xlab = "River km")
isolation_by_distance(dif_river, river_km, habitat = "1D", stat = "D")   # Jost's D
isolation_by_distance(dif_river, river_km, habitat = "1D", exclude = "site6")
```

`habitat` is required: `"1D"` for populations along a line (a river, a
coastline), `"2D"` for populations spread over an area (Rousset 1997). The
Mantel test is one-sided and, with 7 or fewer populations, exact. The slope is
a point estimate (no SE). D is used as it is, not as D/(1 − D);
`vignette("rationale")` explains why.

**The distance file** is a CSV (or tab-separated) file in either layout below.
Names must match the popmap; their order does not matter. Any unit works (the
slope is per that unit). Use the distance the organisms travel: along the
river for stream animals, not a straight line.

A square matrix (one half may be left blank):

```
,site1,site2,site3
site1,0,15,27
site2,,0,12
site3,,,0
```

Or one row per pair:

```
pop1,pop2,distance
site1,site2,15
site1,site3,27
site2,site3,12
```

Extra populations in the file are ignored. A missing one stops the test: add
it, or leave it out with `exclude =`.

## Before you start: Stacks `populations` flags

```bash
populations --in-path ./stacks_out --popmap popmap.tsv -O ./out \
            -p <number of populations> -r 0.8 --min-gt-depth 6 \
            --max-obs-het 0.70 --vcf --vcf-all -t 8
```

* **Use `-r` with `-p` set to your number of populations, not a global `-R`.**
  With `-R`, a large, well-typed population can carry a small one, which then
  has more missing data. That creates a difference between exactly the samples
  you mean to compare.
* **`--min-gt-depth 6`** (Stacks ≥ 2.67) sets genotypes backed by few reads to
  missing. A heterozygote seen in few reads can look homozygous. In the
  author's tests this made little difference, so it is a precaution, not a
  requirement.
* **`--min-mac` or `--min-maf`** are common, but they remove rare variants,
  which are real diversity: about 21% of π at 10 diploids and 10% at 20. If you
  use one, report the threshold, and compare only with data filtered the same
  way.
* **`--vcf-all`** (Stacks ≥ 2.62) writes an all-sites VCF for `pi_allsites()`.
* **Per-site values** use `Sites` from `populations.sumstats_summary.tsv`. That
  count fits only data not filtered after `populations` (see
  `?diversity_stats`).

Data filtered in other ways work too. `vignette("workflow")` has a table of
what common filters change and what to report.

## What each function answers

| question | function |
|---|---|
| Ho, He, F<sub>IS</sub>, allelic richness, rarefied private alleles | `diversity_stats()` |
| one table for a manuscript | `diversity_table()` |
| nucleotide diversity (π), d<sub>xy</sub>, heterozygosity per sequenced site | `pi_allsites()` (all-sites VCF), or `diversity_stats(sites = )` |
| do populations differ in heterozygosity? in inbreeding? | `het_between_pops()` |
| each individual's inbreeding coefficient | `individual_inbreeding()` |
| do individuals differ in inbreeding (g2)? | `identity_disequilibrium()` |
| F<sub>ST</sub>, Jost's D, Weir & Goudet's beta | `differentiation_stats()` |
| do populations farther apart differ more (isolation by distance)? | `isolation_by_distance()`, `read_distances()` |
| are there close relatives in the sample? | `kinship_check()` |
| Hardy–Weinberg departures (a report, never a filter) | `hwe_test()` |
| read and filter | `read_stacks_vcf()`, `read_popmap()`, `filter_samples()`, `filter_call_rate()`, `filter_genotype_depth()`, `filter_maf()`, `filter_mac()`, `filter_max_het()`, `filter_thin_one_snp()`, `filter_low_conf_alt()`, `locus_allele_stats()`, `low_conf_alt_calls()`, `low_conf_alt_sensitivity()` |
| export | `write_vcf()`, `write_plink()`, `write_structure()`, `write_genepop()`, `write_fstat()`, `write_radpainter()`, `as_genind()` and `as_genlight()` (adegenet, dartR) |
| read Stacks' `sumstats_summary.tsv` | `read_sumstats_summary()` |
| the individual estimators | `?estimators` |
| check the installation | `diversity_core_selftest()`, `het_between_pops_selftest()` |

## Which file for which statistic

With `--vcf`, Stacks' `populations` writes two VCFs: one record per SNP
(`populations.snps.vcf`) and one record per RAD locus
(`populations.haps.vcf`). With `--write-single-snp` it writes only the SNP VCF.

| statistic | file | why |
|---|---|---|
| H<sub>o</sub>, H<sub>e</sub>, per-site values | `.snps` | a per-SNP value means the same in every study |
| F<sub>IS</sub> | `.haps` | a ratio, so the scale cancels; multi-allelic loci are more precise |
| allelic richness, rarefied private alleles | `.haps` | on a SNP, richness can only be 1 or 2 |
| tests between populations | either; say which | |

Haplotype allelic richness and rarefied private alleles depend on how many SNPs
each locus has. So compare them only among populations analyzed together (same
run, filters and `g`). Which population is richer is reliable; how much richer is
not. H<sub>o</sub> and H<sub>e</sub> cannot be compared between the two files
(a haplotype is heterozygous if *any* of its SNPs is); F<sub>IS</sub> can.
Each printed result says which numbers to take from it.

## Which standard error to report

A population's value is uncertain for two reasons: you typed some of the
genome's RAD loci, and you caught some of the population's individuals.

* `_se`, `_lo`, `_hi` resample **RAD loci** and hold the individuals fixed.
* `_se_ind` is a jackknife over **individuals** that holds the loci fixed.
* `_se_combined` puts both together: `sqrt(_se^2 + _se_ind^2)`.

**Report `_se_combined` for H<sub>o</sub>, H<sub>e</sub>, F<sub>IS</sub> and
the per-site values, and `_se` for A<sub>r</sub> and rarefied private alleles.** The
short `summary()` and `diversity_table()` show exactly these. In simulations
with a known truth (2 × 15 individuals, 1,000 loci, new loci and individuals
in every repeat), 95% intervals held the true value this often (95% is
perfect):

| | locus SE | individual SE | combined SE |
|---|---|---|---|
| H<sub>o</sub> | 41–96% | 88–91% | **93–98%** |
| H<sub>e</sub> | 95–96% | 72–77% | **98–99%** |
| F<sub>IS</sub> | 31–94% | 92–93% | **93–98%** |
| A<sub>r</sub> | **93–94%** | 98–100% | 100% (too wide) |
| private alleles, per locus | **92–93%** | 31–36% | 92–94% |

The low locus-SE values are for populations whose individuals differ in
inbreeding. The rarefied private alleles are the per-locus value times the
number of loci it is added up over, so their SE is the per-locus SE times that
number and holds the truth as often. No simulation has checked the
standard errors of F<sub>ST</sub>, D, beta, π, d<sub>xy</sub> or
d<sub>a</sub>; they count RAD loci only, and the summaries say so. Never
compare populations by overlapping intervals: use `het_between_pops()`.
`vignette("rationale")`, section 4, has the details.

## Glossary

* **H<sub>o</sub>**: observed heterozygosity, the share of typed individuals
  that are heterozygous, averaged over loci.
* **H<sub>e</sub>** (H<sub>s</sub>, gene diversity): the chance that two gene
  copies drawn from the population differ (Nei & Chesser 1983).
* **π** (nucleotide diversity): H<sub>e</sub> averaged over every sequenced
  site, invariant sites included, so it is much smaller than H<sub>e</sub>.
  `pi_allsites()` gives `pi_nc`, which stays unbiased when individuals are
  inbred (report it), and `pi`, as Stacks, pixy and VCFtools compute it (use it
  to compare with them).
* **F<sub>IS</sub>**: 1 − H<sub>o</sub>/H<sub>e</sub> for a population, summed
  over loci. **F** is the same for one individual.
* **gene copies**: a diploid carries two, so 10 individuals are 20 gene copies.
* **RAD locus vs record**: a locus is one RAD tag. A SNP VCF has one record per
  SNP, so several records can sit on one locus.
* **A<sub>r</sub>, rarefied private alleles**: see "Reading and saving results" above
  (Kalinowski 2004).
* **isolation by distance**: populations farther apart differ more, because
  individuals mostly mate near where they were born (Wright 1943).
* **g2**: identity disequilibrium, how much more often an individual is
  heterozygous at two loci at once than two different individuals are. It is 0
  when individuals do not differ in inbreeding (David et al. 2007).
* **jackknife / bootstrap**: recompute a statistic with one unit left out / with
  units drawn again at random. Here the unit is a RAD locus (or, for
  `_se_ind`, an individual).
* **standard error (SE)**: how much a number would change if the study were
  repeated. The estimate ± 1.96 SE is a 95% confidence interval.

## What this package does not do

| left out | why |
|---|---|
| HWE filtering | a heterozygote deficit is the signal. `hwe_test()` reports departures; its results can be used to make a list of loci to exclude in another program (e.g., *populations*). |
| null-allele correction | restriction-site null alleles cannot be removed by depth filtering; `diversity_stats()` reports `fis_by_call_rate` to detect them. |
| paralog detection beyond an excess-heterozygosity screen | `filter_max_het()` is that screen. Dedicated tools (e.g., HDplot; McKinney et al. 2017) go further. |
| partial Mantel tests | not valid for spatially structured data (Guillot & Rousset 2013). |
| N<sub>e</sub>, AMOVA, neutrality tests | different questions. |

## How it relates to other tools

* **Stacks** reports H<sub>o</sub>, π and F<sub>IS</sub> per site and per
  population. This package recomputes F<sub>IS</sub> (ratio of sums,
  Nei–Chesser H<sub>e</sub>), and adds rarefied allelic richness and rarefied
  private alleles, standard errors over RAD loci and individuals, and individual-level
  tests. `vignette("rationale")` maps every Stacks column onto its counterpart
  here.
* **hierfstat** (Goudet 2005) uses the same H<sub>s</sub> and Weir & Cockerham
  estimators; the tests check this package against it. One deliberate
  difference: F<sub>ST</sub> skips records where only one of the compared
  populations is typed, as VCFtools and Stacks do; `hierfstat::wc()` keeps
  them, which pulls F<sub>ST</sub> toward 0.
* **pixy** (Korunes & Samuk 2021) computes π and d<sub>xy</sub> from all-sites
  VCFs. `pi_allsites()` uses the same estimator, adds standard errors over RAD
  loci, and adds `pi_nc`.
* **adegenet** and **dartR**: `as_genind()` and `as_genlight()` hand them the
  (filtered) data. **snpR** (Hemstrom & Jones 2023) computes many SNP
  statistics across categorical metadata.
* **inbreedR** (Stoffel et al. 2016) computes g2. `identity_disequilibrium()`
  estimates it from its definition, uses exactly the locus pairs typed, and
  leaves out pairs of SNPs on the same RAD locus.

## Citation and further reading

`citation("RADdiversity")`. Please also cite the methods you report; the
references are on each function's help page and in `vignette("rationale")`.

Questions, problems and suggestions: open an issue on GitHub. The
repository's `.github/CONTRIBUTING.md` says what to include, and how to run
the tests and simulations.

* `vignette("workflow", package = "RADdiversity")`: one analysis from start to
  finish.
* `vignette("rationale", package = "RADdiversity")`: why these filters,
  estimators and tests; formulas; the Stacks column mapping.

## Key references

David, P., Pujol, B., Viard, F., Castella, V. & Goudet, J. (2007) Reliable
selfing rate estimates from imperfect population genetic data. *Molecular
Ecology* 16:2474–2487.

Goudet, J. (2005) HIERFSTAT, a package for R to compute and test hierarchical
F-statistics. *Molecular Ecology Notes* 5:184–186.

Guillot, G. & Rousset, F. (2013) Dismantling the Mantel tests. *Methods in
Ecology and Evolution* 4:336–344.

Hemstrom, W. & Jones, M. (2023) snpR: user friendly population genomics for
SNP data sets with categorical metadata. *Molecular Ecology Resources*
23:962–973.

Jombart, T. (2008) adegenet: a R package for the multivariate analysis of
genetic markers. *Bioinformatics* 24:1403–1405.

Kalinowski, S.T. (2004) Counting alleles with rarefaction: private alleles
and hierarchical sampling designs. *Conservation Genetics* 5:539–543.

Korunes, K.L. & Samuk, K. (2021) pixy: Unbiased estimation of nucleotide
diversity and divergence in the presence of missing data. *Molecular Ecology
Resources* 21:1359–1368.

Mantel, N. (1967) The detection of disease clustering and a generalized
regression approach. *Cancer Research* 27:209–220.

McKinney, G.J., Waples, R.K., Seeb, L.W. & Seeb, J.E. (2017) Paralogs are
revealed by proportion of heterozygotes and deviations in read ratios in
genotyping-by-sequencing data from natural populations. *Molecular Ecology
Resources* 17:656–669.

Nei, M. & Chesser, R.K. (1983) Estimation of fixation indices and gene
diversities. *Annals of Human Genetics* 47:253–259.

Rochette, N.C., Rivera-Colón, A.G. & Catchen, J.M. (2019) Stacks 2.
*Molecular Ecology* 28:4737–4754.

Rousset, F. (1997) Genetic differentiation and estimation of gene flow from
F-statistics under isolation by distance. *Genetics* 145:1219–1228.

Schmidt, T.L., Jasper, M.-E., Weeks, A.R. & Hoffmann, A.A. (2021) Unbiased
population heterozygosity estimates from genome-wide sequence data. *Methods
in Ecology and Evolution* 12:1888–1898.

Stoffel, M.A., Esser, M., Kardos, M., Humble, E., Nichols, H., David, P. &
Hoffman, J.I. (2016) inbreedR: an R package for the analysis of inbreeding
based on genetic markers. *Methods in Ecology and Evolution* 7:1331–1339.

Weir, B.S. & Cockerham, C.C. (1984) Estimating F-statistics for the analysis
of population structure. *Evolution* 38:1358–1370.

Wright, S. (1943) Isolation by distance. *Genetics* 28:114–138.
