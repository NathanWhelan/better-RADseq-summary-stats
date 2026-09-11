# RADdiversity

Population-genetic statistics from RAD-seq data, built around Stacks 2
`populations` output: heterozygosity, gene and nucleotide diversity,
F<sub>IS</sub>, allelic richness, differentiation, individual inbreeding,
and tests of whether populations differ -- with the choices that most often
go wrong made carefully, and explained.

## Why use it

Three things routine RAD-seq summaries tend to get wrong:

1. **F<sub>IS</sub>.** Stacks averages per-locus F<sub>IS</sub> ratios and
   divides by an expected heterozygosity that assumes F<sub>IS</sub> = 0.
   Here, H<sub>e</sub> is Nei & Chesser's (1983) estimator, unbiased at any
   F<sub>IS</sub>, and F<sub>IS</sub> = 1 − ΣH<sub>o</sub>/ΣH<sub>e</sub>
   over loci.
2. **Comparing populations.** Bootstrapping or testing over loci treats loci
   as independent replicates of a population mean. They are repeated
   measurements on the same animals: when individuals differ in inbreeding,
   such tests reject a true null 55–76% of the time at a nominal 5%
   (`het_between_pops_selftest()`). Here, populations are compared on one
   number per individual -- heterozygosity, or individual F.
3. **Uncertainty.** SNPs on one RAD tag are linked, so standard errors and
   intervals resample whole RAD loci, never SNP rows. And when individuals
   differ in inbreeding, even that is too optimistic for population-level
   inference: the package measures this (g2) and also gives standard errors
   over individuals.

## Install

```r
# install.packages("remotes")
remotes::install_github("NathanWhelan/better-RADseq-summary-stats")
```

R ≥ 3.5; base R is enough. Optional: `hierfstat` (Weir & Goudet's beta in
`differentiation_stats()`, and `hierfstat_check = TRUE` cross-checks).

## Quick start

On the toy dataset that ships with the package (80 RAD loci; populations of
4 and 3 individuals -- far too few for real inference):

```r
library(RADdiversity)
haps   <- system.file("extdata", "small.haps.vcf",   package = "RADdiversity")
snps   <- system.file("extdata", "small.snps.vcf",   package = "RADdiversity")
popmap <- system.file("extdata", "small_popmap.tsv", package = "RADdiversity")
sumstats <- system.file("extdata", "populations.sumstats_summary.tsv", package = "RADdiversity")

div_haps <- diversity_stats(haps, popmap, g = 4)                     # FIS, Ar, privAr
div_snps <- diversity_stats(snps, popmap, g = 4, sites = sumstats)   # Ho, He, per-site values
div_haps                                    # printing shows the full report
het <- het_between_pops(haps, popmap, min_call = 0.5)                # do populations differ?
het$pairwise_F_tests
```

Results are returned as tables (`div_haps$per_population`,
`div_haps$richness`, ...); add `outdir = "results"` to also write them as TSV
files. On your own data, `g` is the rarefaction size in **gene copies**
(10 diploids = 20), at most twice your smallest population.

From the shell, the same analyses run with the scripts installed alongside
the package:

```bash
RD=$(Rscript -e 'cat(system.file("scripts", package = "RADdiversity"))')
Rscript $RD/diversity_stats.R out/populations.haps.vcf popmap.tsv --g=20
Rscript $RD/diversity_stats.R out/populations.snps.vcf popmap.tsv --g=20 \
        --sites=out/populations.sumstats_summary.tsv
Rscript $RD/het_between_pops.R out/populations.snps.vcf popmap.tsv --min-call=0.9 --outdir=het_out
```

Run any script with no arguments for its options. `vignette("workflow",
package = "RADdiversity")` walks through a whole analysis.

## Before you start: two `populations` flags

```bash
populations --in-path ./stacks_out --popmap popmap.tsv -O ./out \
            --min-gt-depth 10 -r 0.8 -p 2 --min-mac 3 \
            --max-obs-het 0.70 --fstats --vcf -t 8
```

* `--min-gt-depth 10` (Stacks ≥ 2.67): a heterozygote seen in only a few
  reads is easily called a confident homozygote -- 25% of the time at 3
  reads, 0.2% at 10. The flag blanks such genotypes instead.
* `-r 0.8 -p 2`, **not** `-R 0.8`: a missing-data filter applied within each
  population. A pooled filter lets the larger population's coverage carry the
  smaller one, which then carries more missing data -- an artificial
  difference between exactly the two samples you mean to compare.
* No `--hwe`: a heterozygote deficit is what F<sub>IS</sub> measures.

Add `--vcf-all` for an all-sites VCF (Stacks ≥ 2.62), which `pi_allsites()`
uses to compute nucleotide diversity directly.

## What each function answers

| question | function |
|---|---|
| Ho, He, F<sub>IS</sub>, % polymorphic, allelic and private allelic richness | `diversity_stats()` |
| nucleotide diversity π, d<sub>xy</sub> and per-individual heterozygosity per sequenced site | `pi_allsites()` (all-sites VCF), or `diversity_stats(sites = ...)` (approximation) |
| do populations differ in diversity? in inbreeding? | `het_between_pops()` |
| each individual's inbreeding coefficient | `individual_inbreeding()` |
| do individuals differ in inbreeding (g2)? | `identity_disequilibrium()` |
| F<sub>ST</sub>, Jost's D, Weir & Goudet's beta | `differentiation_stats()` |
| are there close relatives in the sample? | `kinship_check()` |
| Hardy–Weinberg departures (a report, never a filter) | `hwe_test()` |
| read and filter | `read_stacks_vcf()`, `read_popmap()`, `filter_call_rate()`, `filter_maf()`, `filter_mac()`, `filter_max_het()`, `filter_thin_one_snp()`, `filter_low_conf_alt()`, `locus_allele_stats()`, `low_conf_alt_sensitivity()` |
| export | `write_vcf()`, `write_plink()`, `write_structure()`, `write_genepop()`, `write_fstat()`, `write_radpainter()` |
| read Stacks' `sumstats_summary.tsv` | `read_sumstats_summary()` |
| the individual estimators | `?estimators` |
| check the installation | `diversity_core_selftest()`, `het_between_pops_selftest()` |

## Which file for which statistic

Stacks writes one record per SNP (`populations.snps.vcf`) and one per RAD
locus (`populations.haps.vcf`). They answer different questions:

| statistic | file | why |
|---|---|---|
| H<sub>o</sub>, H<sub>e</sub>, per-site values (π) | `.snps` | per-site values are on a scale other studies can compare |
| F<sub>IS</sub> | `.haps` | a ratio, so the scale cancels; multi-allelic loci are more precise |
| allelic richness, private allelic richness | `.haps` | on biallelic SNPs richness can only be 1 or 2 |
| between-population tests | either -- say which | |

H<sub>o</sub> and H<sub>e</sub> are **not** comparable between the two files
(haplotype H<sub>o</sub> asks "is this tag heterozygous?", per-site
H<sub>o</sub> "is this site?"); F<sub>IS</sub> is. Each printed report says
which of its numbers to take.

## Which standard error to report

* `_se`, `_lo`, `_hi` resample **RAD loci**. They hold the sampled individuals
  fixed: "what if I had typed different loci in these same animals?".
* `_se_ind` (`diversity_stats(se_individuals = TRUE)`) is a
  delete-one-**individual** jackknife: the uncertainty from which individuals
  were sampled.

With many loci and exchangeable individuals the locus-based intervals are
enough. When individuals differ in inbreeding they are not
(`inst/sims/coverage_se.R`; 2 × 15 individuals, 1,000 loci, mean F = 0.1):

| SD of F among individuals | coverage of 95% interval, H<sub>o</sub>: locus SE / individual SE | F<sub>IS</sub>: locus SE / individual SE |
|---|---|---|
| 0    | 98.8% / 94.2% | 95.5% / 94.2% |
| 0.10 | 65.5% / 94.5% | 54.2% / 93.0% |
| 0.25 | 37.0% / 93.0% | 27.5% / 93.5% |

So: check g2 (`identity_disequilibrium()`, also printed by
`het_between_pops()`). If its interval excludes 0, report `_se_ind`. Either
way, never compare populations by overlapping intervals -- use
`het_between_pops()`.

## Glossary

* **H<sub>o</sub>** -- observed heterozygosity: the fraction of typed
  individuals that are heterozygous, averaged over loci.
* **H<sub>e</sub> (H<sub>s</sub>, gene diversity)** -- the chance that two
  gene copies drawn from the population differ (Nei & Chesser 1983).
* **π (nucleotide diversity)** -- the same quantity averaged over every
  sequenced site, invariant ones included.
* **F<sub>IS</sub>** -- 1 − H<sub>o</sub>/H<sub>e</sub> for a population,
  summed over loci; **F** -- the same for one individual.
* **gene copies** -- a diploid carries two: 10 individuals are 20 gene copies.
* **RAD locus vs record** -- a locus is one RAD tag; a SNP VCF has one record
  per SNP, so several records can share a locus.
* **A<sub>r</sub>, privA<sub>r</sub>** -- allelic and private allelic
  richness rarefied to the same number of gene copies (Kalinowski 2004).
* **g2** -- identity disequilibrium: how much more often an individual is
  heterozygous at two loci at once than two different individuals are; 0 when
  individuals do not differ in inbreeding (David et al. 2007).
* **jackknife / bootstrap** -- recompute a statistic with one unit left out /
  with units resampled; here the unit is a RAD locus (or, for `_se_ind`, an
  individual).

## What this package does not do

| left out | why |
|---|---|
| HWE filtering | a heterozygote deficit is the signal. `hwe_test()` reports departures; nothing removes loci on them. |
| null-allele correction | restriction-site null alleles cannot be removed by depth filtering; `vignette("reviewer-faq")` has the sentence to write instead. |
| paralog detection beyond an excess-heterozygosity screen | `filter_max_het()` is that screen; dedicated tools (e.g. HDplot) go further. |
| N<sub>e</sub>, AMOVA, neutrality tests | different questions. |

## How it relates to other tools

* **Stacks** reports H<sub>o</sub>, π and F<sub>IS</sub> per site and per
  population. This package recomputes F<sub>IS</sub> (ratio of sums,
  Nei–Chesser H<sub>e</sub>), adds rarefied and private allelic richness,
  standard errors over RAD loci and individuals, and individual-level tests.
  `vignette("rationale")` maps every Stacks column onto its counterpart here.
* **hierfstat** (Goudet 2005) uses the same H<sub>s</sub> and Weir &
  Cockerham estimators; the package's own implementations are checked
  against it in the tests. Its bootstrap resamples loci as independent rows,
  ignoring linkage within RAD tags.
* **dartR** can correct heterozygosity for invariant sites (Schmidt et al.
  2021); **snpR** (Hemstrom & Jones 2023) computes many SNP statistics across
  categorical metadata.
* **pixy** (Korunes & Samuk 2021) computes π and d<sub>xy</sub> from
  all-sites VCFs; `pi_allsites()` uses the same estimator, with standard
  errors over RAD loci.
* **inbreedR** (Stoffel et al. 2016) computes g2; `identity_disequilibrium()`
  estimates it from its definition, using exactly the locus pairs typed when
  data are missing.

## Citation

`citation("RADdiversity")`. Please also cite the methods you report -- the
references are in each function's help page and in `vignette("rationale")`.

## Further reading

* `vignette("workflow", package = "RADdiversity")` -- one analysis from start
  to finish.
* `vignette("rationale", package = "RADdiversity")` -- why these flags,
  estimators and tests; formulas; the Stacks column mapping.
* `vignette("reviewer-faq", package = "RADdiversity")` -- methods text to
  adapt, and answers to the usual reviewer objections.

## Key references

David, P., Pujol, B., Viard, F., Castella, V. & Goudet, J. (2007) Reliable
selfing rate estimates from imperfect population genetic data. *Molecular
Ecology* 16:2474–2487.

Goudet, J. (2005) HIERFSTAT, a package for R to compute and test hierarchical
F-statistics. *Molecular Ecology Notes* 5:184–186.

Hemstrom, W. & Jones, M. (2023) snpR: user friendly population genomics for
SNP data sets with categorical metadata. *Molecular Ecology Resources*
23:962–973.

Kalinowski, S.T. (2004) Counting alleles with rarefaction: private alleles
and hierarchical sampling designs. *Conservation Genetics* 5:539–543.

Korunes, K.L. & Samuk, K. (2021) pixy: Unbiased estimation of nucleotide
diversity and divergence in the presence of missing data. *Molecular Ecology
Resources* 21:1359–1368.

Nei, M. & Chesser, R.K. (1983) Estimation of fixation indices and gene
diversities. *Annals of Human Genetics* 47:253–259.

Rochette, N.C., Rivera-Colón, A.G. & Catchen, J.M. (2019) Stacks 2.
*Molecular Ecology* 28:4737–4754.

Schmidt, T.L., Jasper, M.-E., Weeks, A.R. & Hoffmann, A.A. (2021) Unbiased
population heterozygosity estimates from genome-wide sequence data. *Methods
in Ecology and Evolution* 12:1888–1898.

Stoffel, M.A., Esser, M., Kardos, M., Humble, E., Nichols, H., David, P. &
Hoffman, J.I. (2016) inbreedR: an R package for the analysis of inbreeding
based on genetic markers. *Methods in Ecology and Evolution* 7:1331–1339.

Van Dongen, S. (1995) How should we bootstrap allozyme data? *Heredity*
74:445–447.

Weir, B.S. & Cockerham, C.C. (1984) Estimating F-statistics for the analysis
of population structure. *Evolution* 38:1358–1370.
