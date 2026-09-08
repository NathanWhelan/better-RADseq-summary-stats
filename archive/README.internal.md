# Minimal RADseq diversity workflow

H<sub>o</sub>, H<sub>e</sub>, π, F<sub>IS</sub>, rarefied allelic richness and
rarefied private allelic richness from Stacks output, plus a correctly
calibrated test of whether two populations differ.

Read **`MINIMAL_WORKFLOW.md`** first. It is the documentation; the scripts print
numbers and point back to it rather than repeating it.

## Quick start

```bash
# one Stacks run
populations --in-path ./stacks_out --popmap popmap.tsv -O ./out \
            --min-gt-depth 10 -r 0.8 -p 2 --min-mac 3 \
            --max-obs-het 0.70 --fstats --vcf -t 8

# Fis, allelic richness, private allelic richness
Rscript diversity_stats.R out/populations.haps.vcf popmap.tsv --g=<gene_copies> --nboot=10000

# Ho, He, % polymorphic, and the per-sequenced-site (pi) values
Rscript diversity_stats.R out/populations.snps.vcf popmap.tsv --g=<gene_copies> --nboot=10000 --sites=<sites>

# does mean heterozygosity differ between populations?
Rscript het_between_pops.R out/populations.snps.vcf popmap.tsv 0.9 het_out
```

`--g` is required: the rarefaction size in GENE COPIES (10 diploids = 20),
no more than 2x your smallest population. `<sites>` is the `Sites` column of
the **All positions (variant and fixed)** block of
`populations.sumstats_summary.tsv`. `popmap.tsv` is two tab-separated
columns, no header: sample ID, population.

Each run prints a `TAKE FROM THIS RUN` box naming the statistics that file is
the right source for. Do not mix them up — see "Which file for which statistic".

## Files

| file | what it is |
|---|---|
| `MINIMAL_WORKFLOW.md` | the documentation: what to run, why, and how to write it up |
| `diversity_stats.R` | diversity statistics with block-bootstrap intervals |
| `het_between_pops.R` | between-population test, individual as the replicate |
| `diversity_core.R` | estimators, rarefaction, self-tests (`--selftest`) |
| `haps_common.R` | VCF and popmap readers |
| `test/` | test infrastructure: `run_tests.sh` (full suite), fixtures, golden values |
| `archive/` | superseded/abandoned material, kept for history (see `archive/README.md`) |

## Requirements

R. Base R is sufficient. `hierfstat` is used as the engine when installed and
the scripts fall back to validated internal code otherwise; each run says which.

## Testing

```bash
bash test/run_tests.sh                   # everything
Rscript diversity_core.R --selftest      # estimators, rarefaction, Stacks formulas
Rscript het_between_pops.R --selftest    # type I error of each candidate test
```

The self-tests check estimators against brute-force Monte Carlo, reproduce the
type I error table that motivates the individual-level test, and verify the
Stacks formulas against real published Stacks output.

## What this deliberately does not do

Null-allele screening, paralog screening, N<sub>e</sub> estimation, HWE
filtering, population structure. See "Deliberately not here" in the workflow for
why, and what to reach for instead.
