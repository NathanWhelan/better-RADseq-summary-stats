# Minimal workflow — H<sub>o</sub>, H<sub>e</sub>, F<sub>IS</sub>, and whether two populations differ

Four changes to an ordinary Stacks run. Together they took the reported
F<sub>IS</sub> values from 0.14 and 0.04 to about 0.07 and 0.05, and turned a
significant difference between populations into a non-significant one.

| # | change | where it lives |
|---|---|---|
| 1 | `--min-gt-depth 10` | a `populations` flag |
| 2 | `-r 0.8 -p 2` instead of `-R` | a `populations` flag |
| 3 | Nei & Chesser H<sub>e</sub>, F<sub>IS</sub> as a ratio of sums | `diversity_stats.R` |
| 4 | Welch's *t* on per-individual heterozygosity | `het_between_pops.R` |

Two flags and two scripts. Three commands total. Everything else in the larger
package is diagnosis, and none of it changes these numbers.


**Contents.** Step 0 (what you need) · Step 1 (`populations`, and the two flags
that matter) · Step 2 (the statistics: terminology, which file for which
statistic, the estimators, denominators) · Step 3 (the between-population test)
· Step 4 (what to write) · Formulas and the Stacks mapping · Reviewer
objections · What is deliberately not here.

---

## Step 0 — what you need

- Stacks 2, already run through `gstacks` (or `ustacks`/`cstacks`/`sstacks`).
- `popmap.tsv`: two columns, no header, `sample_id <TAB> population`.
- R. Base R is enough. `hierfstat` is used as the engine when installed
  (`basic.stats()` for per-locus H<sub>o</sub>/H<sub>s</sub>,
  `allelic.richness()`, `wc()` for Weir & Cockerham). Without it the scripts
  fall back to internal code checked against hierfstat's own source over 30
  random datasets: per-locus H<sub>s</sub> agreed to 5.0e-13, per-population
  F<sub>IS</sub> to 3.2e-13, rarefied allelic richness exactly. Each run states
  which engine it used. Rarefied *private* allelic richness, the block bootstrap
  over RAD loci, and the autosomal conversion are in no package and are
  implemented here, validated against brute-force Monte Carlo in
  `diversity_core.R --selftest`.
- These four scripts, **in one directory** (they `source()` each other):
  `haps_common.R`, `diversity_core.R`, `diversity_stats.R`, `het_between_pops.R`.

Check both scripts run before you point them at real data:

```bash
Rscript diversity_core.R    --selftest
Rscript het_between_pops.R  --selftest
```

The first checks the estimators against brute-force Monte Carlo. The second
reproduces the type I error table in Step 3 on simulated data, so you can see
the problem it fixes rather than take it on faith.

---

## Step 1 — run `populations` once

```bash
populations --in-path ./stacks_out --popmap popmap.tsv -O ./out \
            --min-gt-depth 10 \
            -r 0.8 -p 2 \
            --min-mac 3 \
            --max-obs-het 0.70 \
            --fstats --vcf -t 8
```

This writes `out/populations.snps.vcf` and `out/populations.haps.vcf`. Stacks 2
exports haplotypes by default (that is what `--no-hap-exports` turns *off*), so
one run gives you both files. Older builds spelled some options with
underscores (`--min_gt_depth`, `--max_obs_het`) and needed `-H` for the
haplotype VCF — run `populations -h` and match your own build rather than
copying this verbatim.

### The two flags that matter

**`--min-gt-depth 10`** is the single most important flag in the pipeline.

A genotype is called from the reads that happened to be sequenced. If a true
heterozygote is covered by only a few reads, there is a real chance that every
read came from the same chromosome, and the genotyper then calls a **confident
homozygote** — not a low-quality call, a confident wrong one. The probability is
2<sup>1−d</sup> at depth *d*:

| reads backing the genotype | chance a true heterozygote is called homozygous |
|---|---|
| 3 | 25% |
| 5 | 6% |
| 10 | 0.2% |
| 15 | 0.006% |

The flag turns any genotype below the threshold into a blank instead. It will
**raise** your missing-data rate, and that is the point: a visible blank is
better than an invisible mistake. This is the entire cure for stochastic
dropout — at 10 reads the problem is essentially gone.

**Version requirement.** `--min-gt-depth` was added to `populations` in Stacks
2.67 (July 2024). This is not a spelling difference like the underscore/hyphen
change noted below — earlier builds have no equivalent flag in `populations` at
all. On Stacks < 2.67, apply the depth filter downstream instead (e.g.
`vcftools --minDP 10 --recode` on `populations.snps.vcf`, then rebuild the
haplotype-level calls accordingly) rather than assuming a differently-named
flag will do it. Check `populations -h` for the flag before assuming it exists.

**`-r 0.8 -p 2`, and never `-R`.** These look interchangeable and are not.

- `-r 0.8 -p 2` — keep a locus if it is genotyped in ≥80% of individuals
  **within each population**, in at least 2 populations.
- `-R 0.8` (`--min-samples-overall`) — keep a locus if it is genotyped in ≥80%
  of individuals **pooled across all populations**.

With a pooled threshold, a locus can pass by being well genotyped in your larger
sample while being mostly blank in your smaller one. The bigger sample absorbs
the allowed gaps and the smaller one ends up carrying more missing data. That
asymmetry is exactly what makes one population's F<sub>IS</sub> look higher than
the other's for reasons that have nothing to do with biology, and it is why one
population in this project originally carried roughly twice the other's
missingness. Since a between-population comparison is the whole question, a
filter that treats the two populations differently is disqualifying.

### The other flags

`--min-mac 3` (minor allele count ≥3 — prefer a count to a frequency, because a
fixed MAF demands different numbers of copies from differently sized samples),
`--max-obs-het 0.70` (a crude paralog screen; users have reported Stacks
retaining sites above the stated threshold, so spot-check against your own
`sumstats.tsv`), `--fstats`, `--vcf`. These are conventional and are **not** part
of the four changes. Set them how you like; just set them the same way for both
populations.

Note there is no `--hwe` filtering here. A heterozygote deficit is exactly what
inbreeding produces, so an HWE filter deletes the signal you are trying to
measure.

### The filter you did not set, and why `-r` matters more than it looks

**DEFAULT as of 2026-09: available-data, not complete-case.** `diversity_stats.R`
now includes a locus in a population's statistics once that population has
`min_n` (default 2 — the floor below which H<sub>s</sub>/F<sub>IS</sub> are
undefined) typed individuals there. This is decided **per population,
independently, locus by locus** — a population missing a few individuals at a
locus does not cost any *other* population that locus. The aggregation this
package already used (F<sub>IS</sub> = 1 − Σ(H<sub>o</sub>)/Σ(H<sub>s</sub>), a
ratio of sums, `diversity_core.R`) combines loci with different per-locus
sample sizes correctly by construction — a locus with fewer typed individuals
simply contributes proportionally less to the sums — so relaxing the inclusion
rule required no change to the estimators themselves. `diversity_core.R
--selftest` includes a Monte Carlo check (#11) confirming this stays within 5%
of the true F<sub>IS</sub> under locus-driven (not just missing-at-random)
dropout, the pattern real RAD data actually shows.

**`--complete-case` reproduces the old rule**: a record is used only if EVERY
individual of EVERY population is genotyped there (Schmidt et al. 2021
recommendation (b); this was the only default before). That is a JOINT
condition across all populations on the same record, so one population's
missing individual costs every other population that record too. It interacts
hard with `-r 0.8`: if genotypes were missing independently at rate *m*, the
chance a record is complete across 25 individuals is (1−*m*)<sup>25</sup>:

| missing rate | records surviving |
|---|---|
| 2% | 60% |
| 5% | 28% |
| 10% | 7% |

On a test file with 10% missing genotypes, 154 of 2,316 records survived under
`--complete-case`. Real missingness in RAD data is typically clustered by
**locus**, not spread evenly across individuals or concentrated in a few bad
libraries — a locus's own depth/assembly behavior fails it broadly — which is
exactly the case `--complete-case` handles worst and the default handles well
(see the "available-data" self-test above, and the real-data comparison in
"Formulas, and what maps onto what in Stacks": `He`/`He_autosomal` moved from
roughly half of Stacks' native values under `--complete-case`'s ~6–9% retention
to within ~1% of them once available data was used).

Reasons you might still want `--complete-case`: reproducing a previously
published number from this package, or a reviewer specifically asking for the
Schmidt et al. (2021) rule. Otherwise, run Step 2 with the default and read the
per-population retention lines it prints (loci usable, mean/min typed-*n*)
rather than the old single retention percentage — coverage is now a
per-population question, not one dataset-wide number.

### Write down one number before you move on

Open `out/populations.sumstats_summary.tsv`. It has two blocks. Record the
**`Sites` column from the block headed "All positions (variant and fixed)"** —
the total number of nucleotides the run examined. You need it in Step 2.

---

## Step 2 — compute the statistics

```bash
# FIS, rarefied allelic richness, rarefied private allelic richness
Rscript diversity_stats.R out/populations.haps.vcf popmap.tsv --g=<gene_copies> --nboot=10000

# Ho, He, % polymorphic, and He per sequenced site -- SNP VCF, pass --sites=<sites>
Rscript diversity_stats.R out/populations.snps.vcf popmap.tsv --g=<gene_copies> --nboot=10000 --sites=<sites>
```

`<vcf>` and `<popmap>` are positional; everything else is a flag:
`--g=N` (required), `--nboot=N`, `--sites=N`, `--min-n=N`, `--complete-case`.
Each run prints a `TAKE FROM THIS RUN` box naming the statistics it is the
right file for, so you do not have to remember the split.

`--g` is the rarefaction size in **gene copies**, not individuals: 10 diploids
is g = 20. It is required and has no default -- pick it deliberately (no more
than twice your smallest population's size) rather than let it silently track
whichever population happens to be smallest.

**Do not pass `--sites` on the haplotype run.** See below.

### Two words that are not synonyms

Get these apart before anything else, because the rest of Step 2 depends on it
and mixing them produces a number that is wrong rather than merely mislabelled.

**Gene diversity** (H<sub>e</sub>, H<sub>s</sub>) = 1 − Σp². A property of **one
locus**: the chance that two gene copies drawn from the population are different
alleles.

**Nucleotide diversity** (π) = gene diversity computed **at a single site** and
averaged over **every sequenced site**, monomorphic ones included.

So π is not a different parameter. It is gene diversity under two extra
conditions — one record must be one *site*, and the denominator must be *all*
sites. Break either and what you have is gene diversity, not π:

- a haplotype VCF breaks the first: a record is a whole tag
- averaging over variant sites only breaks the second

**And the estimator is part of the label.** Published π values, Stacks' `Pi`
column and VCFtools all use the (2n/(2n−1)) form. This package's H<sub>e</sub> is
Nei–Chesser, which stays unbiased when F<sub>IS</sub> ≠ 0 where that form does
not. They estimate the same parameter and differ by F/(2n−1) — under 0.5% here —
but they are not the same estimator. `diversity_stats.R` prints both, as
`He_2n_corr` and `He_NeiChesser`. Quote `He_2n_corr` when comparing against
another study, `He_NeiChesser` when you want the better estimate, and say which.

### One parameter, three names

Before anything else, because the literature uses three words for two ideas and
this document would otherwise be unreadable.

There is **one parameter**: the probability that two randomly drawn gene copies
differ at a site. Everything below is that number, and the words differ only by
what you average it over and how you estimate it.

| what changes | names it goes by |
|---|---|
| **denominator: variant sites only** | expected heterozygosity, gene diversity, H<sub>e</sub> |
| **denominator: every sequenced site** | nucleotide diversity, π |
| **estimator: (2n/(2n−1))(1 − Σp²)** | what Stacks prints as `Pi`; unbiased only if F<sub>IS</sub> = 0 |
| **estimator: Nei & Chesser (1983)** | what `hierfstat` and this package use; unbiased at any F<sub>IS</sub> |

Two traps follow directly.

**π and H<sub>e</sub> are not different quantities, but they are different
numbers**, because the denominators differ by orders of magnitude. Ours were
0.34 per variant site and 0.00079 per sequenced site — the same parameter.

**Stacks' `Pi` column is a per-variant-site number**, despite the name. Its
per-sequenced-site counterpart is the `Pi` in the *All positions* block of
`sumstats_summary.tsv`. Those two also differ by orders of magnitude, and papers
routinely quote the first while calling it nucleotide diversity.

In this document and in the script output: **H<sub>e</sub>** always means per
variant record, **He_autosomal** always means per sequenced site, and π is used
only as a synonym for the second. Neither the script nor this document has a
column called `pi`.

Formulas for all of these, with each Stacks column matched to its counterpart
here, are in **"Formulas, and what maps onto what in Stacks"** near the end.
Every one is verified against real Stacks output by
`diversity_core.R --selftest` (section 10).

### Which file for which statistic

Step 1 gave you two files, and `populations` can also emit a third
(`--write-single-snp`, one SNP per RAD tag). Here is which to use, and why.

| statistic | file | reason in one line |
|---|---|---|
| H<sub>o</sub>, H<sub>e</sub>, and H<sub>e</sub> per sequenced site (= π) | **`.snps`** (all SNPs) | the only file whose numbers are on a scale another study can read |
| F<sub>IS</sub> | **`.haps`** | a ratio, so scale-free; unbiased and the most precise of the three |
| rarefied allelic richness | **`.haps`** | on biallelic SNPs it can only be 1 or 2 |
| rarefied private alleles | **`.haps`** | same, though it stays usable on SNPs |
| between-population test | either, but say which | it works on both |
| `--write-single-snp` | **not needed for any of this** | see "the thinning question" below |

Two runs, two jobs. The split is not arbitrary: **F<sub>IS</sub> and allelic
richness are the statistics that benefit from multi-allelic markers, and
H<sub>e</sub>/H<sub>o</sub>/π are the ones that need a fixed, interpretable
scale.** The rest of this section is why.

#### The reasoning, worked through

**Why haplotypes for F<sub>IS</sub>.** A RAD tag with three SNPs becomes one
marker with up to eight alleles instead of three markers with two alleles each.
More alleles means more information per locus, and it dodges a specific
biallelic pathology: at a SNP where the minor allele is present in a single
copy, that copy must sit in a heterozygote, so H<sub>o</sub> and H<sub>s</sub>
are forced equal and F<sub>IS</sub> is exactly 0 by arithmetic rather than by
biology. Multi-allelic loci have no equivalent trap.

Simulated with a known F<sub>IS</sub> of 0.10 (1,200 tags, n = 15 and 10, six
replicates), all three files recover it — the difference is precision:

| file | F<sub>IS</sub> | bias | width of the 95% CI | records |
|---|---|---|---|---|
| `.haps` | 0.0993 | −0.0007 | **0.0287** | 1,200 |
| `.snps`, all SNPs | 0.1010 | +0.0010 | 0.0321 | 2,247 |
| `.snps`, one per tag | 0.0996 | −0.0004 | 0.0339 | 1,200 |

One modelling note, because it is the thing that would have rigged this
comparison if done wrong: inbreeding makes the two chromosomes identical by
descent across the **whole tag**, not independently at each SNP. A simulator
that inbreeds each SNP separately understates haplotype homozygosity and makes
the haps file look worse than it is. IBD is drawn once per individual per tag
above.

**Why all SNPs for π.** Nucleotide diversity is a per-**site** quantity: sum
H<sub>e</sub> over sites, divide by every sequenced base. That only works when
one record is one site, which is true of the all-SNP file and false of the other
two. On `.haps` a record is a whole tag (the script now refuses the conversion
and says so). On a single-SNP file a record *is* a site, but you have thrown
away most of the variable sites while still dividing by all sequenced sites, so
you would understate π by roughly the fraction discarded.

**Why SNPs for H<sub>e</sub>, H<sub>o</sub> and π.** Start with the honest part:
for the *comparison between your populations*, the haplotype file is fine. With
popB simulated as a founder subset of popA, five replicates, pooled denominator:

| | H<sub>e</sub> popA | H<sub>e</sub> popB | ratio A:B |
|---|---|---|---|
| per-SNP (`.snps`) | 0.3050 | 0.1560 | 1.96 |
| haplotype (`.haps`) | 0.4893 | 0.2346 | 2.09 |

Both say popA is about twice as diverse. If that ratio were the only thing you
reported, either file would do.

The problem is the absolute number, and it is a problem of **scale**. Varying
SNP density — a purely technical property, set by read length, enzyme choice and
your SNP-calling filters — while holding the population model fixed:

| max SNPs/tag | mean SNPs/tag | per-SNP H<sub>e</sub> | haplotype H<sub>e</sub> |
|---|---|---|---|
| 2 | 1.64 | 0.2821 | 0.3776 |
| 4 | 2.72 | 0.3050 | 0.4893 |
| 6 | 3.75 | 0.3122 | 0.5361 |
| 8 | 4.75 | 0.3104 | 0.5575 |

Haplotype H<sub>e</sub> swings **48%**; per-SNP H<sub>e</sub> swings 10%. (The
simulation couples SNP count and haplotype count, so part of that 10% is genuine
diversity rather than artefact — read the ~5× difference in sensitivity, not the
absolute numbers.)

Haplotype gene diversity has no fixed ceiling: a tag with four SNPs can carry
eight haplotypes, a tag with one can carry two. So "H<sub>e</sub> = 0.49" from a
haplotype file is uninterpretable without knowing your read length, your enzyme
and your filters, and it is not comparable to another study — or to your own
data reprocessed with different Stacks parameters. Per-site H<sub>e</sub> is
bounded at 0.5 for a biallelic site and means the same thing in everyone's paper.

A secondary point, smaller than you might expect: a single base error anywhere
in a tag creates a spurious haplotype, whereas at SNP level it corrupts one
site. At a 1% per-SNP error rate H<sub>e</sub> inflated **+1.6% per-SNP versus
+2.5% haplotype**. Real, but not on its own a reason to choose.

And π is not a preference at all: nucleotide diversity is per-**site**, so only
a file where one record is one site can estimate it.

**Why F<sub>IS</sub> escapes the scale problem.** F<sub>IS</sub> is a ratio of
H<sub>o</sub> to H<sub>e</sub>. Both numerator and denominator inflate together
when a locus carries more alleles, so the ratio is scale-free — which is exactly
why the same argument that rules haplotypes out for H<sub>e</sub> leaves them as
the better choice for F<sub>IS</sub>.

**What this rules out.** Running `.haps` for everything. It would give you
H<sub>e</sub>, H<sub>o</sub> and π on a scale no reader can interpret. Running
`.snps` for everything is defensible if you want one file — you lose about 11%
precision on F<sub>IS</sub> and cannot report allelic richness, but every number
reconciles with every other and with Stacks.

#### The trap that will actually bite you

H<sub>o</sub> and H<sub>e</sub> are **not comparable between the two files.** In
the same simulated dataset, popA's H<sub>o</sub> was **0.366 from `.haps` and
0.299 from `.snps`**. Neither is wrong — they are different quantities.
Haplotype H<sub>o</sub> asks "is this tag heterozygous?", which is true if *any*
SNP on it is; per-site H<sub>o</sub> asks "is this site heterozygous?"

F<sub>IS</sub> is a ratio, so it is comparable across the two. H<sub>o</sub> and
H<sub>e</sub> are not. Put the file in the table caption or a column, because a
reader who checks 1 − H<sub>o</sub>/H<sub>e</sub> against your F<sub>IS</sub>
across mismatched files will find a discrepancy far larger than any of the
estimator issues in this document.

**Free diagnostic:** F<sub>IS</sub> from the two files should agree closely. If
it does not, that is informative — the usual suspects are sequencing error
creating spurious rare haplotypes (inflates haplotype H<sub>o</sub>, pushes haps
F<sub>IS</sub> down) or dropout removing whole tags. Run both and check.

#### The thinning question

`--write-single-snp` exists because SNPs on one RAD tag are physically linked
and are not independent replicates. That is a real problem, and it is the
objection a reviewer is most likely to raise.

It is already handled. `diversity_stats.R` bootstraps over **RAD loci**, not
over SNP rows — it takes the per-locus H<sub>o</sub> and H<sub>s</sub> and
resamples tags, so linked SNPs move together as a block. This is exactly why it
does not simply call `hierfstat::boot.ppfis()`, which resamples rows and gave
intervals about 1.8× too narrow on this project's data.

So thinning buys independence you already have, and charges you half your
records for it. The table above shows the result: same estimate, widest interval
of the three.

If you do thin for some other reason, prefer `--write-random-snp` to
`--write-single-snp`. The latter takes the first SNP on every tag, a fixed
position relative to the cut site rather than a random draw.

**Where a single-SNP file genuinely is required** — all of it outside this
workflow: STRUCTURE, DAPC and `find.clusters()`, PCA, F<sub>ST</sub> outlier
scans, and LD-based N<sub>e</sub>. Those either assume unlinked markers or, for
LD-N<sub>e</sub>, are actively biased downward by within-tag linkage. Generate
that file when you need it, and do not use it for diversity statistics.

### Why not just read Stacks' F<sub>IS</sub> column

F<sub>IS</sub> is a fraction — observed heterozygotes over expected — and you
have thousands of loci. There are two ways to average a fraction.

**Average the fractions.** Compute F<sub>IS</sub> at each locus, then average
those. **One big fraction.** Add up all observed heterozygotes, add up all
expected, divide once at the end.

Think of a team batting average. The second way adds every hit and every at-bat
and divides — correct. The first averages each player's personal average, so the
pitcher who batted twice counts as much as the guy who batted 600 times. Loci
differ enormously in how much information they carry, so the two answers are not
close.

**Stacks reports F<sub>IS</sub> the first way and H<sub>o</sub>/H<sub>e</sub> the
second way, so its own columns do not reconcile:**

| | Stacks H<sub>o</sub> | Stacks H<sub>e</sub> | F<sub>IS</sub> those imply | Stacks' F<sub>IS</sub> column |
|---|---|---|---|---|
| population 1 | 0.2488 | 0.2995 | **0.169** | **0.1375** |
| population 2 | 0.1707 | 0.1891 | **0.097** | **0.0446** |

A reviewer with a calculator finds that in ten seconds. Use the ratio of sums
(Weir & Cockerham 1984; Bhatia et al. 2013), which is what this script does. A
locus that is monomorphic in a population contributes 0 to both sums, so the
ratio doesn't care whether you include such loci — which is the deeper reason to
prefer it.

### Why not Stacks' `Pi` for H<sub>e</sub> either

#### Which Stacks column is which

`populations.sumstats.tsv` carries three columns that all sound like the same
thing. They are not interchangeable.

| column | formula | what it is |
|---|---|---|
| `Obs_Het` | count of heterozygotes ÷ individuals | a count. No estimator, nothing to get wrong. Safe to report. |
| `Exp_Het` | 1 − Σp̂² | the plug-in estimate, **biased low**. Don't report it and don't compare against it. |
| `Pi` | (2n/(2n−1))(1 − Σp̂²) | `Exp_Het` with the standard unbiased-π correction. **This is the one to report.** |
| `Fis` | mean of per-locus ratios | doesn't reconcile with this file's own `Obs_Het` and `Pi`. Don't report it. |

So `Pi` = `Exp_Het` × 2n/(2n−1) — not a rounding difference:

| n (diploids) | `Pi` ÷ `Exp_Het` |
|---|---|
| 2 | 1.333 |
| 10 | 1.053 |
| 15 | 1.034 |
| 25 | 1.020 |

**Verify it on your own build** rather than trusting this table: take any row of
`sumstats.tsv`, multiply `Exp_Het` by 2N/(2N−1) using that row's N, and it should
land on `Pi`. (Checked against a published Stacks summary file: `Exp_Het` 0.40000
at `Num_Indv` 2 is reported as `Pi` 0.53333, and 0.4 × 4/3 = 0.53333 exactly.)

**One naming trap.** At a *biallelic* site, π and expected heterozygosity are the
same parameter — the probability two randomly drawn gene copies differ. At a
*multi-allelic* locus they are not: gene diversity treats every pair of distinct
alleles as equally different, while π weights each pair by how many nucleotides
actually differ. Stacks keeps them apart in `populations.hapstats.tsv`, which has
**Gene Diversity** and **Haplotype Diversity** as separate columns. So on
haplotype data, nothing here is nucleotide diversity whatever it is called.

Finally: **never take F<sub>IS</sub> from the "All positions (variant and fixed)"
block.** Monomorphic sites have F<sub>IS</sub> = 0 and Stacks averages ratios, so
every fixed site drags the mean toward zero and you get a number an order of
magnitude too small. (The flip side of why a ratio of sums is right: a
monomorphic locus adds 0 to both sums and cannot move the answer at all.)

---

Returning to the estimator. `Exp_Het`'s sample-size correction — the one in `Pi` —
assumes the 2n
gene copies are an independent sample of gametes — i.e. that F<sub>IS</sub> = 0,
the very thing you are measuring. It therefore runs low whenever there is a
heterozygote deficit and drags F<sub>IS</sub> toward zero: about 3% at n = 15 and
5% at n = 10.

**The bias depends on sample size, so it does not cancel when you compare a
sample of 15 against a sample of 10.** That is what makes it fatal here rather
than merely untidy. The correct estimator subtracts one more term (Nei & Chesser
1983) and is what `hierfstat` uses.

If you have older numbers computed the other way, expect F<sub>IS</sub> to rise
by roughly 3–5% when you re-run. Both values are printed so old drafts stay
traceable.

### Denominators, and where nucleotide diversity fits

Averaged over your markers, H<sub>e</sub> is heterozygosity *per marker your run
happened to call*. Call more markers and the number changes. It is fine for
comparing your populations **to each other** and useless for comparing to anyone
else's paper. The fix is one multiplication and the script does it if you pass
`<sites>`:

> H<sub>e</sub>(autosomal) = H<sub>e</sub>(per marker) × markers ÷ sites sequenced

**On a SNP file, that per-site number is nucleotide diversity** — because both
conditions from the vocabulary box now hold: each record is a site, and the
denominator is every sequenced base. It is still gene diversity; π is the name
for gene diversity under exactly those two conditions.

The script's `He_autosomal` column is this Nei–Chesser `He` rescaled by
markers ÷ sites — it does not separately recompute the (2n/(2n−1)) estimator
at the autosomal scale. If you need that comparison, it is the
`He_NeiChesser` vs `He_2n_corr` table further up (see "One parameter, three
names"): it is a per-record comparison, but since the autosomal figure is just
that same `He` multiplied through by a constant, the same F/(2n−1) gap between
the two estimators carries over to the autosomal scale unchanged. Report one
estimator and name it.

**This only works on the SNP VCF**, because only there are markers the same
thing as sites. On the haplotype VCF, H<sub>e</sub> is per-tag gene diversity,
so markers ÷ sites gives "gene diversity per sequenced nucleotide" -- not π, and
an underestimate of it that grows with your SNPs-per-tag. The script detects
the file type and refuses the conversion rather than printing a number that
looks quotable. Hence the two invocations above.

Compare the SNP-VCF autosomal figure against Stacks' `Pi` from the **"All
positions (variant and fixed)"** block of `sumstats_summary.tsv`, not the
variant-positions block. Expect the script's value to sit slightly higher when
F<sub>IS</sub> > 0. Quote the autosomal figure next to published values (Schmidt
et al. 2021).

#### Reporting per polymorphic site instead

Most of the RADseq literature reports H<sub>e</sub>, H<sub>o</sub> and π per
*variant* site rather than per sequenced site, and there are good reasons to do
the same — comparability with prior work on your system, co-author expectations,
or an unreliable `Sites` count. That is fine. The script gives you both from the
same run. Two things to get right.

**Use the pooled denominator.** "Per polymorphic site" is ambiguous, and the two
readings give opposite answers:

- **Pooled** — sites variable somewhere in the dataset. Every population scored
  on the same denominator. This is what Stacks' "Variant positions" block does
  and what `diversity_stats.R` does.
- **Within-population** — for each population, average only over the sites where
  *that* population is variable.

Simulated with popB genuinely about 2.3× less diverse than popA:

| denominator | H<sub>e</sub> popA | H<sub>e</sub> popB | ratio A:B |
|---|---|---|---|
| pooled | 0.3512 | 0.1496 | **2.35** |
| within-population | 0.3591 | 0.3731 | **0.96** |

The within-population denominator does not merely shrink the difference, it
**reverses** it — popB comes out marginally more diverse. Across five seeds:
2.23–2.45 pooled, 0.96–1.02 within-population.

The mechanism is simple once seen. popB is monomorphic at 60% of sites, and that
*is* its reduced diversity. Dropping those sites from popB's denominator
discards the evidence and leaves only the sites where it happened to retain
variation, which look normal. It conditions on the outcome.

**Report `pct_poly` beside it.** The two denominators are related by an exact
identity:

> H<sub>e</sub>(pooled) = H<sub>e</sub>(within-population) × fraction of sites
> polymorphic in that population

0.3731 × 0.401 = 0.1496, exactly. So a table giving the pooled H<sub>e</sub>
*and* `pct_poly` — which the script already prints — lets a reader recover
either version and cannot mislead anyone. It also hands you a second,
independent diversity statistic for free.

**And do not thin.** Once the denominator is "SNPs used," one-random-SNP-per-tag
looks like it should be safe: a random subset estimates the mean of the whole.
It is not a random subset of *SNPs* — it gives every tag equal weight regardless
of how many SNPs it carries, and SNP-rich tags are more diverse tags:

| | H<sub>e</sub> popA | H<sub>e</sub> popB | ratio |
|---|---|---|---|
| all SNPs | 0.3512 | 0.1496 | 2.35 |
| one random SNP per tag | 0.3511 | 0.1317 | 2.67 |

popA barely moves (0–2% across seeds); popB drops 11–15%. The bias is
**asymmetric between populations**, so it inflates the headline ratio by about
14%. A symmetric bias could be argued away; this one cannot. Its size depends on
how strongly SNP density correlates with diversity in your data, but the
direction is generic.

A caption that settles all of this: *"H<sub>e</sub> averaged over the N SNPs
variable in at least one population; % polymorphic gives the fraction variable
within each population."*

### Reading the output

You get H<sub>o</sub>, H<sub>e</sub>, F<sub>IS</sub>, rarefied allelic richness,
rarefied private allelic richness, and the dataset-wide private-allele total
(`priv_total`) per population, each with a 95% interval from a bootstrap over
RAD loci AND a jackknife standard error (the `_se` and `_lo`/`_hi` columns
alongside each estimate, in both the printed table and the
`diversity_per_population*.tsv` / `diversity_richness*.tsv` files — see
"Bootstrap mode" and "Standard errors" below for what each one is and how
they were cross-checked against each other). % polymorphic is reported
alongside but has no interval or SE of its own.
`priv_total` is `NA` (not zero) for a population with no defined locus at all
— that is "no information," not "no private alleles"; a population that DOES
have defined loci and simply carries zero private alleles at every one of them
gets a `priv_total` of 0, correctly distinguished from the `NA` case.

**Those intervals answer "what if I had genotyped different markers in these same
animals?"** That is a legitimate question and it is not the one you are asking.
Use them to describe uncertainty in a single population's estimate. Do **not**
use them to test whether two populations differ — that is Step 3, and the reason
is the whole of Step 3.

### Bootstrap mode: `--boot=loci` (default), `individuals`, `both`

Every one of the intervals above comes from resampling RAD loci with
replacement (`--boot=loci`, the default and the only mode enabled by
default). Two other modes exist for explicit comparison:

- `--boot=individuals` resamples individuals within each population instead
  (the convention used by `diveRsity::divBasic`, and the scheme Van Dongen
  1995 and Petit & Pons 1998 argue for on the grounds that loci measured on
  the same individuals are not independent replicates).
- `--boot=both` resamples both axes at once, independently, per replicate —
  the crossed-factor "product weight" bootstrap of Owen & Eckles (2012,
  *Annals of Applied Statistics* 6:895–927), built to handle exactly this
  kind of two-way (locus × individual) structure. McCullagh (2000,
  *Bernoulli* 6:285–301) proved no *exact* bootstrap exists for a crossed
  design like this; Owen & Eckles show that independently reweighting each
  axis and multiplying the two weights together is provably conservative
  (never anti-conservative) for a linear mean, and cite Hall (1992) and
  Mammen (1992) for extending that guarantee, via the delta method, to
  smooth functions of a resampled mean.

**Neither of those two is the default, and neither should be treated as a
valid confidence interval, because we checked and they aren't.** Before
shipping `both` as the new default (the original plan for this feature), an
explicit coverage simulation was run: synthetic populations with a KNOWN
true H<sub>e</sub>, F<sub>IS</sub> and A<sub>r</sub> were bootstrapped under
each mode across 150–300 independently simulated datasets (n = 10 and 15
individuals, matching this project's typical population sizes; L = 1,500
loci), and what was measured was not "does the theory sound reasonable" but
"does the resulting 95% CI actually contain the true value 95% of the time."

```
                          loci      individuals   both
He  coverage (n=10)       93.7%       0.0%        0.0%
He  coverage (n=15)       92.7%       0.0%        0.0%
Fis coverage (n=10)       94.0%       4.0%       13.7%
Fis coverage (n=15)       94.0%       4.7%       22.7%
Ar  coverage (g = 2n)      ~95%*      0.0%        0.0%
Ho  coverage (n=10,15)     n/a       88–92%  (Ho has no correction factor;
                                      point estimate inside its own CI 100%)
```
*The `loci`-mode Ar number is coverage of its own SAMPLE-CONDITIONAL
estimand (Kalinowski's formula asks "if I drew g gene copies from the SAME
2n gene copies I already observed, without replacement, how many distinct
alleles would I expect" — not "g draws from an infinite source
population"). Those two questions coincide only when the observed 2n gene
copies are themselves 2n independent draws, which real inbreeding breaks
(identical-by-descent copies are correlated, not independent) — checked
directly by re-running the simulation with and without inbreeding in the
generative model at the identical g = 2n case: without it, coverage was
~95%; with it, only 52–71%. Neither number is "wrong" — they're coverage
of different questions. This is a limitation of what A<sub>r</sub> means
under inbreeding, not a defect in resampling loci, and it was true before
this session's changes (`--boot=loci`'s A<sub>r</sub>/privA<sub>r</sub>
arithmetic is unchanged, pre-existing code).

**Why**: H<sub>e</sub> (via Nei & Chesser's `(n/(n-1))` term) and rarefied
richness (via a hypergeometric formula keyed on total gene copies) both
assume the individuals contributing to a replicate are *distinct*. A
with-replacement resample of individuals typically contains only about 63%
distinct individuals — duplicates take the place of others — and the
correction factor has no way to know that. This pushes the resampled
H<sub>e</sub>/F<sub>IS</sub>/A<sub>r</sub> systematically LOW. Confirmed
directly, not just theorized: forcing the resampling weights to be all-1
(no duplication) reproduces the point estimate to 13 decimal places, so the
arithmetic is correct and the effect is a genuine property of resampling
with replacement; H<sub>o</sub>, which carries no such correction factor, is
unaffected (88–92% coverage, and its point estimate always falls inside its
own interval) — direct confirmation that the mechanism is specific to
estimators with a built-in finite-sample correction, not to individual
resampling generally.

**This gets WORSE, not better, with more loci — a real problem for RAD-seq
data specifically.** The bias above is a property of the individual axis
alone and does not shrink as locus count grows; the bootstrap's spread,
however, does shrink as locus count grows (more loci means a more precise
estimate of wherever the resampling distribution is centered). So the
bias-to-spread ratio grows, not shrinks, with more loci — meaning a typical
RAD-seq dataset (hundreds to thousands of loci) sits in the worst part of
this problem, not a forgiving corner of it. `--boot=both` does not rescue
this: it resamples individuals too, and inherits the same bias, which the
added locus-resampling variance does not reliably cover (0–23% in the table
above). Owen & Eckles' guarantee is real, but it is proved for a linear
sample mean; H<sub>e</sub> and rarefied richness are not smooth functions of
a resampled mean alone; they have an explicit, separate dependence on
sample size that the delta-method extension does not reach. Earlier drafts
of this feature claimed that extension covered F<sub>IS</sub>/A<sub>r</sub>
— that claim did not survive the coverage simulation and is corrected here.

**What each mode is actually for**: `--boot=loci` is the only mode to
quote in a methods section. `--boot=individuals` and `--boot=both` exist so
you can directly compare this pipeline's numbers against a diveRsity- or
Owen-Eckles-style analysis someone else ran, or reproduce what those tools
would report — a legitimate thing to want, just not a legitimate substitute
for the default when reporting your own confidence intervals. Filenames get
a `.boot-individuals`/`.boot-both` suffix only when the flag is passed
explicitly, so a comparison run cannot silently overwrite the default
output.

### Per-metric summary, with `hierfstat` and `diveRsity` as reference points

Every metric in this pipeline gets the same answer — `--boot=loci` — but
the reasoning differs by metric, and it's worth seeing what the two most
common alternative tools actually do, verified directly from their own
source code (both are re-verifiable: `hierfstat` is on CRAN and its
functions can be printed directly in R; `diveRsity`'s source is on GitHub
at `kkeenan02/diveRsity`).

`boot.ppfis()` is `hierfstat`'s ONLY bootstrap function that touches any of
these statistics, and it returns an F<sub>IS</sub> confidence interval
only — it has no H<sub>o</sub>/H<sub>e</sub> bootstrap at all (confirmed by
printing the function: its output is `1 - colSums(Ho[x,])/colSums(Hs[x,])`
and nothing else is returned).

| Metric | Correction factor baked in? | This pipeline | `hierfstat` | `diveRsity::divBasic` |
|---|---|---|---|---|
| H<sub>o</sub> | No | `--boot=loci` | No H<sub>o</sub> bootstrap exists (`boot.ppfis()` returns F<sub>IS</sub> only) | Resamples **individuals** (`sample(dim(pasub)[1], replace=TRUE)`) |
| H<sub>e</sub> | Yes (Nei-Chesser `n/(n-1)`) | `--boot=loci` | No H<sub>e</sub> bootstrap exists, same reason | Resamples individuals, but on the **uncorrected** `1-sum(p^2)` estimator — no `n/(n-1)` term anywhere, confirmed in both its point estimate and its bootstrap code. This means diveRsity's individual-resampling isn't actually a counterexample to the finding above: it never had a correction factor for duplication to interact badly with. Its H<sub>e</sub> is a different, already-known-biased-low quantity regardless of resampling axis — this document already says to compare against Stacks' `Pi`, never an uncorrected `Exp_Het`-style number like diveRsity reports. |
| F<sub>IS</sub> | Yes (built from H<sub>e</sub>) | `--boot=loci` | `boot.ppfis()` resamples loci, but per-SNP **row**, not per-RAD-**block** (`x <- sample(nloc, replace=TRUE)`) — the exact mechanism behind this project's own "~1.8x too narrow" finding above | Resamples individuals; same uncorrected-H<sub>e</sub> caveat as above |
| A<sub>r</sub> / privA<sub>r</sub> / `priv_total` | Yes (hypergeometric, keyed on total gene copies) | `--boot=loci` | `allelic.richness()` has **no bootstrap at all** — point estimate only | `rarefactor()` computes point-estimate rarefied richness only — **no bootstrap, no private-allele richness at all** |

Two things this table means in practice:
- "Matches hierfstat" (as this document says elsewhere) means matching its
  resampling *axis* (loci), not its resampling *unit* — hierfstat itself
  resamples individual SNP rows, which is the within-tag-linkage problem
  this pipeline's block bootstrap exists to avoid.
- Neither tool bootstraps A<sub>r</sub>/privA<sub>r</sub> at all, so there
  is no external tool to benchmark this pipeline's richness CI against.

### Standard errors: the delete-one-block jackknife over loci

Alongside every `_lo`/`_hi` pair, both tables also report a `_se` column
(`Ho_se`, `He_se`, `Fis_se`, `Ar_se`, `privAr_se`, `priv_total_se`). This
is **not** derived from the bootstrap replicates — there is no simple
closed-form variance formula for a ratio-of-sums statistic over linked RAD
loci (that's the whole reason population genetics uses resampling here
instead of a formula in the first place). It's a second, independent
method: the **delete-one-block jackknife**, the standard population-
genetics alternative to bootstrapping (Weir 1996, *Genetic Data Analysis
II* — the same method behind GENEPOP's and FSTAT's own F<sub>ST</sub>/F<sub>IS</sub>
standard errors), computed here over whole RAD-locus blocks for the same
linkage reason the bootstrap does.

```
SE_jackknife = sqrt( (nL-1)/nL × Σ_b (θ̂₍₋b₎ − θ̄)² )
```

where `θ̂₍₋b₎` is the statistic recomputed with RAD-locus block *b* dropped
entirely, and `nL` is the number of RAD loci. Unlike a with-replacement
bootstrap, a jackknife replicate never duplicates anything — it only ever
removes one block — so it does **not** carry the compositional-duplication
bias documented above for individual-resampling; it uses the same `loci`
axis already validated by the coverage simulation. It needs no random
numbers and no `--nboot`: with `nL` RAD loci there are exactly `nL`
delete-one replicates, so the `_se` columns are populated even at
`--nboot=0` and are identical run to run.

**Cross-checked against the bootstrap's own spread**, at two very
different locus counts, because a jackknife's reliability is known to
depend on how many blocks it has: on the 304-RAD-locus test fixture used
throughout this section, `1.96 × SE_jackknife` matched the `--boot=loci`
bootstrap CI's half-width to within about 1% for every metric
(H<sub>e</sub>, F<sub>IS</sub>, A<sub>r</sub>, privA<sub>r</sub>,
`priv_total`); re-checked on a much smaller 60-RAD-locus fixture with an
n=2 population (about as few blocks and as few individuals as this
pipeline allows), the same comparison held within about 2–6% for every
metric — still close agreement, not a divergence, though naturally less
tight than at 304 loci. Two independently-computed quantities agreeing
closely at both ends of the locus-count range this pipeline supports is
the actual reason to trust either of them.

**One more thing this jackknife work caught, fixed here rather than left
for later**: `Ar`/privA<sub>r</sub>'s point estimate (and everything built
from it, including the new jackknife SE) printed the literal text `NaN`
instead of a clean `NA` when a population had it undefined at every single
retained locus (the same degenerate case Bug A's fix, earlier in this
project's history, already made a NOTE about) — `colMeans(..., na.rm=TRUE)`
on an all-`NA` column is `0/0`, which R evaluates to `NaN`, not `NA`.
`priv_total` already had this guarded (a population with no defined locus
gets `NA`, not `NaN` — see above); `Ar`/privA<sub>r</sub> did not, and now
do, using the identical guard. This was pre-existing behavior this session
did not introduce, caught only because validating the new jackknife SE's
NA-handling required constructing exactly the input that triggers it.

### How does this compare to Stacks' own `StdErr` column?

Stacks' `populations.sumstats_summary.tsv` reports its own `StdErr` for
`Obs_Het`, `Exp_Het`, `Fis` and others. Verified directly from Stacks
2.68's own source (`src/populations.cc`,
`SumStatsSummary::final_calculation()`), not inferred from its
documentation: it computes the classical

```
StdErr = sqrt(sample variance) / sqrt(n)
```

with the variance accumulated (via Welford's online algorithm) over
**every variant SNP site** — not RAD-locus blocks. That is the identical
per-row, not-per-block treatment already identified above as the mechanism
behind hierfstat's `boot.ppfis()` giving intervals "~1.8x too narrow" (see
"The thinning question"): Stacks' own `StdErr` should be expected to be
optimistically narrow for the same reason, on the same kind of data — not
a coincidentally different number from a different bug. Stacks' `Fis`
mean is also accumulated as a running average of the **per-site Fis
ratio** — the mean-of-per-locus-ratios aggregation this document already
argues against in favor of the ratio-of-sums (see "Formulas", and Weir &
Cockerham 1984; Bhatia et al. 2013), independent of the SE question.

**A further, separate finding, reported factually rather than as an
opinion**: in that same Stacks 2.68 source, the `Fis` `StdErr` output line
reads

```cpp
<< sqrt(_num_indv_var[j]) / _sq_n[j] << "\n";   // should be _fis_var[j]
```

— a copy-paste error, so the `Fis` `StdErr` column in
`populations.sumstats_summary.tsv` is actually `Num_Indv`'s standard
error, not Fis's. This was checked in one specific version (2.68); it is
worth confirming against whichever version is actually installed rather
than assumed to hold for every release, but the site-vs-block and
mean-of-ratios issues above are structural, not version-specific.

**So: this project's jackknife SE is not "similar" to Stacks' `StdErr`.**
It differs in the two dimensions that matter — delete-one-**block**
(RAD locus) replication vs delete-one-**site**, and a delete-one recompute
of the actual ratio-of-sums estimator vs a classical SD/√n applied to raw
per-site ratios — for the same reason the rest of this project's bootstrap
machinery already differs from the field's per-row conventions.

---

## Step 3 — test whether the populations actually differ

```bash
Rscript het_between_pops.R out/populations.snps.vcf popmap.tsv 0.9 het_out
```

Arguments: `<vcf> <popmap> [min_call] [outdir]`. `min_call` is the minimum
per-individual genotyping rate for a locus to be included (0.9 is a sensible
default).

### Why the obvious test is wrong

The standard approaches — bootstrap over loci, or a paired Wilcoxon on per-locus
values — treat **loci** as the replication unit. Both are common in the
literature. Both are badly wrong for this question.

Measured type I error, two populations with **identical** true heterozygosity, so
every "significant" result is a false alarm. n = 15 and 10, 2000 loci. Target 5%:

| how much individuals differ from each other | bootstrap over loci | paired Wilcoxon over loci | **Welch *t* on individuals** |
|---|---|---|---|
| not at all | 5% ✓ | 4% ✓ | 5% ✓ |
| a little (SD of *F* = 0.10) | **55%** ✗ | **54%** ✗ | 4% ✓ |
| a fair amount (SD = 0.25) | **76%** ✗ | **75%** ✗ | 5% ✓ |

The mechanism is pseudoreplication. The quantity you care about is a *population
mean*, and its uncertainty is dominated by **which animals you caught**, not
which markers you typed. Loci are repeated measurements on the same individuals.
Adding loci does not reduce the real uncertainty at all — but a locus bootstrap
interval keeps shrinking as 1/√(number of loci), so with thousands of loci it
converges on a width of nearly zero and rejects almost always.

And individual heterozygosity being correlated across loci is not an exotic
condition. It is exactly what inbreeding produces. It is also what one
poor-quality library produces. H<sub>e</sub> fails the same way through
relatedness rather than inbreeding: 5% type I error with unrelated individuals,
19% with 3–4 per family, 66% with 5 per family.

**This is not a new objection.** Van Dongen (1995, *Heredity* 74:445–447)
concluded that resampling over loci does not satisfy the bootstrap's assumptions
and should be avoided, and suggested resampling over individual genotypes
instead. Petit & Pons (1998) repeat it; Neff (2010) wrote software for both. Most
papers ignore it.

### The fix

Give every individual **one** number — its proportion of heterozygous loci — and
you have an ordinary two-sample problem with no bootstrap at all. Welch's *t*
held its nominal rate at every level tested and is the most powerful of the
options; Wilcoxon is reported alongside as a distribution-free check.

This needs the VCF. Per-individual heterozygosity cannot be recovered from
`sumstats.tsv`, which holds only per-population counts.

### Reading the output — three things, in order

**1. The missingness confound check.** The script correlates each individual's
call rate with its heterozygosity. If that correlation is positive and
significant, individuals with more missing data look less heterozygous, which is
the allele-dropout signature — you would be measuring library quality, not
biology. Re-run with `min_call 1.0` and see whether the difference survives. It
also warns if mean call rate differs between populations by more than 2
percentage points.

**2. The overdispersion factor.** How much more individuals vary than pure
chance predicts. Near 1 means a locus bootstrap would have been correctly
calibrated on *your* data and the two approaches should agree. Above about 5
means a locus bootstrap understated the standard error by roughly √(that factor)
and its p-values were far too small. This is how you answer a reviewer who asks
why you didn't do it the usual way — you measured it on their data.

**3. The test itself.** Welch *t* with a difference, a confidence interval and
Hedges' *g* (≈0.2 small, 0.5 medium, 0.8 large), plus the Wilcoxon check.

Two caveats the script prints and you should carry into the text:

- **Power is limited by the number of individuals, not the number of loci.** At
  n = 15 vs 10 with realistic individual variation, simulation gave 52% power for
  a moderate difference. So a non-significant result is *weak evidence of no
  difference*, not evidence of no difference. Write "no difference was detected."
- **If some of your individuals are relatives they are not independent units**
  and even this test is anti-conservative. Check relatedness first. At small n,
  one bad library or one outlier can drive the entire result — look at the
  per-individual table before you trust anything.

Outputs: `individual_heterozygosity.tsv` and `het_between_pops_tests.tsv`, in
the directory you named.

`diversity_stats.R` tags its output files with the input stem — each run
writes `diversity_per_population.<stem>.tsv` and `diversity_richness.<stem>.tsv`
(so `.haps.tsv` and `.snps.tsv` pairs, four files total across both runs) —
so the second run does not overwrite the first.

---

## Step 4 — what to write

### Methods (adapt the numbers)

> Loci were exported with Stacks 2 `populations` requiring a minimum genotype
> depth of 10 reads (`--min-gt-depth 10`), retaining loci genotyped in ≥80% of
> individuals within each population and present in both populations
> (`-r 0.8 -p 2`); a per-population rather than pooled threshold was used so that
> the two samples were filtered symmetrically. A minor allele count of ≥3 was
> required. No Hardy–Weinberg filtering was applied, as a heterozygote deficit is
> the quantity under study.
>
> Gene diversity was estimated with the sample-size-unbiased estimator of Nei &
> Chesser (1983), and F<sub>IS</sub> as 1 − ΣH<sub>o</sub>/ΣH<sub>s</sub> summed
> over loci (Weir & Cockerham 1984), rather than as a mean of per-locus ratios.
> H<sub>e</sub> is additionally reported per sequenced site (Schmidt et al. 2021).
>
> Because heterozygosity is correlated across loci within an individual, loci are
> not independent replicates for a between-population comparison (Van Dongen
> 1995). Populations were therefore compared using Welch's *t*-test on
> per-individual heterozygosity, with a Wilcoxon rank-sum test as a
> distribution-free check. The overdispersion of individual heterozygosity
> relative to binomial expectation was [X], indicating that a locus-based
> resampling test would have understated the standard error by roughly [√X]-fold.

### Results

State the two F<sub>IS</sub> values, the *t* statistic, the p-value and the
confidence interval on the difference. If it is not significant, say **"no
difference in heterozygosity was detected between populations (Welch's t = …,
p = …, 95% CI on the difference … to …)"** — not "the populations do not differ."
With n = 15 and 10 you have roughly 50% power for a moderate effect and you
should say so.

### Which number to take from where

You do not have to recompute everything. Only F<sub>IS</sub> and the rarefied
statistics actually need this package.

| statistic | which file | Stacks OK? | source |
|---|---|---|---|
| H<sub>o</sub> | `.snps` | **yes** — it's a count, no estimator involved | `sumstats_summary.tsv`, `Obs_Het` |
| H<sub>e</sub> (gene diversity, per variant site) | `.snps` | **yes** — off by only F/(2n−1), under 0.5% here | `sumstats_summary.tsv`, `Pi`, Variant-positions block |
| π (per sequenced site) | `.snps` | **yes**, and it is the estimator other papers used | `sumstats_summary.tsv`, `Pi`, All-positions block |
| % polymorphic | `.snps` | yes | `sumstats_summary.tsv` |
| F<sub>IS</sub> | **`.haps`** | **no** — two separate errors, 3–5% and larger | `diversity_stats.R` |
| rarefied allelic richness | **`.haps`** | not computed by Stacks | `diversity_stats.R` |
| rarefied private allelic richness | **`.haps`** | **no** — Stacks' `Private` is a raw count, not rarefied, so it scales with sample size | `diversity_stats.R` |

If you are taking some numbers from Stacks and some from here, the
column-by-column equivalences are in "Formulas, and what maps onto what in
Stacks".

Put the file in the table caption or a column. H<sub>o</sub> and H<sub>e</sub>
are **not comparable between the two files** — 0.366 vs 0.299 for the same
population in simulation — so a reader checking 1 − H<sub>o</sub>/H<sub>e</sub>
against a F<sub>IS</sub> from the other file will find a discrepancy far larger
than any estimator issue in this document. F<sub>IS</sub> itself is a ratio and
is comparable across the two, which is the free cross-check below.

Use the **"All positions (variant and fixed)"** block for H<sub>o</sub> and `Pi`
if you want values comparable to published estimates, and the **"Variant
positions"** block for per-ascertained-SNP values — but say which, because they
differ by a factor of hundreds and papers routinely quote the second while
implying the first.

**Guard against the reconciliation problem.** You came here because Stacks' own
columns don't reconcile; mixing sources can recreate that. A reader computing
1 − H<sub>o</sub>/`Pi` from your Stacks numbers should land near your reported
F<sub>IS</sub>. It will, to within 3–5%, **provided the marker sets match** — and
it won't if you print Stacks' per-SNP numbers beside a F<sub>IS</sub> from the
haplotype run. Either report everything from the SNP VCF and use the haplotype
run only for allelic richness (labelled as such), or report the script's own
H<sub>o</sub> and H<sub>e</sub> and note that they agree with Stacks' `Obs_Het`
and `Pi` to within X%. `diversity_stats.R` prints the records it used; that
should be close to Stacks' `Variant_Sites`. If it isn't, your discrepancy is the
marker set, not the estimator.

The script prints a `He_NeiChesser` vs `He_2n_corr` comparison so you can quote
the actual percentage rather than the one in this table.

### One sentence on null alleles

RAD data carry restriction-site null alleles, which inflate F<sub>IS</sub>
because a chromosome with a mutation in the cut site produces no reads at all and
its carrier is called homozygous. Depth filtering cannot fix this. Add:

> F<sub>IS</sub> estimates should be read as upper bounds, since restriction-site
> null alleles inflate apparent homozygosity in RAD data and cannot be removed by
> depth filtering.

And if you are reporting Stacks values for everything else:

> Observed heterozygosity and nucleotide diversity are reported as calculated by
> Stacks. F<sub>IS</sub> was recalculated as 1 − ΣH<sub>o</sub>/ΣH<sub>s</sub>
> using the estimator of Nei & Chesser (1983), because the correction underlying
> Stacks' `Pi` assumes F<sub>IS</sub> = 0 and its F<sub>IS</sub> column is a mean
> of per-locus ratios rather than a ratio of sums. The two H<sub>e</sub>
> estimators differ by less than 1% at these sample sizes, so nucleotide
> diversity is unaffected; F<sub>IS</sub> is not.

That is sufficient, because null alleles inflate F<sub>IS</sub> in **both**
populations. A shared upward bias cannot manufacture a difference between them —
it can only make both numbers too high together. Since the comparison is your
claim, an upper bound supports it fully. (If you had found a *significant*
difference, this would need more work: the two populations could carry different
null loads. You did not.)

---

## Formulas, and what maps onto what in Stacks

Everything in this section is confirmed against real Stacks output, not
inferred. `diversity_core.R --selftest` reproduces the checks; the source rows are
from the Galaxy stacks2 test data, and the F<sub>IS</sub> definition is
Catchen's own, posted to the stacks-users list ("the calculation for Fis uses pi
for expected heterozygosity: Fis = (pi − obs_het) / pi", pointing to equation 2
of Hohenlohe et al. 2010).

### Variables

Everything below is *per site, per population*, unless stated.

| symbol | meaning | where to see it |
|---|---|---|
| *N* | diploid individuals genotyped here | `N` column of `sumstats.tsv` |
| *n* | **gene copies** = 2*N* | — |
| *n<sub>i</sub>* | copies of allele *i*, so Σ*n<sub>i</sub>* = *n* | — |
| *p<sub>i</sub>* | *n<sub>i</sub>* / *n* | `P` column is the commoner allele's *p* |
| *H*<sub>obs</sub> | heterozygous individuals / *N* | `Obs Het` |
| *L* | records used (SNPs or tags) | printed by the script |
| *S* | total sequenced sites | `Sites`, All-positions block |

The single most common slip is reading *n* as individuals. Every correction
below is in **gene copies**: at *N* = 15 the factor is 30/29, not 15/14.

### Per site

| quantity | formula | Stacks | this package |
|---|---|---|---|
| observed heterozygosity | *H*<sub>obs</sub> | `Obs Het` | `Ho` |
| gene diversity, plug-in | 1 − Σ*p<sub>i</sub>*² | `Exp Het` | not reported |
| gene diversity, *n*-corrected | (*n*/(*n*−1))(1 − Σ*p<sub>i</sub>*²) | **`Pi`** | `He_2n_corr` |
| gene diversity, Nei & Chesser | (*N*/(*N*−1))(1 − Σ*p<sub>i</sub>*² − *H*<sub>obs</sub>/*n*) | not available | **`He`** |
| inbreeding coefficient | 1 − *H*<sub>obs</sub>/(denominator) | `Fis`, using **`Pi`** | `Fis`, using `He` |

Two identities worth knowing:

**`Pi` = `Exp Het` × *n*/(*n*−1).** Catchen's combinatorial form,
1 − Σ C(*n<sub>i</sub>*,2)/C(*n*,2), is algebraically identical — both verified
on real rows. The factor is 1.034 at *N* = 15 and 1.053 at *N* = 10.

**Stacks' per-SNP `Fis` = (`Pi` − `Obs Het`)/`Pi`.** Not `Exp Het`. Confirmed:
at *N* = 2, *p* = 0.5, `Obs Het` = 1, Stacks reports `Fis` = −0.5, and
1 − 1/0.66667 = −0.5, whereas 1 − 1/0.5 would be −1.0.

So Stacks' F<sub>IS</sub> denominator is the estimator that assumes
F<sub>IS</sub> = 0. That is the first of its two errors.

### Across sites

| quantity | formula | Stacks | this package |
|---|---|---|---|
| mean gene diversity | (1/*L*) Σ<sub>*l*</sub> He<sub>*l*</sub> | `Pi`, Variant-positions block | `He` |
| **nucleotide diversity** | (1/*S*) Σ<sub>all sites</sub> He<sub>*l*</sub> | **`Pi`, All-positions block** | `He_autosomal` |
| F<sub>IS</sub> | (1/*L*) Σ<sub>*l*</sub> F<sub>IS,*l*</sub> — **mean of ratios** | `Fis` | — |
| F<sub>IS</sub> | 1 − (Σ<sub>*l*</sub> Ho<sub>*l*</sub>)/(Σ<sub>*l*</sub> He<sub>*l*</sub>) — **ratio of sums** | — | **`Fis`** |

Monomorphic sites contribute 0 to both sums, so the ratio of sums is unchanged
by whether you include them — which is why the variant/all-sites question does
not touch F<sub>IS</sub> at all. A mean of ratios has no such property: every
fixed site adds a 0 to the average and drags it toward zero, which is why
F<sub>IS</sub> must never be read from the All-positions block.

The script's autosomal figure is `He_autosomal` = (mean `He` over the loci
actually used) × *n*<sub>rec</sub> / *S*, where *n*<sub>rec</sub> is **all**
variant records, not just the used subset. Scaling by the subset instead would
understate π by the retention fraction. Under the default available-data mode
the used subset is close to all *n*<sub>rec</sub> records for most real
datasets, which is why `He_autosomal` now tracks Stacks' `Pi` closely (see
below); under `--complete-case` it can be a small, potentially unrepresentative
sliver, in which case this multiplication is a correction on top of a shakier
mean and the two numbers can diverge substantially (observed ~2× on real data
at ~6% retention).

### So: which numbers are actually interchangeable

| | Stacks | this package | same? |
|---|---|---|---|
| H<sub>o</sub> | `Obs Het` | `Ho` | **yes** in default (available-data) mode — a count either way, using the same available individuals; **not necessarily** under `--complete-case`, whose smaller record subset can be unrepresentative |
| gene diversity per variant site | `Pi` (Variant block) | `He` | ~0.3% apart in default mode: same parameter, different estimator; can diverge much further under `--complete-case` for the same reason as `Ho` |
| **nucleotide diversity** | **`Pi` (All-positions block)** | `He_autosomal` | ~0.3–1% apart in default mode (confirmed on real virgata data); can be off by ~2× under `--complete-case` if retention is low |
| F<sub>IS</sub> | `Fis` | `Fis` | **no, in either mode** — two independent errors in Stacks' `Fis` (biased denominator + mean-of-ratios aggregation), unrelated to missing-data handling |
| allelic richness | not computed | `Ar` | — |
| private alleles | `Private` (raw count) | `privAr` (rarefied) | **no** — the raw count scales with *N* |

### If you report π from Stacks

Reasonable, and this is how to do it cleanly.

Take **`Pi` from the "All positions (variant and fixed)" block** of
`populations.sumstats_summary.tsv` — not the Variant-positions block, which is a
per-variant-site number that will be hundreds of times larger and is not
nucleotide diversity. Take `Obs_Het` from the same block so the two share a
denominator. Take F<sub>IS</sub>, A<sub>r</sub> and private A<sub>r</sub> from
`diversity_stats.R`.

The one cost: the numbers no longer reconcile arithmetically, because Stacks'
π uses the *n*-corrected estimator and the script's F<sub>IS</sub> uses Nei &
Chesser. The gap is F/(*n*−1) — about 0.3% at your sample sizes, and the script
prints it for your data in the `He_NeiChesser` vs `He_2n_corr` table. Quote that
figure once in the methods and the discrepancy is documented rather than
mysterious:

> Nucleotide diversity and observed heterozygosity are reported as calculated by
> Stacks (`Pi` and `Obs_Het`, all positions), where π is estimated as
> (*n*/(*n*−1))(1 − Σ*p<sub>i</sub>*²) over *n* = 2*N* gene copies. F<sub>IS</sub>
> was recalculated as 1 − ΣH<sub>o</sub>/ΣH<sub>s</sub> using the estimator of
> Nei & Chesser (1983), because Stacks divides by π — which assumes
> F<sub>IS</sub> = 0 — and averages per-locus ratios rather than taking a ratio
> of sums. The two gene-diversity estimators differ by X% here, so π is
> unaffected; F<sub>IS</sub> is not.

### Check it on your own build

Version differences are real, so verify rather than trust this page. Take any
row of `populations.sumstats.tsv` and confirm all three:

```
Pi   ==  Exp Het * 2N/(2N - 1)
Fis  ==  (Pi - Obs Het) / Pi
Ho + Obs Hom  ==  1
```

If any fails, your Stacks differs from the version documented here and nothing
above should be trusted until you know how.

## If a reviewer pushes back

Every one of these is a reasonable-sounding objection with a short, measured
answer. Keep the numbers, not just the argument.

**"SNPs within a RAD locus are linked. You should have used one SNP per
locus."**

> Agreed that they are linked, which is why the confidence intervals are
> bootstrapped over RAD loci rather than over SNPs — linked SNPs are resampled
> together as a block. Thinning to one SNP per locus was tested and returns the
> same point estimate with a wider interval (0.0339 vs 0.0321 in simulation at a
> known F<sub>IS</sub> of 0.10), because it discards roughly half the data to
> buy independence the block bootstrap already provides.

**"Why doesn't your F<sub>IS</sub> match the F<sub>IS</sub> Stacks reports?"**

> Stacks computes F<sub>IS</sub> as a mean of per-locus ratios while computing
> H<sub>o</sub> and H<sub>e</sub> as ratios of sums, so its own columns do not
> reconcile with each other: its H<sub>o</sub> and `Pi` for these populations
> imply 0.169 and 0.097, while its F<sub>IS</sub> column reads 0.1375 and
> 0.0446. We use the ratio of sums throughout (Weir & Cockerham 1984), which is
> also what `hierfstat` does.

**"Your H<sub>e</sub> doesn't match Stacks' `Exp_Het`."**

> `Exp_Het` is the uncorrected 1 − Σp², and `Pi` applies the standard
> sample-size correction to it, so the two estimate the same parameter but
> differ by a factor of 2n/(2n−1) — 3.4% at n = 15. The correct comparison is against `Pi`, where the difference is under
> 0.5% and reflects only that our estimator does not assume F<sub>IS</sub> = 0.

**"Why not use Stacks' π directly?"**

> We do, for π. The correction underlying `Pi` assumes F<sub>IS</sub> = 0, which
> biases H<sub>e</sub> low by F/(2n−1) — under 0.5% here, negligible. But
> F<sub>IS</sub> divides by that same small quantity, so the identical error is
> amplified by roughly (1−F)/F, about ninefold, giving a 3–5% bias that does not
> cancel between samples of unequal size. Hence π from Stacks, F<sub>IS</sub>
> recomputed. `diversity_stats.R` prints both side by side.

**"Why report H<sub>e</sub> from the SNP dataset but F<sub>IS</sub> from the
haplotype dataset?"**

> Because they are affected differently by marker type. Haplotype gene diversity
> has no fixed ceiling — a tag with four SNPs can carry eight haplotypes, one
> with a single SNP can carry two — so it varies with read length, enzyme and
> SNP-calling parameters and is not comparable across studies. In simulation,
> varying SNP density alone moved haplotype H<sub>e</sub> by 48% and per-site
> H<sub>e</sub> by 10%. F<sub>IS</sub> is a ratio, so numerator and denominator
> inflate together and the scale cancels; there we prefer haplotypes because
> multi-allelic loci give more information per locus and narrower intervals. The
> two datasets give the same between-population ratio (1.96 vs 2.09), so the
> conclusion does not depend on the choice.

**"Why per polymorphic site rather than per sequenced site?"**

> Both are reported. Per-variant-site values are given for comparability with
> the existing literature on this system; per-sequenced-site (autosomal) values
> are given because they are the sample-size-unbiased quantity (Schmidt et al.
> 2021) and are what should be compared against whole-genome estimates.

**"Your H<sub>e</sub> is averaged over sites where the population is
monomorphic, which deflates it."**

> Deliberately, and both populations are scored on the same pooled denominator
> so the comparison is symmetric. Restricting each population to its own
> polymorphic sites conditions on the outcome: in simulation it erased a
> genuine 2.3-fold difference in diversity and marginally reversed its
> direction, because a population's monomorphic sites *are* its reduced
> diversity. The percentage of polymorphic sites is reported alongside, and
> H<sub>e</sub>(pooled) = H<sub>e</sub>(within-population) × % polymorphic
> exactly, so either version can be recovered from the table.

**"How much of your data actually went into these numbers?"**

> In the default (available-data) mode: a locus was retained for a population
> once that population had at least `min_n` (default 2) individuals genotyped
> there, decided independently per population and per locus — the run prints,
> per population, how many of N total loci cleared that bar. Per-sequenced-site
> values are scaled using the full variant-record count regardless, since the
> used subset is a sample of all variant records rather than the whole of them.
>
> If `--complete-case` was used instead: sites were retained only where every
> individual of every population was genotyped, following Schmidt et al. (2021)
> recommendation (b), because H<sub>o</sub> and H<sub>e</sub> diverge in
> proportion to the missingness allowed. N of M variant records met that
> criterion.

**"An allelic richness of 1.98 is uninformative."**

> Correct, and that is why allelic richness is reported from the haplotype
> dataset, where a RAD tag is one multi-allelic locus. On biallelic SNPs
> rarefied richness is bounded at 2 by construction.

**"Why not report the `Private` column from Stacks?"**

> That column is a raw count of private alleles and scales with sample size, so
> it cannot be compared between samples of 15 and 10. The values reported here
> are rarefied to a common number of gene copies (Kalinowski 2004; Szpiech
> et al. 2008).

**"Your H<sub>o</sub> is inconsistent with your π."**

> They come from different marker sets and are different quantities. Haplotype
> H<sub>o</sub> is the proportion of RAD tags at which an individual is
> heterozygous — true if any SNP on the tag is heterozygous — while π is
> per-site. The file each statistic comes from is given in the table.

**"Schmidt et al. (2021) recommend analysing populations in independent
`populations` runs to avoid SNP-ascertainment bias from joint calling across
differentiated populations. You called SNPs once, jointly. Why?"**

> Deliberately. Recommendations (a) and (b) — autosomal heterozygosity, and
> omitting sites with any missing data — are followed here; (c) is not, because
> it would break the thing a between-population comparison needs most: a shared
> marker set scored under identical rules. Separate runs can call a different
> variant set in each population, at which point Ho, He and F<sub>IS</sub> are
> no longer being computed over the same loci and are not directly comparable.
> The asymmetric-missingness problem that (c) guards against is instead handled
> here by scoring both populations under one joint, symmetric filter (`-r 0.8
> -p 2`, never `-R`, see above) — a pooled call set with a filter that treats
> both populations identically, rather than two independent call sets that
> might not agree on which sites are even variant.

**"You should have filtered for Hardy–Weinberg equilibrium."**

> A heterozygote deficit is the quantity under study, so an HWE filter removes
> the signal rather than noise. Applied to this dataset it cut F<sub>IS</sub>
> from 0.0721 to 0.0551, a 24% reduction created entirely by the filter. See
> also Pearman et al. (2022) on how HWE filtering schemes affect downstream
> inference.

**"Why a t-test on individuals rather than the standard bootstrap over loci?"**

> Because loci are repeated measurements on the same individuals, not
> independent replicates of the population mean. With two populations of
> identical true heterozygosity, a locus bootstrap rejected 55–76% of the time
> against a nominal 5% as soon as individuals varied in heterozygosity — which
> is what inbreeding produces. Welch's t on per-individual heterozygosity held
> its nominal rate at every level tested. The overdispersion factor reported in
> the results quantifies how badly the locus bootstrap would have failed on
> *these* data. The objection is not new (Van Dongen 1995).

**"Non-significant just means you were underpowered."**

> Correct, and we say so. Simulation at n = 15 versus 10 with realistic
> individual variation gave roughly 52% power for a moderate difference, so this
> is reported as no difference detected, not as evidence of no difference.

**"RAD data have null alleles; your F<sub>IS</sub> is inflated."**

> Agreed, which is why F<sub>IS</sub> is reported as an upper bound.
> Restriction-site null alleles inflate apparent homozygosity in both
> populations, so a shared upward bias cannot generate a difference between
> them; since the comparison is the claim, an upper bound supports it. If a
> difference had been detected, differential null load between populations would
> need ruling out.

## Deliberately not here

| left out | why |
|---|---|
| Null-allele / paralog screening | diagnosis, not analysis. Neither moves these numbers; the "One sentence on null alleles" above already covers what to write instead of measuring it directly. |
| HWE filtering | a heterozygote deficit is the signal, not the noise. Filtering on HWE deleted 24% of F<sub>IS</sub> in testing — created entirely by the filter. |
| Rarefaction, private allelic richness | `diversity_stats.R` reports them anyway; they are not part of the F<sub>IS</sub> question. |
| N<sub>e</sub> estimation | a different question on a different timescale. |
| A separate one-SNP-per-locus run | only needed for STRUCTURE, DAPC, outlier scans and N<sub>e</sub>. Not for this. |

`README.md` is the quick-start entry point; this document is the part that
produces the numbers in the results section and the place to go if a reviewer
challenges something.

---

## References

Bhatia, G. et al. (2013) Estimating and interpreting F<sub>ST</sub>. *Genome
Research* 23: 1514–1521.

Nei, M. & Chesser, R.K. (1983) Estimation of fixation indices and gene
diversities. *Annals of Human Genetics* 47: 253–259.

Neff, B.D. (2010) *Molecular Ecology Resources* 10: 546–550.

Pearman, W.S., Urban, L. & Alexander, A. (2022) Commonly used Hardy–Weinberg
equilibrium filtering schemes impact population structure inferences using
RADseq data. *Molecular Ecology Resources* 22: 2599–2613.

Petit, R.J. & Pons, O. (1998) Bootstrap variance of diversity and
differentiation estimators. *Heredity* 80: 56–61.

Hohenlohe, P.A., Bassham, S., Etter, P.D., Stiffler, N., Johnson, E.A. &
Cresko, W.A. (2010) Population genomics of parallel adaptation in threespine
stickleback using sequenced RAD tags. *PLoS Genetics* 6: e1000862.

Hall, P. (1992) *The Bootstrap and Edgeworth Expansion*. Springer.

Kalinowski, S.T. (2004) Counting alleles with rarefaction: private alleles and
hierarchical sampling designs. *Conservation Genetics* 5: 539–543.

Keenan, K., McGinnity, P., Cross, T.F., Crozier, W.W. & Prodöhl, P.A. (2013)
diveRsity: An R package for the estimation and exploration of population
genetics parameters and their associated errors. *Methods in Ecology and
Evolution* 4: 782–788.

Mammen, E. (1992) *When Does Bootstrap Work? Asymptotic Results and
Simulations*. Springer.

McCullagh, P. (2000) Resampling and exchangeable arrays. *Bernoulli* 6:
285–301.

Owen, A.B. & Eckles, D. (2012) Bootstrapping data arrays of arbitrary order.
*Annals of Applied Statistics* 6: 895–927. (The crossed-factor "product
weight" bootstrap cited for `--boot=both`; see "Bootstrap mode" in Step 2
for what its guarantee does and does not cover for this project's
estimators, checked by simulation rather than assumed.)

Szpiech, Z.A., Jakobsson, M. & Rosenberg, N.A. (2008) ADZE: a rarefaction
approach for counting alleles private to combinations of populations.
*Bioinformatics* 24: 2498–2504.

Rochette, N.C., Rivera-Colón, A.G. & Catchen, J.M. (2019) Stacks 2. *Molecular
Ecology* 28: 4737–4754.

Schmidt, T.L., Jasper, M.-E., Weeks, A.R. & Hoffmann, A.A. (2021) Unbiased
population heterozygosity estimates from genome-wide sequence data. *Methods in
Ecology and Evolution* 12: 1888–1898.

Van Dongen, S. (1995) How should we bootstrap allozyme data? *Heredity* 74:
445–447.

Waples, R.S. (2015) Testing for Hardy–Weinberg proportions: have we lost the
plot? *Journal of Heredity* 106: 1–19.

Weir, B.S. (1996) *Genetic Data Analysis II: Methods for Discrete
Population Genetic Data*. Sinauer Associates. (Delete-one jackknife over
loci for F<sub>ST</sub>/F<sub>IS</sub> standard errors, cited for the
`_se` columns' jackknife method.)

Weir, B.S. & Cockerham, C.C. (1984) Estimating F-statistics for the analysis of
population structure. *Evolution* 38: 1358–1370.
