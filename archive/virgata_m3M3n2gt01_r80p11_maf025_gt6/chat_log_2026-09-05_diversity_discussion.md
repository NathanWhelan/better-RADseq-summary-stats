# Discussion log: diversity stats interpretation
Folder: `/home/whelan.105/RADseq-Claude/results/virgata_m3M3n2gt01_r80p11_maf025_gt6`
Date: 2026-09-05

---

## Q: What is Hs in the context of the results in this folder?

In this pipeline's output, **Hs** is the within-population expected heterozygosity (gene diversity) as defined by Nei & Chesser (1983) — it's `hierfstat`'s name for the same quantity this results folder reports as **He**.

Concretely, per the log (`diversity_stats_haps.log`, line 29):
> `He = Nei & Chesser (1983) Hs, the hierfstat quantity`

So:
- Per-population, in `diversity_per_population.haps.tsv`, the `He` column *is* Hs for that population — the small-sample-corrected expected heterozygosity, averaged across all 1,559 haplotype loci (including ones monomorphic in that population, per the "pooled denominator" convention explained near the end of the log).
- In the overall `hierfstat::basic.stats()` summary line at the bottom of the log:
  ```
  Ho      Hs      Ht     Dst     Htp    Dstp     Fst    Fstp     Fis    Dest 
  0.1330  0.1116  0.2152  0.1036  0.2255  0.1139  0.4814  0.5052 -0.1919  0.1283 
  ```
  `Hs = 0.1116` is the **average within-subpopulation gene diversity** across all 11 populations (as opposed to `Ht`, the total/pooled gene diversity across the whole metapopulation). It's the standard hierarchical-Fstat decomposition: `Fst = Dst/Ht` where `Dst = Ht - Hs`.

And `FIS` in this run is computed directly from Hs and Ho as a ratio of sums: `FIS = 1 - sum(Ho)/sum(Hs)` (log line 30), not as a mean of per-locus ratios.

One caveat flagged explicitly in the log's "TAKE FROM THIS RUN" section: since this is the **haplotype** VCF, Ho/He(Hs) values aren't comparable across studies (they depend on RAD read length/enzyme/filters) — only `Fis`, `Ar`, and `privAr` are portable; for a "real" Hs/He worth comparing across studies you'd want the `populations.snps.vcf` per-site run instead.

---

## Q: What is Ht?

Following directly from that same overall `hierfstat::basic.stats()` line:

```
Ho      Hs      Ht     Dst     Htp    Dstp     Fst    Fstp     Fis    Dest 
0.1330  0.1116  0.2152  0.1036  0.2255  0.1139  0.4814  0.5052 -0.1919  0.1283 
```

**Ht** is the **total gene diversity** — the expected heterozygosity computed as if all 11 populations were pooled into one big random-mating population, using the average allele frequencies across the whole sample. It's `0.2152` here, versus `Hs = 0.1116` (the average diversity *within* populations).

The gap between the two is exactly what population structure produces:
- **Dst = Ht − Hs = 0.1036** — the between-population component of diversity.
- **Fst = Dst / Ht = 0.4814** — the fraction of total gene diversity attributable to differences among populations.

So the relationship is Hs ≤ Ht always, and how much smaller Hs is than Ht tells you how differentiated the populations are. Here Fst ≈ 0.48 is quite high — the populations are strongly differentiated relative to their within-population diversity.

The `p` variants (`Htp`, `Dstp`, `Fstp`) are the same quantities with Nei & Chesser's small-sample bias correction applied (matches the "p" = "prime" convention in the Fstat/hierfstat literature) — `Fstp = 0.5052` is the corrected estimate you'd typically report over the uncorrected `Fst`.

---

## Q: That Fst is huge. But it's the same species. What's going on?

Two things are going on—one mathematical, one biological—and the log already hands you the tool to tell them apart.

**1. Fst is mathematically inflated when Hs is low (the bigger factor here)**

Fst = Dst/Ht is bounded above by `1 − Hs`. With `Hs = 0.1116` (fairly low, typical for haplotype-based RAD loci), the *maximum possible* Fst for this dataset is ~0.89 — so even modest actual allele-frequency divergence among populations gets stretched into a large standardized Fst. This is the well-known Hedrick (2005)/Jost (2008) critique of Fst as a measure of differentiation.

The log gives you the fix already, in the same summary line:
```
Ho      Hs      Ht     Dst     Htp    Dstp     Fst    Fstp     Fis    Dest 
0.1330  0.1116  0.2152  0.1036  0.2255  0.1139  0.4814  0.5052 -0.1919  0.1283 
```
**Dest = 0.1283** is Jost's D — it's not normalized against Hs the way Fst is, so it isn't inflated by low within-population heterozygosity. It's saying the actual degree of allelic differentiation is moderate (~0.13), not the ~0.48–0.51 that Fst/Fstp report. When Fst and Dest disagree this much, Dest is the more honest number, and Fst ≈ 0.48 alone would be a misleading headline stat.

**2. The population labels are separate stream drainages**

LookoutCreek, FlintRiver, Sequatchie, PaintRock, Dry_Creek, Pidgeon_River, BentCreek, SpringCreek — these read like distinct creek/river systems, not sites within one drainage. For a low-dispersal, stream-restricted species, genuine near-zero gene flow *between drainages* (versus free mixing within one) is completely ordinary and can produce real Fst in the 0.3–0.6 range even intraspecifically — this is a very standard pattern for headwater stream endemics (salamanders, crayfish, freshwater snails, small-bodied fish). So a real, biologically-driven high value isn't implausible for "the same species" if dispersal is drainage-limited.

**One thing worth double-checking:** every single population has a fairly strongly *negative* Fis (−0.06 to −0.33 — consistent heterozygote excess everywhere). That's unusual enough to be worth a sanity check before trusting the Fst/Dest numbers at face value — it can be genuine (small local Ne, avoidance of inbreeding) or it can be a Stacks/RAD artifact (e.g., over-merged paralogous loci inflating apparent heterozygosity uniformly). Since it's the same sign and rough magnitude in all 11 populations, I'd lean toward it being either a real shared biological pattern or a systematic calling/assembly effect rather than population-specific noise — worth checking against the SNP-based `populations.snps.vcf` run (which the log recommends anyway for comparable Ho/He) to see if the negative Fis persists there too.

**Bottom line:** report Dest (~0.13) as your primary differentiation statistic, cite Fst/Fstp (~0.48–0.51) only alongside it with the caveat that it's amplified by low Hs, and treat the drainage-level structure as plausibly real for a dispersal-limited stream species — but verify the universal negative Fis isn't a paralog artifact before leaning too hard on any of these numbers.

---

## Q: I really don't think these are paralogs, but it is odd. What could be going on with such high Hobs vs Hexp (and negative Fis)? Can you find other examples in low-dispersing freshwater taxa? Find citations that would support your conclusions.

*(Context established via project directory: sample sites correspond to `Leptoxis virgata`, a Pleuroceridae freshwater snail endemic to the Tennessee River drainage in Alabama/Tennessee, per `/fs/project/whelan.105/virgata/README`.)*

Good instinct to rule out paralogs — nothing here needs that. Two mechanisms in the literature produce the same signature (Ho > Hs, negative Fis), and both fit dam-fragmented Pleuroceridae specifically well. Then some close taxonomic comparators showing your Fst isn't an outlier for this group.

### Why heterozygote excess / negative FIS, without invoking paralogs

**1. Transient heterozygosity excess after a bottleneck (Cornuet & Luikart 1996).**
When Ne drops sharply, rare alleles are purged fast but heterozygosity barely moves — losing a rare allele removes very little expected heterozygosity, so allelic richness collapses faster than gene diversity does. For several generations afterward, observed/realized heterozygosity sits *above* what mutation–drift equilibrium (as inferred from the reduced allele count) would predict — the classic "heterozygosity excess" bottleneck signature, which shows up as negative Fis in exactly this kind of dataset.

This fits your system unusually well: *L. virgata* is documented as having been "drastically diminished" by 20th-century Tennessee River impoundment (Guntersville and related dams), and its congeners show the same history. Wright et al. (2020, PeerJ) explicitly test and discuss a bottleneck hypothesis for *L. compacta* after 20th-century range collapse in the same drainage system — same demographic story, same taxon.

- Cornuet JM, Luikart G (1996). Description and Power Analysis of Two Tests for Detecting Recent Population Bottlenecks From Allele Frequency Data. *Genetics* 144(4):2001–2014. https://doi.org/10.1093/genetics/144.4.2001
- Wright AD, Garrison NL, Williams AS, Johnson PD, Whelan NV (2020). Range reduction of Oblong Rocksnail, *Leptoxis compacta*, shapes riverscape genetic patterns. *PeerJ* 8:e9789. https://doi.org/10.7717/peerj.9789

**2. Small effective number of breeders / family-correlated sampling (the "Hedgecock effect").**
In organisms with high fecundity and strongly skewed reproductive success (a few individuals contributing disproportionately to a cohort), a sample dominated by full/half-sibs shows heterozygote excess relative to the population's allele frequencies, because sib groups are less variable than random mating would predict. This is well documented in oysters and other aquatic invertebrates with sweepstakes-style recruitment, and it's directly testable: if your creek samples skew toward a narrow size/age class (recent recruits from one spawning event) rather than a mixed-age assemblage, this mechanism alone can generate the pattern in every site independently — no need for it to be the same cause everywhere, just the same demographic tendency (patchy recruitment, few effective spawners per site).

- Pudovkin AI, Zaykin DV, Hedgecock D (1996). On the Potential for Estimating the Effective Number of Breeders from Heterozygote-Excess in Progeny. *Genetics* 144(1):383–387. https://pmc.ncbi.nlm.nih.gov/articles/PMC1207510/

Both mechanisms predict the *same* sign of effect you see (universal, not locus-specific, roughly consistent in magnitude across your 11 sites) — which is a point in their favor over something locus-specific like paralogy.

### Is Fst ~0.48–0.51 actually anomalous for this genus? No.

Congeneric/confamilial pleurocerid snails in the *same* fragmented Alabama/Tennessee drainage systems, studied by essentially the same lineage of work, routinely report Fst in this range or higher:

- **Whelan NV, Galaska MP, Sipley BN, Weber JM, Johnson PD, Halanych KM, Helms BS (2019).** Riverscape genetic variation, migration patterns, and morphological variation of the threatened Round Rocksnail, *Leptoxis ampla*. *Molecular Ecology* — pairwise Fst **0.377–0.773** among Cahaba River sites.
- **Redak CA, Williams AS, Garner JT, Halanych KM, Whelan NV (2021).** Assessing genomic diversity, connectivity, and riverscape genetics hypotheses in the endangered rough hornsnail, *Pleurocera foremani*, following habitat disruption. *Journal of Heredity* 112(7):635–645. https://doi.org/10.1093/jhered/esab065 — Fst up to **0.719** between Yellowleaf Creek and Lake Mitchell sites.
- **Minton RL et al. (2017).** Genetic structuring in the Pyramid Elimia, *Elimia potosiensis* (Gastropoda, Pleuroceridae), with implications for pleurocerid conservation. *Zoosystematics and Evolution* — significant structure even over a few hundred meters. https://zse.pensoft.net/article/14856/

Pleurocerids are direct-developing, crawl-only, no planktonic stage — about as dispersal-limited as a freshwater animal gets — so drift-driven Fst in the 0.4–0.7+ range between drainages/reaches is essentially the expected outcome for this family, not a red flag.

One useful contrast: Wright et al. (2020) found the *opposite* sign in *L. compacta* — heterozygote **deficit** (positive Fis) — which they attribute to a Wahlund effect from cryptic subpopulation mixing within sites. That your data go the other direction (consistent excess) argues against simple within-site lumping of divergent subpopulations as the explanation, and points back toward bottleneck/family-structure mechanisms instead.

**A concrete way to adjudicate between the two mechanisms in your own data:** check whether the negative-Fis magnitude tracks site isolation/small effective size (favors bottleneck) versus whether it tracks apparent size/age-class skew in your specimen collections (favors the few-breeders/sib-sample effect) — the two aren't mutually exclusive and could well both be operating given the dam history.

Sources: [Cornuet & Luikart 1996](https://academic.oup.com/genetics/article-abstract/144/4/2001/6017142) · [Pudovkin, Zaykin & Hedgecock 1996](https://pmc.ncbi.nlm.nih.gov/articles/PMC1207510/) · [Wright et al. 2020, PeerJ](https://pmc.ncbi.nlm.nih.gov/articles/PMC7473045/) · [Redak et al. 2021, J Hered](https://academic.oup.com/jhered/article/112/7/635/6407498) · [Whelan et al. 2019, Leptoxis ampla](https://www.researchgate.net/publication/330730129_Riverscape_genetic_variation_migration_patterns_and_morphological_variation_of_the_threatened_Round_Rocksnail_Leptoxis_ampla) · [Minton et al. 2017, Elimia potosiensis](https://zse.pensoft.net/article/14856/)

---

## Q: But these populations are locally abundant. There are thousands or millions at any given site.

Fair pushback, and it's the right question — a bottleneck story alone doesn't sit well with "thousands to millions of individuals per site right now." The resolution is that **census size (Nc) and genetic effective population size (Ne) are not the same thing**, and for high-fecundity aquatic invertebrates they can differ by orders of magnitude — which explains both the negative Fis and the high Fst *without* needing the population to actually be small or recently crashed.

### The core idea: Ne << Nc in high-fecundity aquatic invertebrates

A freshwater snail can lay huge numbers of eggs, but if reproductive success is highly skewed — a small subset of females/egg masses contribute disproportionately to the surviving cohort each generation, while most offspring die before reproducing — then the *genetic* population size each generation is set by the handful of successful parents, not by the huge number of adults crawling around. This is "sweepstakes reproductive success," and it routinely produces Ne/Nc ratios far below what you'd guess from headcounts.

- Hedgecock D (1994). Does variance in reproductive success limit effective population sizes of marine organisms? In: *Genetic and Evolutionary Aspects of Marine Organisms*, ed. Beaumont. — the original sweepstakes hypothesis: effective size can be orders of magnitude below census size.
- Hedgecock D, Pudovkin AI (2011). Sweepstakes Reproductive Success in Highly Fecund Marine Fish and Shellfish: A Review and Commentary. *Bulletin of Marine Science* 87(4):971–1002.
- Palstra FP, Fraser DJ (2012). Effective/census population size ratio estimation: a compendium and appraisal. *Ecology and Evolution* 2(9):2357–2365. https://doi.org/10.1002/ece3.329 — median Ne/Nc for invertebrates specifically is **~0.09**, i.e., roughly a tenth or less of census size is typical, and it can go much lower under sweepstakes recruitment.

### This directly produces both signatures you're seeing, every generation, independent of current abundance

**Negative Fis** falls out mechanistically: if a site's snails sampled in one collecting trip are disproportionately drawn from a small number of successful egg-laying events, they're a set of large sibships rather than a random draw from millions of unrelated adults — and a sample dominated by sibships shows heterozygote excess relative to population allele frequencies (Pudovkin, Zaykin & Hedgecock 1996, cited earlier). This doesn't require the population to be small or bottlenecked at any point — it recurs every generation as long as recruitment is skewed.

**High Fst** falls out the same way: if each site's gene pool is effectively re-founded each generation by a small number of successful spawners, per-generation drift is fast even though the standing adult census is enormous — because drift depends on Ne, not Nc. This is exactly the phenomenon named **"chaotic genetic patchiness"** — strong, spatially/temporally unstable genetic differentiation among sites that are demonstrably abundant and even physically connected, first described in an intertidal limpet and since documented repeatedly in molluscs:

- Johnson MS, Black R (1982). Chaotic genetic patchiness in an intertidal limpet, *Siphonaria* sp. *Marine Biology* 70:157–164.
- Riquet F, Le Cam S, Fonteneau E, Viard F (2016). Moderate genetic drift is driven by extreme recruitment events in the invasive mollusk *Crepidula fornicata*. *Heredity* 117:42–50. https://doi.org/10.1038/hdy.2016.24 — a direct gastropod case: abundant, well-connected populations still showing drift-driven genetic patchiness traced to skewed recruitment events.

### Why this fits your data better than a pure historical-bottleneck story

A past dam-era bottleneck predicts a genetic signature that should be *decaying* as census recovers — allelic richness eventually catches back up via mutation/migration, even if that takes a long time. Sweepstakes/chaotic-patchiness instead predicts the signature to be a **standing, recurring feature** of the system regardless of current census — which matches "abundant now, still negative Fis, still high Fst" much more comfortably than "still recovering from a 90-year-old bottleneck" does. The two aren't mutually exclusive (a historical crash could have set the current drainage-level divergence, while ongoing sweepstakes recruitment keeps regenerating heterozygote excess within sites every generation), but given current abundance, sweepstakes/chaotic patchiness should be treated as the primary explanation, with the dam history as a secondary contributor to among-site divergence rather than the driver of the within-site Fis pattern.

**A testable prediction this gives you:** if you have any temporal replication (same site, different years/cohorts) or can bin samples by size/age class, sweepstakes recruitment predicts Fis and even fine-scale allele-frequency composition should shift between collecting events at the same site, whereas a pure legacy-bottleneck story predicts a stable signature over time.

Sources: [Hedgecock 1994 summary](https://www.researchgate.net/publication/245970620_Does_variance_in_reproductive_success_limit_effective_population_sizes_of_marine_organisms_In_A) · [Hedgecock & Pudovkin 2011](https://www.researchgate.net/publication/233542350_Sweepstakes_Reproductive_Success_in_Highly_Fecund_Marine_Fish_and_Shellfish_A_Review_and_Commentary) · [Palstra & Fraser 2012](https://onlinelibrary.wiley.com/doi/full/10.1002/ece3.329) · [Johnson & Black 1982](https://link.springer.com/article/10.1007/BF00397680) · [Riquet et al. 2016](https://pmc.ncbi.nlm.nih.gov/articles/PMC4901356)

---

## Q: What's the test in PLINK to test for sibship or family structure?

Two options in PLINK, and given this data (structured populations, negative Fis, non-HWE), one is clearly the better choice.

### Option A — classic PLINK1.9 IBD (`--genome`)

```bash
# LD-prune first — IBD estimation assumes independent markers
plink --bfile mydata --indep-pairwise 50 5 0.2 --out pruned

# pairwise IBD estimation
plink --bfile mydata --extract pruned.prune.in --genome --out ibd
```

Output is `ibd.genome`, with `Z0`, `Z1`, `Z2` (probability of sharing 0/1/2 alleles IBD) and `PI_HAT = Z1/2 + Z2`. Rough expectations:

| Relationship | PI_HAT |
|---|---|
| Unrelated | ~0 |
| Parent–offspring | ~0.5 (Z1≈1) |
| Full sibs | ~0.5 (Z0≈.25, Z1≈.5, Z2≈.25) |
| Half-sib / avuncular / grandparent | ~0.25 |
| First cousins | ~0.125 |

**Caveat that matters a lot here:** `--genome`'s method-of-moments estimator assumes the sample is drawn from a single randomly-mating population in HWE with accurately estimated allele frequencies. These 11 populations are exactly the opposite of that — differentiated demes each showing departures from HWE. Run on the pooled dataset, shared ancestry/structure gets misread as relatedness and inflates PI_HAT for everyone, which would make a sibship hypothesis look artificially true. If using this method, **run it separately within each population** (not pooled across all 11), on a population-specific pruned SNP set.

### Option B — KING-robust estimator (recommended for this case)

```bash
plink2 --bfile mydata --king-cutoff 0.0884 --out king
# or, to get the full pairwise table instead of just a filtered/pruned set:
plink2 --bfile mydata --make-king-table --king-table-filter 0.0884 --out king
```

The KING-robust estimator (Manichaikul et al. 2010, *Bioinformatics* 26(22):2867–2873) is specifically designed to be **robust to unknown population structure and departures from HWE** — it doesn't need accurate population allele frequencies the way `--genome` does. That property is exactly what's wanted when the populations under test are the ones showing the anomalous Fis in the first place. Run it per-population.

Kinship-coefficient cutoffs (same paper): >0.354 duplicate/MZ twin, 0.177–0.354 first-degree (parent-offspring or full-sib — would need `Z0`/IBS patterns or a separate tool like KING's `--related` or `SNPRelate`/`sequoia` to split those two), 0.0884–0.177 second-degree (half-sib, avuncular, grandparent), 0.0442–0.0884 third-degree.

### What to actually check against the Fis result

Within each population, tabulate the pairwise kinship/PI_HAT distribution and see whether it's dominated by an excess of first/second-degree pairs relative to what random mating in a population that size would predict — i.e., is the sample effectively a small number of large sibships rather than ~20 unrelated individuals? If a site's 20 samples resolve into 3–4 family clusters, that directly supports the sweepstakes/few-effective-breeders explanation for the heterozygote excess at that site, and re-running the diversity stats after keeping one representative per sibship would show whether Fis moves back toward zero.
