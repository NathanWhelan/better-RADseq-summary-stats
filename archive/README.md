# Archive

Superseded or abandoned material, kept for history but not part of the
maintained pipeline. Nothing here is sourced or invoked by the current
scripts (`diversity_stats.R`, `diversity_core.R`, `haps_common.R`,
`het_between_pops.R`).

## `simpler-script/`

`fis_and_diversity_robust.R` is the pre-refactor draft that
`diversity_stats.R`/`diversity_core.R`/`haps_common.R` superseded. It works
on genepop files (via `adegenet`/`hierfstat`) rather than VCF, has its own
independent genotype-decoding and rarefaction code, and is hardcoded to one
dataset (`~/snails/L_virgata/...`, `popmap_virgata.txt`). It found a real bug
in an earlier version of this analysis (invalid locus-as-replicate
statistical tests) that motivated the current, tested implementation. Kept
for reference; do not extend it — any fix belongs in the current scripts.

## `gt6_comparison/`

A one-off comparison of `main_pipeline` (the current scripts) against
`simpler_script` (the draft above) on the same `gt6` dataset, run to check
the two implementations agreed before retiring the draft.

**Known issue with this comparison:** `main_pipeline/` was produced via
`results/run_virgata.sh` with the correct rarefaction target (`g = 28`,
derived from `-r 0.8` and the smallest population's raw size of 17), while
`simpler_script/fis_and_diversity_robust.R` independently computed its own
target as `2 * min(pop_sizes) = 34`. So the allelic-richness figures
(`Ar`/`privAr`) in this comparison were computed at two different
rarefaction depths and are not a clean apples-to-apples check — any
richness-value discrepancy here partly reflects that mismatch, not
necessarily a disagreement between the two implementations. Fis is
unaffected (it doesn't depend on `g`). Treat this directory as historical
context only, not as a validated cross-check.

## `virgata_m3M3n2gt01_r80p11_maf025_gt6/`

An abandoned/exploratory run using a different upstream `populations`
parameterization (`m3M3n2gt01`, vs. the `m3M3n3` used by every other
`results/virgata_*` dataset) and missing the SNP-level outputs and
`het_out/` directory its siblings have. Includes a stray chat-log file
(`chat_log_2026-09-05_diversity_discussion.md`) that was left inside it.
