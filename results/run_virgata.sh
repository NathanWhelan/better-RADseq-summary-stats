#!/bin/bash
set -uo pipefail

SRC="/fs/project/whelan.105/virgata/virgata-assembly/post-submission-filtering-tests"
POPMAP="/fs/project/whelan.105/virgata/virgata-assembly/popmap_virgata.txt"
REPO="/home/whelan.105/RADseq-Claude"
DIVSTATS="$REPO/diversity_stats.R"
HETPOPS="$REPO/het_between_pops.R"
## diversity_stats.R now requires --g explicitly (it used to default to twice
## the smallest population). Twice the smallest population's raw size (2*17
## = 34) overstates what's actually usable: `populations` was run with
## -r 0.8, so only ~80% of a population's individuals are typically
## genotyped at a given locus. g = 2 * ceil(r * n_min) = 2 * ceil(0.8 * 17)
## = 2 * 14 = 28, using r = 0.8 (the -r call-rate flag) and n_min = 17
## (Pidgeon_River, the smallest population in popmap_virgata.txt).
## Revisit if the popmap or the -r value changes.
GENE_COPIES=28

declare -A SITES=(
  [virgata_m3M3n3_r80p11_maf025]=225175
  [virgata_m3M3n3_r80p11_maf025_gt10]=222142
  [virgata_m3M3n3_r80p11_maf025_gt10_H]=222139
  [virgata_m3M3n3_r80p11_maf025_gt10_single]=222111
  [virgata_m3M3n3_r80p11_maf025_gt6]=222861
  [virgata_m3M3n3_r80p11_maf025_gt6_H]=222857
  [virgata_m3M3n3_r80p11_maf025_gt6_single]=222825
)

for ds in "${!SITES[@]}"; do
  #echo $ds
  #echo "TEST!"
  outdir="$REPO/results/$ds"
  sites="${SITES[$ds]}"
  echo "=== $ds (Sites=$sites) starting $(date) ==="
  mkdir $outdir
  ( cd "$outdir" && Rscript "$DIVSTATS" "$SRC/$ds/populations.haps.vcf" "$POPMAP" \
      --g="$GENE_COPIES" --nboot=10000 \
      > diversity_stats_haps.log 2>&1 )
  echo "  haps diversity_stats.R exit=$? $(date)"

  ( cd "$outdir" && Rscript "$DIVSTATS" "$SRC/$ds/populations.snps.vcf" "$POPMAP" \
      --g="$GENE_COPIES" --nboot=10000 --sites="$sites" \
      > diversity_stats_snps.log 2>&1 )
  echo "  snps diversity_stats.R exit=$? $(date)"

  ( cd "$outdir" && Rscript "$HETPOPS" "$SRC/$ds/populations.snps.vcf" "$POPMAP" 0.9 "het_out" \
      > het_between_pops.log 2>&1 )
  echo "  het_between_pops.R exit=$? $(date)"

  echo "=== $ds done $(date) ==="
done

echo "ALL DATASETS COMPLETE $(date)"
