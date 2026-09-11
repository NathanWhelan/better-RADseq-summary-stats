#!/usr/bin/env python3
"""
filter_low_conf_alt.py

Identify genotype calls where the ALT allele is called with very few
supporting reads (from the VCF AD field), and optionally filter them
in one of two ways:

  1. --mask-genotypes OUT.vcf
       Write a new VCF where ONLY the flagged genotypes are set to
       missing (./.). The locus itself, and every other individual's
       genotype at that locus, are left untouched. Use this when you
       want to keep a locus but can't trust one or two individuals'
       calls there.

  2. --drop-loci OUT.vcf
       Write a new VCF where entire records are removed if they meet
       a record-level criterion. By default a record is dropped if it
       contains ANY flagged genotype; raise --drop-loci-min-frac to
       only drop records where a larger fraction of called genotypes
       are flagged (useful for distinguishing "one bad individual"
       from "this whole locus looks unreliable").

Both options can be given in the same run (they produce two separate
output VCFs from the same flagging pass). Neither is required if you
only want the summary CSVs.

Flagging rule
-------------
A genotype call is flagged if:
  - it is fully called (no missing alleles), AND
  - it includes at least one ALT allele (e.g. 0/1, 1/1, 1/2), AND
  - the read depth supporting the ALT allele(s) present in that
    genotype (summed, for multi-allelic genotypes) is <= --min-alt-reads

--min-alt-reads defaults to 2. This is an absolute-count threshold, not
a depth-normalized one -- see --sensitivity-table, and note that "how
low is too low" is a judgment call, not a statistically derived cutoff.

Notes on VCF edge cases this script handles
--------------------------------------------
- Records sharing the same CHROM:POS (e.g. from split multiallelic
  sites) are tracked independently -- each physical VCF record gets
  its own stats and its own keep/drop decision. They are NOT merged.
- Per the VCF spec, trailing FORMAT subfields may be dropped entirely
  when missing (e.g. "0/1:20" for a FORMAT of GT:DP:AD:GQ). Such a
  genotype is still counted as an ALT-containing call; it just can't
  be flagged (no AD to check). This is different from a fully-missing
  genotype ("./." with no other subfields), which is not counted at all.
- If a data row has a different number of sample columns than the
  #CHROM header declares, this is a malformed VCF; the script warns
  (rather than silently misaligning or dropping data) and skips
  per-genotype processing for that row.
- The script never loads the whole VCF into memory. It streams the
  input file multiple times (once to compute stats, once per requested
  output) so peak memory stays roughly constant regardless of file size.

Usage
-----
    python filter_low_conf_alt.py input.vcf \\
        --min-alt-reads 2 \\
        --mask-genotypes masked.vcf \\
        --drop-loci dropped.vcf \\
        --drop-loci-min-frac 0.0 \\
        --summary-calls low_confidence_ALT_calls.csv \\
        --summary-loci low_alt_support_loci_summary.csv

All output arguments are optional -- only the outputs you request are written.
"""

import argparse
import csv
import sys


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("vcf", help="Input VCF file (must have GT, AD and DP in FORMAT)")
    p.add_argument("--min-alt-reads", type=int, default=2,
                    help="Flag genotype calls with ALT-supporting reads <= this value (default: 2)")
    p.add_argument("--mask-genotypes", metavar="OUT.vcf", default=None,
                    help="Write a VCF with flagged genotypes set to missing (./.), locus kept")
    p.add_argument("--drop-loci", metavar="OUT.vcf", default=None,
                    help="Write a VCF with whole records removed per --drop-loci-min-frac")
    p.add_argument("--drop-loci-min-frac", type=float, default=0.0,
                    help="Drop a record if the fraction of ALT-containing genotypes that are "
                         "flagged exceeds this value (default: 0.0, i.e. drop on ANY flagged "
                         "genotype in that record)")
    p.add_argument("--summary-calls", metavar="OUT.csv", default=None,
                    help="Write a CSV of every flagged individual call")
    p.add_argument("--summary-loci", metavar="OUT.csv", default=None,
                    help="Write a CSV summarizing flagged calls per record")
    p.add_argument("--sensitivity-table", action="store_true",
                    help="Print flagged-call counts across a range of thresholds and exit "
                         "(no other outputs written)")
    return p.parse_args()


def parse_genotype(gt_str):
    """Return list of allele ints, or None if genotype has any missing allele."""
    gt_str = gt_str.replace("|", "/")
    alleles = gt_str.split("/")
    if any(a == "." for a in alleles):
        return None
    return [int(a) for a in alleles]


def iter_header_and_data(path):
    """Yield ('header', line) for header lines, then ('data', line) for data lines."""
    with open(path) as f:
        for line in f:
            if line.startswith("##"):
                yield ("header", line)
            elif line.startswith("#CHROM"):
                yield ("header", line)
            else:
                yield ("data", line)


def get_samples(path):
    with open(path) as f:
        for line in f:
            if line.startswith("#CHROM"):
                return line.rstrip("\n").split("\t")[9:]
    return []


def evaluate_genotype(raw, gt_i, ad_i, dp_i, min_alt_reads):
    """
    Parse one sample's colon-delimited genotype field.
    Returns a dict:
      {"is_alt_call": bool, "flagged": bool, "ad_ref": int, "ad_alt": int,
       "dp": int, "alt_fraction": float, "gt": str}
    or None if the genotype is entirely missing / GT itself unavailable.
    """
    parts = raw.split(":")
    if len(parts) <= gt_i:
        return None  # can't even recover GT -- e.g. bare "."

    gt_str = parts[gt_i]
    alleles = parse_genotype(gt_str)
    if alleles is None or max(alleles) == 0:
        return None  # missing or homozygous reference -- not an ALT call

    result = {"is_alt_call": True, "flagged": False, "gt": gt_str,
              "ad_ref": None, "ad_alt": None, "dp": None, "alt_fraction": float("nan")}

    # AD may be legitimately absent (trailing FORMAT fields dropped per VCF spec)
    if len(parts) <= ad_i:
        return result
    ad_field = parts[ad_i]
    if ad_field in (".", ""):
        return result
    ad = [int(x) if x != "." else 0 for x in ad_field.split(",")]

    if len(parts) <= dp_i or parts[dp_i] in (".", ""):
        dp = sum(ad)
    else:
        dp = int(parts[dp_i])

    alt_allele_idxs = sorted(set(a for a in alleles if a >= 1))
    ad_alt = sum(ad[i] for i in alt_allele_idxs if i < len(ad))
    ad_ref = ad[0] if ad else 0

    result.update({
        "ad_ref": ad_ref, "ad_alt": ad_alt, "dp": dp,
        "alt_fraction": (ad_alt / dp) if dp else float("nan"),
        "flagged": ad_alt <= min_alt_reads,
    })
    return result


def scan_vcf(path, min_alt_reads):
    """
    Single streaming pass. Returns:
      samples       : list of sample names
      flagged_calls : list of dicts, one per flagged genotype call (includes CHROM/POS/record_idx)
      records       : list of dicts, one per VCF data record IN FILE ORDER:
                       {chrom, pos, n_alt_calls, n_flagged}
                       (records sharing the same CHROM:POS are kept as separate entries)
      n_mismatched_rows : count of data rows whose sample-column count didn't match the header
    """
    samples = get_samples(path)
    n_samples_expected = len(samples)

    flagged_calls = []
    records = []
    n_mismatched_rows = 0
    warned_examples = []

    for kind, line in iter_header_and_data(path):
        if kind == "header":
            continue
        fields = line.rstrip("\n").split("\t")
        chrom, pos = fields[0], fields[1]
        sample_fields = fields[9:]

        if len(sample_fields) != n_samples_expected:
            n_mismatched_rows += 1
            if len(warned_examples) < 5:
                warned_examples.append(
                    f"{chrom}:{pos} has {len(sample_fields)} sample columns, "
                    f"expected {n_samples_expected}"
                )
            records.append({"chrom": chrom, "pos": pos, "n_alt_calls": 0, "n_flagged": 0})
            continue

        fmt = fields[8].split(":")
        try:
            gt_i, ad_i, dp_i = fmt.index("GT"), fmt.index("AD"), fmt.index("DP")
        except ValueError:
            records.append({"chrom": chrom, "pos": pos, "n_alt_calls": 0, "n_flagged": 0})
            continue

        n_alt_calls = 0
        n_flagged = 0
        record_idx = len(records)
        for sample, raw in zip(samples, sample_fields):
            res = evaluate_genotype(raw, gt_i, ad_i, dp_i, min_alt_reads)
            if res is None:
                continue
            n_alt_calls += 1
            if res["flagged"]:
                n_flagged += 1
                flagged_calls.append({
                    "CHROM": chrom, "POS": pos, "record_idx": record_idx, "sample": sample,
                    "GT": res["gt"], "AD_ref": res["ad_ref"], "AD_alt": res["ad_alt"],
                    "DP": res["dp"], "alt_fraction": res["alt_fraction"],
                })

        records.append({"chrom": chrom, "pos": pos,
                         "n_alt_calls": n_alt_calls, "n_flagged": n_flagged})

    if n_mismatched_rows:
        print(f"WARNING: {n_mismatched_rows} data row(s) had a different number of sample "
              f"columns than the header declares ({n_samples_expected}). These rows were "
              f"NOT processed for genotype flagging. Examples:", file=sys.stderr)
        for ex in warned_examples:
            print(f"    {ex}", file=sys.stderr)

    return samples, flagged_calls, records, n_mismatched_rows


def stream_write_masked(in_path, out_path, samples, min_alt_reads):
    """Stream through the input again, writing flagged genotypes as './.', without
    holding the whole file in memory."""
    n_samples_expected = len(samples)
    n_masked = 0
    with open(out_path, "w") as out:
        for kind, line in iter_header_and_data(in_path):
            if kind == "header":
                out.write(line)
                continue
            fields = line.rstrip("\n").split("\t")
            sample_fields = fields[9:]
            if len(sample_fields) != n_samples_expected:
                out.write(line)  # can't safely process; pass through unchanged
                continue
            fmt = fields[8].split(":")
            try:
                gt_i, ad_i, dp_i = fmt.index("GT"), fmt.index("AD"), fmt.index("DP")
            except ValueError:
                out.write(line)
                continue

            new_sample_fields = []
            for raw in sample_fields:
                res = evaluate_genotype(raw, gt_i, ad_i, dp_i, min_alt_reads)
                if res is not None and res["flagged"]:
                    n_fields = len(raw.split(":"))
                    masked = ["."] * n_fields
                    masked[gt_i] = "./."
                    new_sample_fields.append(":".join(masked))
                    n_masked += 1
                else:
                    new_sample_fields.append(raw)
            out.write("\t".join(fields[:9] + new_sample_fields) + "\n")
    return n_masked


def stream_write_dropped(in_path, out_path, records, min_frac):
    """Stream through the input again, writing only records that pass the
    drop-loci criterion. `records` (from scan_vcf, in file order) tells us,
    per physical record, whether to keep it -- avoiding any CHROM:POS collision."""
    keep_flags = []
    for r in records:
        if r["n_alt_calls"] == 0:
            keep_flags.append(True)
            continue
        frac = r["n_flagged"] / r["n_alt_calls"]
        if min_frac == 0.0:
            keep_flags.append(r["n_flagged"] == 0)
        else:
            keep_flags.append(frac <= min_frac)

    n_dropped = sum(1 for k in keep_flags if not k)
    idx = 0
    with open(out_path, "w") as out:
        for kind, line in iter_header_and_data(in_path):
            if kind == "header":
                out.write(line)
                continue
            if keep_flags[idx]:
                out.write(line)
            idx += 1
    return n_dropped


def print_sensitivity_table(path):
    samples = get_samples(path)
    n_samples_expected = len(samples)
    ad_alt_values = []
    for kind, line in iter_header_and_data(path):
        if kind == "header":
            continue
        fields = line.rstrip("\n").split("\t")
        sample_fields = fields[9:]
        if len(sample_fields) != n_samples_expected:
            continue
        fmt = fields[8].split(":")
        try:
            gt_i, ad_i, dp_i = fmt.index("GT"), fmt.index("AD"), fmt.index("DP")
        except ValueError:
            continue
        for raw in sample_fields:
            res = evaluate_genotype(raw, gt_i, ad_i, dp_i, min_alt_reads=10**9)
            if res is not None and res["ad_alt"] is not None:
                ad_alt_values.append(res["ad_alt"])

    total = len(ad_alt_values)
    print(f"Total ALT-containing genotype calls with usable AD: {total}")
    for t in [1, 2, 3, 4, 5, 10]:
        n = sum(1 for aa in ad_alt_values if aa <= t)
        print(f"  AD_alt <= {t:>2}: {n:>7} calls ({100*n/total:.2f}%)" if total else "  (no data)")


def main():
    args = parse_args()

    if args.sensitivity_table:
        print_sensitivity_table(args.vcf)
        return

    samples, flagged_calls, records, n_mismatched = scan_vcf(args.vcf, args.min_alt_reads)

    total_alt_calls = sum(r["n_alt_calls"] for r in records)
    n_records_affected = sum(1 for r in records if r["n_flagged"] > 0)
    print(f"Samples: {len(samples)}", file=sys.stderr)
    print(f"Records: {len(records)}", file=sys.stderr)
    print(f"ALT-containing genotype calls: {total_alt_calls}", file=sys.stderr)
    print(f"Flagged calls (AD_alt <= {args.min_alt_reads}): {len(flagged_calls)}", file=sys.stderr)
    print(f"Records with >=1 flagged call: {n_records_affected}", file=sys.stderr)

    if args.summary_calls:
        with open(args.summary_calls, "w", newline="") as fh:
            writer = csv.DictWriter(fh, fieldnames=[
                "CHROM", "POS", "sample", "GT", "AD_ref", "AD_alt", "DP", "alt_fraction"
            ])
            writer.writeheader()
            for row in flagged_calls:
                writer.writerow({k: row[k] for k in writer.fieldnames})
        print(f"Wrote {args.summary_calls}", file=sys.stderr)

    if args.summary_loci:
        with open(args.summary_loci, "w", newline="") as fh:
            writer = csv.DictWriter(fh, fieldnames=[
                "CHROM", "POS", "n_alt_calls", "n_flagged", "flagged_fraction"
            ])
            writer.writeheader()
            for r in sorted(records, key=lambda r: -r["n_flagged"]):
                if r["n_flagged"] == 0:
                    continue
                frac = r["n_flagged"] / r["n_alt_calls"] if r["n_alt_calls"] else 0
                writer.writerow({
                    "CHROM": r["chrom"], "POS": r["pos"],
                    "n_alt_calls": r["n_alt_calls"], "n_flagged": r["n_flagged"],
                    "flagged_fraction": round(frac, 4),
                })
        print(f"Wrote {args.summary_loci}", file=sys.stderr)

    if args.mask_genotypes:
        n_masked = stream_write_masked(args.vcf, args.mask_genotypes, samples, args.min_alt_reads)
        print(f"Wrote {args.mask_genotypes} ({n_masked} genotypes masked to ./.)", file=sys.stderr)

    if args.drop_loci:
        n_dropped = stream_write_dropped(args.vcf, args.drop_loci, records, args.drop_loci_min_frac)
        print(f"Wrote {args.drop_loci} ({n_dropped} records dropped, "
              f"threshold: flagged_fraction > {args.drop_loci_min_frac})", file=sys.stderr)


if __name__ == "__main__":
    main()
