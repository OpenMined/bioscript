import argparse
import csv
import re
from collections import defaultdict

def parse_args():
    parser = argparse.ArgumentParser(description="Aggregate population allele statistics")
    parser.add_argument("--input", required=True, help="Input TSV from combined classifier output")
    parser.add_argument("--output", required=True, help="Output TSV path")
    return parser.parse_args()

def locus_key_for_row(row):
    chrom = (row.get("chromosome") or "").strip()
    pos = (row.get("position") or "").strip()
    ref = (row.get("ref") or "").strip()
    alt = (row.get("alt") or "").strip()
    if not (chrom and pos and ref and alt):
        return None
    return f"{chrom}-{pos}-{ref}-{alt}"

def to_int(value):
    try:
        return int(value)
    except Exception:
        return 0

def parse_alleles(genotype, alleles):
    if not genotype:
        return []
    genotype = genotype.strip()
    if not genotype:
        return []

    # Split on common genotype separators
    if "/" in genotype or "|" in genotype:
        parts = [p for p in re.split(r"[\/|]", genotype) if p]
        return parts

    # If alleles are same length and genotype is a concat of two alleles, split evenly
    if alleles:
        lengths = {len(a) for a in alleles if a}
        if len(lengths) == 1:
            allele_len = lengths.pop()
            if allele_len > 0 and len(genotype) == 2 * allele_len:
                return [genotype[:allele_len], genotype[allele_len:]]

    # Fallback for single-base alleles like "AA" or "AG"
    if len(genotype) == 2:
        return [genotype[0], genotype[1]]

    return [genotype]

def main():
    args = parse_args()

    # stats keyed by (locus_base, allele)
    stats = defaultdict(lambda: {
        "allele_count": 0,
        "num_homo": 0,
        "num_hetero": 0,
        "rsid": None,
    })
    locus_called_samples = defaultdict(int)

    with open(args.input, "r", encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            match_type = (row.get("match_type") or "").strip()
            if match_type.lower() == "no call":
                continue

            locus_key = locus_key_for_row(row)
            if locus_key is None:
                continue

            alt_count = to_int(row.get("alt_count"))
            ref_count = to_int(row.get("ref_count"))
            allele_number = ref_count + alt_count

            if allele_number <= 0:
                continue

            chrom = (row.get("chromosome") or "").strip()
            pos = (row.get("position") or "").strip()
            ref = (row.get("ref") or "").strip()
            alt = (row.get("alt") or "").strip()
            locus_base = f"{chrom}-{pos}"

            alleles = []
            if ref:
                alleles.append(ref)
            if alt:
                alleles.extend([a.strip() for a in alt.split(",") if a.strip()])

            genotype = (row.get("genotype") or "").strip()
            gt_alleles = parse_alleles(genotype, alleles)
            if not gt_alleles:
                continue

            locus_called_samples[locus_base] += 1

            rsid = row.get("rsid")
            is_homo = len(gt_alleles) == 2 and gt_alleles[0] == gt_alleles[1]

            for allele in gt_alleles:
                key = (locus_base, allele)
                entry = stats[key]
                entry["allele_count"] += 1
                if is_homo and allele == gt_alleles[0]:
                    entry["num_homo"] += 1
                elif not is_homo:
                    entry["num_hetero"] += 1
                if rsid and not entry["rsid"]:
                    entry["rsid"] = rsid

    with open(args.output, "w", encoding="utf-8", newline="") as out_fh:
        fieldnames = [
            "locus_key",
            "allele_count",
            "allele_number",
            "num_homo",
            "num_hetero",
            "allele_freq",
            "rsid",
        ]
        writer = csv.DictWriter(out_fh, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for (locus_base, allele) in sorted(stats.keys()):
            entry = stats[(locus_base, allele)]
            called = locus_called_samples.get(locus_base, 0)
            allele_number = called * 2
            allele_freq = (entry["allele_count"] / allele_number) if allele_number else 0
            writer.writerow({
                "locus_key": f"{locus_base}-{allele}",
                "allele_count": entry["allele_count"],
                "allele_number": allele_number,
                "num_homo": entry["num_homo"],
                "num_hetero": entry["num_hetero"],
                "allele_freq": f"{allele_freq:.6f}",
                "rsid": entry["rsid"],
            })

if __name__ == "__main__":
    main()
