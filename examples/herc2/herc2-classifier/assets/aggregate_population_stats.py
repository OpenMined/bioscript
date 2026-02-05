import argparse
import csv
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

def main():
    args = parse_args()

    stats = defaultdict(lambda: {
        "allele_count": 0,
        "allele_number": 0,
        "num_homo": 0,
        "num_hetero": 0,
        "rsid": None,
    })

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

            entry = stats[locus_key]
            entry["allele_count"] += alt_count
            entry["allele_number"] += allele_number
            if alt_count == 2:
                entry["num_homo"] += 1
            elif alt_count == 1:
                entry["num_hetero"] += 1

            rsid = row.get("rsid")
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
        for locus_key in sorted(stats.keys()):
            entry = stats[locus_key]
            allele_number = entry["allele_number"]
            allele_freq = (entry["allele_count"] / allele_number) if allele_number else 0
            writer.writerow({
                "locus_key": locus_key,
                "allele_count": entry["allele_count"],
                "allele_number": allele_number,
                "num_homo": entry["num_homo"],
                "num_hetero": entry["num_hetero"],
                "allele_freq": f"{allele_freq:.6f}",
                "rsid": entry["rsid"],
            })

if __name__ == "__main__":
    main()
