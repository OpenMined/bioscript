import argparse
import csv
from collections import defaultdict

STATUSES = ["G0", "G1", "G2"]


def parse_args():
    parser = argparse.ArgumentParser(
        description="Aggregate APOL1 status statistics (G0/G1/G2 counts + hetero/homo)"
    )
    parser.add_argument("--input", required=True, help="Input TSV from combined classifier output")
    parser.add_argument("--output", required=True, help="Output TSV path")
    return parser.parse_args()


def parse_status(apol1_status):
    status = (apol1_status or "").strip()
    if not status:
        return []
    parts = [part.strip() for part in status.split("/")]
    if len(parts) != 2:
        return []
    return parts


def main():
    args = parse_args()

    # Track one APOL1 status per participant (all rows for a participant should match).
    participant_status = {}

    with open(args.input, "r", encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            participant_id = (row.get("participant_id") or "").strip()
            apol1_status = (row.get("apol1_status") or "").strip()
            if not participant_id or not apol1_status:
                continue
            if participant_id not in participant_status:
                participant_status[participant_id] = apol1_status

    stats = defaultdict(lambda: {"count": 0, "hetero": 0, "homo": 0})
    total_successful_reads = 0

    for apol1_status in participant_status.values():
        alleles = parse_status(apol1_status)
        called = [allele for allele in alleles if allele in STATUSES]

        # "number" is the total number of successful allele reads across participants.
        total_successful_reads += len(called)

        for allele in called:
            stats[allele]["count"] += 1

        # Track per-status zygosity participant counts when both alleles are callable.
        if len(called) == 2:
            if called[0] == called[1]:
                stats[called[0]]["homo"] += 1
            else:
                stats[called[0]]["hetero"] += 1
                stats[called[1]]["hetero"] += 1

    with open(args.output, "w", encoding="utf-8", newline="") as out_fh:
        fieldnames = ["status", "count", "number", "hetero", "homo", "frequency"]
        writer = csv.DictWriter(out_fh, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()

        for status in STATUSES:
            count = stats[status]["count"]
            frequency = (count / total_successful_reads) if total_successful_reads else 0
            writer.writerow(
                {
                    "status": status,
                    "count": count,
                    "number": total_successful_reads,
                    "hetero": stats[status]["hetero"],
                    "homo": stats[status]["homo"],
                    "frequency": f"{frequency:.6f}",
                }
            )


if __name__ == "__main__":
    main()
