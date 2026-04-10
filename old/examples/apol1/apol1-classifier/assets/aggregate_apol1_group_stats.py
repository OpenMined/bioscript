import argparse
import csv
from collections import Counter

GROUPS = ["G0/G0", "G1/G0", "G1/G1", "G2/G0", "G2/G1", "G2/G2"]


def parse_args():
    parser = argparse.ArgumentParser(description="Count distinct APOL1 diploid classification groups")
    parser.add_argument("--input", required=True, help="Input TSV from combined classifier output")
    parser.add_argument("--output", required=True, help="Output TSV path")
    return parser.parse_args()


def normalize_status(apol1_status):
    status = (apol1_status or "").strip()
    if not status:
        return None
    parts = [p.strip() for p in status.split("/")]
    if len(parts) != 2:
        return None
    order = {"G2": 0, "G1": 1, "G0": 2}
    sorted_parts = sorted(parts, key=lambda p: order.get(p, 99))
    return f"{sorted_parts[0]}/{sorted_parts[1]}"


def main():
    args = parse_args()

    participant_status = {}
    with open(args.input, "r", encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            pid = (row.get("participant_id") or "").strip()
            status = (row.get("apol1_status") or "").strip()
            if not pid or not status:
                continue
            if pid not in participant_status:
                participant_status[pid] = normalize_status(status)

    counts = Counter(s for s in participant_status.values() if s)

    with open(args.output, "w", encoding="utf-8", newline="") as out_fh:
        writer = csv.DictWriter(out_fh, fieldnames=["classification", "count"], delimiter="\t")
        writer.writeheader()
        for group in GROUPS:
            writer.writerow({"classification": group, "count": counts.get(group, 0)})


if __name__ == "__main__":
    main()
