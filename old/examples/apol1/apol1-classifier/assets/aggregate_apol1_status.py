import argparse
import csv


def parse_args():
    parser = argparse.ArgumentParser(description="Aggregate APOL1 participant statuses")
    parser.add_argument("--input", required=True, help="Input TSV from combined classifier output")
    parser.add_argument("--output", required=True, help="Output TSV path")
    return parser.parse_args()


def main():
    args = parse_args()

    participant_rows = {}

    with open(args.input, "r", encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            participant_id = (row.get("participant_id") or "").strip()
            original_file = (row.get("filename") or "").strip()
            apol1_status = (row.get("apol1_status") or "").strip()

            if not participant_id:
                continue
            if participant_id not in participant_rows:
                participant_rows[participant_id] = {
                    "participant_id": participant_id,
                    "original_file": original_file,
                    "apol1_status": apol1_status,
                }

    with open(args.output, "w", encoding="utf-8", newline="") as out_fh:
        fieldnames = ["participant_id", "original_file", "apol1_status"]
        writer = csv.DictWriter(out_fh, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in participant_rows.values():
            writer.writerow(row)


if __name__ == "__main__":
    main()
