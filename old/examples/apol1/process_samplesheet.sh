#!/bin/bash
set -e

# Process samplesheet CSV and run classifiers on each participant
#
# Usage: ./process_samplesheet.sh <samplesheet.csv> <classifier1.py> [classifier2.py ...]

SAMPLESHEET="$1"
shift
CLASSIFIERS="$@"

# Get header from first row
FIRST_LINE=$(tail -n +2 "$SAMPLESHEET" | head -n 1)
IFS=',' read -r PARTICIPANT_ID SNP_FILE <<< "$FIRST_LINE"

# Print TSV header
bioscript classify --participant_id="${PARTICIPANT_ID}" \
    --file="$SNP_FILE" \
    $CLASSIFIERS \
    --out=tsv | head -n 1

# Process each row
tail -n +2 "$SAMPLESHEET" | while IFS=',' read -r PARTICIPANT_ID SNP_FILE; do
    bioscript classify --participant_id="${PARTICIPANT_ID}" \
        --file="$SNP_FILE" \
        $CLASSIFIERS \
        --out=tsv | tail -n +2
done
