#!/usr/bin/env bash
# Regenerate the tiny CRAM fixtures used by bioscript-formats integration tests.
# Run: ./generate_fixtures.sh   (requires samtools)
#
# Produces:
#   mini.fa, mini.fa.fai          — 3 kb synthetic reference (contig chr_test)
#   mini.cram, mini.cram.crai     — 2000 reads covering pos 500..2500
#   mini_bad_ref.fa, .fai          — reference with a single-base mutation at 2800
#                                    (inside the slice span → triggers MD5 mismatch)
#
# The fixtures are deterministic: regenerating produces byte-identical output
# (modulo samtools version differences in the CRAM container encoding).
set -euo pipefail
cd "$(dirname "$0")"

python3 - <<'PY'
import random, pathlib
random.seed(20260414)
bases = "ACGT"
contig = "".join(random.choices(bases, k=3000))
# Keep the reference all-uppercase so MD5 matches after noodles normalization.
with open("mini.fa", "w") as fh:
    fh.write(">chr_test\n")
    for i in range(0, len(contig), 60):
        fh.write(contig[i:i+60] + "\n")

# Tampered reference: flip one base at pos 2800 (0-indexed 2799).
tampered_list = list(contig)
idx = 2799
original = tampered_list[idx]
tampered_list[idx] = next(b for b in "ACGT" if b != original)
tampered = "".join(tampered_list)
with open("mini_bad_ref.fa", "w") as fh:
    fh.write(">chr_test\n")
    for i in range(0, len(tampered), 60):
        fh.write(tampered[i:i+60] + "\n")

# 2000 reads, 50bp each, uniformly covering pos 500..2500 (1-based inclusive start).
# Sequences exactly match the reference (no mismatches) so depth parity is clean.
reads = []
for i in range(2000):
    pos = 500 + (i * (2000 // 2000)) % 2000  # spread across 500..2499
    pos = 500 + (i * 2000 // 2000)           # i.e. 500..2499 in order
    pos = 500 + i                             # 1-based; 500..2499
    seq = contig[pos-1:pos-1+50]
    reads.append((f"r{i:05d}", pos, seq))
# Stable sort by pos (already sorted since pos = 500+i).
sam = []
sam.append("@HD\tVN:1.6\tSO:coordinate")
sam.append("@SQ\tSN:chr_test\tLN:3000")
for name, pos, seq in reads:
    qual = "I" * 50
    sam.append(f"{name}\t0\tchr_test\t{pos}\t60\t50M\t*\t0\t0\t{seq}\t{qual}")
pathlib.Path("mini.sam").write_text("\n".join(sam) + "\n")
PY

samtools faidx mini.fa
samtools faidx mini_bad_ref.fa

samtools view -C --no-PG -T mini.fa -o mini.cram mini.sam
samtools index mini.cram

rm -f mini.sam

echo "fixtures written:"
ls -la mini.* mini_bad_ref.*
