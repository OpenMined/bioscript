# VNtyper Anonymized Test Data - Multi-Reference Dataset

**Version**: 2.1 | **Generated**: 2026-03-23 | **Status**: ✅ Ready

---

## Overview

Fully anonymized MUC1 VNTR test data for VNtyper, aligned to all six supported reference assemblies (hg19, hg38, GRCh37, GRCh38, hg19_ensembl, hg38_ensembl) using BWA-MEM. Includes regression guard samples for parameter stability testing.

**Total Files**: 116 | **8 Samples** (7 multi-reference + 1 hg38 regression guard)

---

## Directory Structure

```
tests/data/
├── example_XXXX_hg19_subset.bam (×7)          # Original hg19 subsets + indexes
├── example_40cf_hg38_subset.bam               # hg38 regression guard (Issue #156)
├── fastqs/
│   └── example_XXXX_hg19_subset_R{1,2}.fastq.gz (×14)
└── remapped/bwa/
    ├── hg19/          example_XXXX_hg19_bwa.bam + .bai (×7)
    ├── hg38/          example_XXXX_hg38_bwa.bam + .bai (×7)
    ├── GRCh37/        example_XXXX_GRCh37_bwa.bam + .bai (×7)
    ├── GRCh38/        example_XXXX_GRCh38_bwa.bam + .bai (×7)
    ├── hg19_ensembl/  example_XXXX_hg19_ensembl_bwa.bam + .bai (×7)
    └── hg38_ensembl/  example_XXXX_hg38_ensembl_bwa.bam + .bai (×7)
```

---

## Samples

| Pseudonym | Size | Reads | Type |
|-----------|------|-------|------|
| `example_6449` | 16M | ~167K | MUC1 mutant |
| `example_b178` | 3.7M | ~34K | MUC1 mutant |
| `example_6c28` | 16M | ~120K | MUC1 mutant |
| `example_dfc3` | 6.1M | ~68K | MUC1 mutant |
| `example_66bf` | 4.2M | ~40K | MUC1 mutant |
| `example_7a61` | 81M | ~985K | Negative control |
| `example_a5c1` | 4.8M | ~43K | MUC1 mutant + adVNTR |
| `example_40cf` | 3.3M | ~39K | Negative (GDP inflation guard, hg38) |

---

## Reference Assemblies

| Assembly | Type | Chromosome | MUC1 Region |
|----------|------|------------|-------------|
| **hg19** | UCSC | chr1 | chr1:155158000-155163000 |
| **hg38** | UCSC | chr1 | chr1:155184000-155194000 |
| **GRCh37** | NCBI | NC_000001.10 | NC_000001.10:155158000-155163000 |
| **GRCh38** | NCBI | NC_000001.11 | NC_000001.11:155184000-155194000 |
| **hg19_ensembl** | ENSEMBL | 1 | 1:155158000-155163000 |
| **hg38_ensembl** | ENSEMBL | 1 | 1:155184000-155194000 |

---

## Usage

### With Original BAMs (hg19)
```bash
vntyper pipeline --bam tests/data_anonymized/example_6449_hg19_subset.bam \
  --reference hg19 --output results/
```

### With FASTQs
```bash
vntyper pipeline \
  --fastq1 tests/data_anonymized/fastqs/example_6449_hg19_subset_R1.fastq.gz \
  --fastq2 tests/data_anonymized/fastqs/example_6449_hg19_subset_R2.fastq.gz \
  --reference hg38 --output results/
```

### With Remapped BAMs
```bash
# UCSC naming (chr1)
vntyper pipeline --bam tests/data_anonymized/remapped/bwa/hg38/example_6449_hg38_bwa.bam \
  --reference hg38 --output results/

# NCBI naming (NC_000001.11)
vntyper pipeline --bam tests/data_anonymized/remapped/bwa/GRCh38/example_6449_GRCh38_bwa.bam \
  --reference GRCh38 --output results/

# ENSEMBL naming (1)
vntyper pipeline --bam tests/data_anonymized/remapped/bwa/hg38_ensembl/example_6449_hg38_ensembl_bwa.bam \
  --reference hg38_ensembl --output results/
```

---

## Testing

```bash
# Run all integration tests
pytest tests/test_integration.py -v

# Test specific sample
pytest tests/test_integration.py -k "example_6449" -v

# Unit test (FASTQ)
pytest tests/test_vntyper.py::test_fastq_shark -v
```

---

## Verification

```bash
# Check BAM integrity
samtools quickcheck tests/data_anonymized/example_6449_hg19_subset.bam

# View read count
samtools view -c tests/data_anonymized/example_6449_hg19_subset.bam

# Check indexes
samtools idxstats tests/data_anonymized/example_6449_hg19_subset.bam | head
```

---

## File Sizes

| Category | Files | Size |
|----------|-------|------|
| Original Subset BAMs | 16 | ~135 MB |
| FASTQ Files | 14 | ~140 MB |
| Remapped BAMs (6 refs) | 84 | ~420 MB |
| Metadata | 3 | ~1 MB |
| **Total** | **117** | **~696 MB** |

---

## Metadata Files

- `pseudonymization_table.csv` - Original → Pseudonym mapping
- `pseudonymization_output.json` - Complete file manifest with MD5 checksums
- `pseudonymization.log` - Generation log

---

## Notes

- All samples are MUC1 region subsets (±5kb) plus unmapped reads
- Read names anonymized (flowcell IDs hashed)
- BWA-MEM alignment (v0.7.17+)
- Paired-end reads preserved
- Compatible with VNtyper v2.0+
- `example_40cf` is an hg38-native sample added as a regression guard for Kestrel parameter stability (Issue #156: GDP inflation with maxhapstates/maxalignstates > 50). Must remain Negative with current parameters.

---

## License

Test Data: CC-BY-4.0 | Code: MIT License

---

**Documentation**: https://github.com/hassansaei/VNtyper
**Issues**: https://github.com/hassansaei/VNtyper/issues
