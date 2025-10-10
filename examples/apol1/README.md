# APOL1 Classifier Example

Classification of APOL1 genotypes (G0, G1, G2) for kidney disease risk assessment.

## Files

- `classify_apol1.py` - APOL1 classifier implementation
- `samplesheet.csv` - Example samplesheet with 3 participants
- `test_snps*.txt` - Test SNP files for each participant
- `process_samplesheet.sh` - Batch processing script

## Quick Test

From repository root:

```bash
# Test with local virtualenv (using uv)
./test_local.sh

# Test with Docker
./test_docker.sh
```

## Usage

### Single participant

```bash
bioscript classify participant_id=P001 \
    --file=test_snps.txt \
    classify_apol1.py \
    --out=tsv
```

### Batch processing

```bash
# Local (requires bioscript installed)
./process_samplesheet.sh samplesheet.csv classify_apol1.py > results.tsv

# Docker
docker run --rm -v $(pwd):/data -w /data \
    ghcr.io/openmined/bioscript:latest \
    bash process_samplesheet.sh samplesheet.csv classify_apol1.py > results.tsv
```

## Expected Output

```tsv
participant_id	APOL1
P001	G1/G0
P002	G0/G0
P003	G2/G1
```

## Classifier Details

**Variants tested**:
- rs73885319 (A>G) - G1 variant site 1
- rs60910145 (T>C) - G1 variant site 2
- rs71785313 (InDel) - G2 variant

**Genotypes**:
- **G0/G0** - No risk variants (healthy)
- **G1/G0** - Heterozygous G1 (low risk)
- **G1/G1** - Homozygous G1 (moderate risk)
- **G2/G0** - Heterozygous G2 (moderate risk)
- **G2/G1** - Compound heterozygous (high risk)
- **G2/G2** - Homozygous G2 (high risk)

## Creating Your Own Classifier

See `classify_apol1.py` for the implementation pattern. Your classifier must export:

```python
# Required exports
variant_calls = [rs1, rs2, rs3]  # List of VariantCall objects
classifier = MyClassifier()       # Classifier instance

# Optional export
name = "MyTest"  # Column name (defaults to script filename)
```
