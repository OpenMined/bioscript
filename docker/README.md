# BioScript Docker Image

Lightweight Docker image for running BioScript genetic variant classifiers.

## Using Pre-built Images

Pull from GitHub Container Registry:

```bash
docker pull ghcr.io/openmined/bioscript:latest
# or specific version
docker pull ghcr.io/openmined/bioscript:0.1.0
```

Images are automatically built and pushed on every commit to main and on releases.

## Building Locally

From the repository root:

```bash
docker build -f docker/Dockerfile -t bioscript:0.1.0 .
```

Or use the build script:

```bash
cd docker
./build.sh
```

## Usage

### Run a classifier

```bash
docker run --rm -v $(pwd)/examples/apol1:/data -w /data \
  ghcr.io/openmined/bioscript:latest \
  bioscript classify participant_id=P001 \
  --file=test_snps.txt \
  classify_apol1.py \
  --out=tsv
```

### Interactive shell

```bash
docker run --rm -it -v $(pwd)/examples/apol1:/data -w /data \
  ghcr.io/openmined/bioscript:latest /bin/bash
```

### Process a samplesheet

```bash
docker run --rm -v $(pwd)/examples/apol1:/data -w /data \
  ghcr.io/openmined/bioscript:latest \
  bash process_samplesheet.sh samplesheet.csv classify_apol1.py > results.tsv
```

## Image Details

- **Base Image**: `python:3.13-slim`
- **Python Version**: 3.13
- **Package Manager**: uv (fast Python package installer)
- **Size**: ~150MB (slim base + bioscript + dependencies)

## Nextflow Integration

Example Nextflow process:

```groovy
process CLASSIFY_APOL1 {
    container 'ghcr.io/openmined/bioscript:latest'

    input:
    tuple val(participant_id), path(snp_file)
    path(classifier)

    output:
    path("${participant_id}_results.tsv")

    script:
    """
    bioscript classify participant_id=${participant_id} \\
        --file=${snp_file} \\
        ${classifier} \\
        --out=tsv > ${participant_id}_results.tsv
    """
}
```

### Full Nextflow Example

```groovy
#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.samplesheet = "samplesheet.csv"
params.classifiers = ["classify_apol1.py"]

process CLASSIFY_VARIANTS {
    container 'ghcr.io/openmined/bioscript:latest'
    publishDir "results", mode: 'copy'

    input:
    tuple val(participant_id), path(snp_file)
    path(classifiers)

    output:
    path("${participant_id}_results.tsv")

    script:
    def classifier_args = classifiers.collect{ it.getName() }.join(' ')
    """
    bioscript classify participant_id=${participant_id} \\
        --file=${snp_file} \\
        ${classifier_args} \\
        --out=tsv > ${participant_id}_results.tsv
    """
}

workflow {
    // Read samplesheet
    Channel
        .fromPath(params.samplesheet)
        .splitCsv(header: true)
        .map { row -> tuple(row.participant_id, file(row.snp_file)) }
        .set { samples }

    // Prepare classifiers
    classifiers = Channel.fromPath(params.classifiers)

    // Run classification
    CLASSIFY_VARIANTS(samples, classifiers.collect())
}
```
