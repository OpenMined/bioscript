# Variant Schema

This is the compact core schema for a bioscript variant record.

Use it as the canonical stored form. Keep it small enough to:

- validate cleanly
- generate `bioscript.variant(...)` calls
- enrich later with external research if needed

## Schema Identity

```yaml
name: "PROJECT-SOURCE-rsNNNNNN-REF-ALT"
version: "1.0"
schema: "bioscript:variant:1.0"
```

- `name`: a human-readable identifier for this specific variant record, typically `PROJECT-SOURCE-rsid-REF-ALT`
- `version`: record version
- `schema`: schema type and version combined as `bioscript:variant:1.0`

## Minimal Shape

```yaml
name: "APOL1-test-rs73885319-A-G"
version: "1.0"
schema: "bioscript:variant:1.0"
gene: "APOL1"

identifiers:
  rsids:
    - "rs73885319"

coordinates:
  grch37:
    chrom: "22"
    pos: 36661906
  grch38:
    chrom: "22"
    pos: 36265860

alleles:
  kind: "snv"
  ref: "A"
  alts:
    - "G"
```

## Required Fields

- `name`
- `version`
- `schema`
- `alleles.kind`
- `alleles.ref`
- `alleles.alts`

At least one of these must also exist:

- `identifiers`
- `coordinates`

## Core Model

### `identifiers`

External IDs for the variant.

Current fields:

- `rsids`
- `aliases`

Each is a simple string list.

Example:

```yaml
identifiers:
  rsids:
    - "rs71785313"
  aliases:
    - "rs1317778148"
    - "rs143830837"
```

### `coordinates`

Assembly-specific genomic coordinates.

Supported assemblies now:

- `grch37`
- `grch38`

Coordinate shape:

```yaml
coordinates:
  grch38:
    chrom: "22"
    start: 36266000
    end: 36266005
```

For single-base variants, you may use `pos` instead of `start` and `end`:

```yaml
coordinates:
  grch38:
    chrom: "4"
    pos: 88131171
```

Rules:

- use `pos` for single-base sites when you want the compact form
- use `start` and `end` for spans such as deletions
- do not include HGVS here

### `alleles`

The biological allele definition.

Fields:

- `kind`: `snv | deletion | insertion | indel | other`
- `ref`
- `alts`
- `deletion_length` optional
- `insertion_sequence` optional
- `motifs` optional
- `equivalent_coordinates` optional
- `notes` optional

Example SNV:

```yaml
alleles:
  kind: "snv"
  ref: "G"
  alts:
    - "A"
    - "C"
    - "T"
```

Example deletion:

```yaml
alleles:
  kind: "deletion"
  ref: "I"
  alts:
    - "D"
  deletion_length: 6
  motifs:
    - "TTATAA"
    - "ATAATT"
  equivalent_coordinates:
    grch37:
      - "22:36662042-36662047"
      - "22:36662046-36662051"
    grch38:
      - "22:36265996-36266001"
      - "22:36266000-36266005"
```

`alts` lists all known alternate alleles for this variant. Use `findings` to annotate which specific alt is associated with a particular result.

### `findings`

Optional list of alt-specific findings from studies or databases. Each entry ties a specific alternate allele to a result, with a source URL for verification.

Fields per entry:

- `alt` required — the specific alternate allele this finding is about (use `"*"` when the finding applies to all alts or the specific alt is ambiguous, e.g. multiallelic indels)
- `notes` required — free-text description of the finding
- `source` optional — URL to the paper, database entry, or evidence

Example:

```yaml
findings:
  - alt: "T"
    notes: "p.Pro7Leu missense — associated with increased GLP-1 receptor agonist weight-loss efficacy (~0.76 kg per allele, p=2.9e-10)"
    source: "https://www.nature.com/articles/s41586-026-10330-z"
  - alt: "T"
    notes: "PharmGKB level 4 Efficacy annotation for rs10305420 (GLP1R) and liraglutide"
    source: "https://www.pharmgkb.org/variant/PA166157167"
```

Use `findings` when:
- A paper reports an association with a specific allele
- Different alts have different clinical significance
- You want to record which alt matters for a panel or assay

## Optional Metadata

These fields are optional and can be added without changing the core shape:

- `gene`
- `summary`
- `findings`
- `provenance`
- `clinical`

### `provenance`

Tracks where each field's data came from.

```yaml
provenance:
  sources:
    - kind: "paper"
      label: "Nature 2026 GLP-1 response study"
      url: "https://www.nature.com/articles/s41586-026-10330-z"
      fields:
        - "summary"
        - "findings[0]"
    - kind: "database"
      label: "dbSNP / NCBI variation services"
      url: "https://www.ncbi.nlm.nih.gov/snp/rs10305420"
      fields:
        - "identifiers.rsids"
        - "coordinates.grch37"
        - "coordinates.grch38"
        - "alleles"
```

Each source entry has:

- `kind`: `"paper"`, `"database"`, or other descriptor
- `label`: human-readable name
- `url`: verification link
- `fields`: list of dotpath references to fields derived from this source

## `clinical`

Optional clinical interpretation metadata.

Current supported optional domain:

- `clinical.pgx`

Example:

```yaml
clinical:
  pgx:
    gene: "ABCG2"
    source: "ClinPGx"
    variant_page_url: "https://www.clinpgx.org/variant/PA166156544/labelAnnotation"
    drug_labels:
      - source: "HCSC"
        title: "Annotation of HCSC Label for rosuvastatin and ABCG2, SLCO1B1"
        genes:
          - "ABCG2"
          - "SLCO1B1"
        drugs:
          - "rosuvastatin"
        pgx_level: "Actionable PGx"
        actionable: true
        tags:
          - "Dosing Info"
          - "Prescribing Info"
```

## Deliberately Out Of Core Schema

These are intentionally not part of the compact core schema right now:

- `hgvs_genomic`
- transcript HGVS
- protein HGVS

Those can usually be derived or enriched later. They are useful annotation data, but they make the base files noisier and harder to validate.

## Mapping To `bioscript.variant(...)`

This YAML:

```yaml
name: "APOL1-test-G2-deletion"
version: "1.0"
schema: "bioscript:variant:1.0"
gene: "APOL1"

identifiers:
  rsids:
    - "rs71785313"
  aliases:
    - "rs1317778148"
    - "rs143830837"

coordinates:
  grch37:
    chrom: "22"
    start: 36662046
    end: 36662051
  grch38:
    chrom: "22"
    start: 36266000
    end: 36266005

alleles:
  kind: "deletion"
  ref: "I"
  alts:
    - "D"
  deletion_length: 6
  motifs:
    - "TTATAA"
    - "ATAATT"

findings:
  - alt: "D"
    notes: "G2 risk deletion for APOL1-associated nephropathy"
```

maps cleanly to:

```python
bioscript.variant(
    rsid="rs71785313",
    grch37="22:36662046-36662051",
    grch38="22:36266000-36266005",
    ref="I",
    alt="D",
    kind="deletion",
    deletion_length=6,
    motifs=["TTATAA", "ATAATT"],
)
```
