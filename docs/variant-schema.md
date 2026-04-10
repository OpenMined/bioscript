# Variant Schema

This is the compact core schema for a bioscript variant record.

Use it as the canonical stored form. Keep it small enough to:

- validate cleanly
- generate `bioscript.variant(...)` calls
- enrich later with external research if needed

## Schema Identity

```yaml
schema: "bioscript:variant"
version: "1.0"
```

`bioscript:variant` is a reasonable name for this object. It is explicit, stable, and leaves room for later schema types like `bioscript:panel` or `bioscript:report`.

## Minimal Shape

```yaml
schema: "bioscript:variant"
version: "1.0"
variant_id: "APOL1_G1_rs73885319"

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
  canonical_alt: "G"
```

## Required Fields

- `schema`
- `version`
- `variant_id`
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
    - "rs1317778148"
    - "rs143830837"
  aliases: []
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
- `canonical_alt` optional
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
  canonical_alt: "A"
```

Example deletion:

```yaml
alleles:
  kind: "deletion"
  ref: "I"
  alts:
    - "D"
  canonical_alt: "D"
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

`canonical_alt` means: this is the specific alternate allele the file is primarily about.

If the file is only a locus-level record and you do not want to choose one allele, omit `canonical_alt`.

## Optional Metadata

These fields are optional and can be added without changing the core shape:

- `label`
- `gene`
- `summary`
- `research`
- `clinical`

## `research`

Optional instructions or tags for enrichment tooling.

```yaml
research:
  tasks:
    - "Confirm canonical dbSNP aliases"
    - "Find ClinVar and Ensembl links"
  tags:
    - "apol1"
    - "pgx"
```

`tasks` is for downstream research automation, not for runtime execution.

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
schema: "bioscript:variant"
version: "1.0"
variant_id: "APOL1_G2"
label: "APOL1 G2 deletion"
gene: "APOL1"

identifiers:
  rsids:
    - "rs71785313"
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
  canonical_alt: "D"
  deletion_length: 6
  motifs:
    - "TTATAA"
    - "ATAATT"
```

maps cleanly to:

```python
bioscript.variant(
    rsid=["rs71785313", "rs1317778148", "rs143830837"],
    grch37="22:36662046-36662051",
    grch38="22:36266000-36266005",
    ref="I",
    alt="D",
    kind="deletion",
    deletion_length=6,
    motifs=["TTATAA", "ATAATT"],
)
```
