# Variant Schema

This is the compact core schema for a bioscript variant record.

Use it as the canonical stored form. Keep it small enough to:

- validate cleanly
- generate `bioscript.variant(...)` calls
- carry variant identity plus typed findings

## Schema Identity

```yaml
schema: "bioscript:variant:1.0"
version: "1.0"
```

## Minimal Shape

```yaml
schema: "bioscript:variant:1.0"
version: "1.0"
name: "traits-common-rs671-G-A"
tags:
  - "type:trait"

identifiers:
  rsids:
    - "rs671"

coordinates:
  grch37:
    chrom: "12"
    pos: 112241766
  grch38:
    chrom: "12"
    pos: 111803962

alleles:
  kind: "snv"
  ref: "G"
  alts:
    - "A"
```

## Required Fields

- `schema`
- `version`
- `name`
- `alleles.kind`
- `alleles.ref`
- `alleles.alts`

At least one of these must also exist:

- `identifiers`
- `coordinates`

## Top-Level `tags`

Optional free-form classification tags for filtering and compatibility hints.

Example:

```yaml
tags:
  - "type:trait"
  - "validated:23andme:v5"
```

Guidelines:

- use simple lowercase strings
- keep the vocabulary small
- `type:trait` is the main broadly useful tag for common trait records
- only add `validated:*` tags after direct assay test runs against local test data
- absence of a `validated:*` tag does not mean unsupported

## Validation Rules

- `schema` must be `bioscript:variant:1.0`
- `version` must be `1.0`
- `name` is required
- `identifiers.rsids` and `identifiers.aliases`, if present, must look like `rs123`
- `coordinates.*.chrom` must be one of:
  - `1` through `22`
  - `X`
  - `Y`
  - `MT`
- `coordinates.*` must use either:
  - `pos`
  - or `start` and `end`
- `pos`, `start`, and `end` must be integers
- `start` and `end` must be `>= 1`
- `end` must be `>= start`
- `alleles.kind` must be one of:
  - `snv`
  - `deletion`
  - `insertion`
  - `indel`
- stored allele values must be biological alleles, not symbolic `I` / `D`
- for `kind: snv`, `ref` and every `alt` must be single-base `A`, `C`, `G`, or `T`
- each finding must have its own `schema`
- if a finding has `alt`, it must exist in `alleles.alts`
- provenance URLs must be valid `http` or `https` URLs

## Core Model

### `identifiers`

External IDs for the variant.

Current fields:

- `rsids`
- `aliases`

Example:

```yaml
identifiers:
  rsids:
    - "rs71785313"
  aliases: []
```

### `coordinates`

Assembly-specific genomic coordinates.

Supported assemblies now:

- `grch37`
- `grch38`

Coordinate shapes:

```yaml
coordinates:
  grch38:
    chrom: "22"
    pos: 36265860
```

```yaml
coordinates:
  grch38:
    chrom: "22"
    start: 36266000
    end: 36266005
```

Rules:

- use `pos` for single-base sites
- use `start` and `end` for spans such as deletions or complex indels
- do not include HGVS here

### `alleles`

The biological allele definition.

Fields:

- `kind`: `snv | deletion | insertion | indel`
- `ref`
- `alts`
- `deletion_length` optional
- `insertion_sequence` optional
- `motifs` optional
- `equivalent_coordinates` optional
- `notes` optional

Stored YAML should describe the biological allele. Do not use symbolic `I` / `D` allele values in this schema.

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

### `findings`

Optional typed interpretation records attached to the variant.

Each finding is a small object with its own `schema`. Parsers may validate schemas they understand and ignore schemas they do not understand.

Minimal finding envelope:

```yaml
findings:
  - schema: "bioscript:trait:1.0"
    alt: "A"
    summary: "Associated with alcohol flushing after alcohol exposure."
```

Envelope fields:

- `schema` required
- `alt` optional, but required for allele-specific findings; use `"*"` when the finding applies to any alternate allele at a multiallelic locus
- `label` optional
- `summary` optional
- `notes` optional

Unknown finding schemas are allowed.

## Optional Metadata

These fields are optional and can be added without changing the core shape:

- `tags`
- `label`
- `gene`
- `summary`
- `findings`
- `provenance`

## `provenance`

Optional source metadata for the variant record and its findings.

Example:

```yaml
provenance:
  sources:
    - kind: "database"
      label: "dbSNP"
      url: "https://www.ncbi.nlm.nih.gov/snp/rs671"
      fields:
        - "identifiers.rsids"
        - "coordinates.grch37"
        - "coordinates.grch38"
        - "alleles"
```

## Deliberately Out Of Core Schema

These are intentionally not part of the compact core schema right now:

- transcript HGVS
- protein HGVS
- domain-specific nested blocks like `clinical`

Those can be expressed as typed findings instead.
