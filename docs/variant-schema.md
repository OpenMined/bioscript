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
- `observed_alts` optional
- `deletion_length` optional
- `insertion_sequence` optional
- `motifs` optional
- `equivalent_coordinates` optional
- `notes` optional

Stored YAML should describe the biological allele. Do not use symbolic `I` / `D` allele values in this schema.

`alts` is the curated set of alternate alleles that matter for this catalogue entry and should drive app flagging. `observed_alts` is the full set of source-reported alternate alleles observed at the same locus, usually from dbSNP. If `observed_alts` is omitted, tools treat `alts` as both curated and observed.

When a source such as dbSNP reports multiple alternates but only one has the clinical or PGx evidence being catalogued, keep that allele in `alts` and put the full dbSNP set in `observed_alts`.

Example SNV:

```yaml
alleles:
  kind: "snv"
  ref: "G"
  alts:
    - "T"
  observed_alts:
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
- `binding` optional; use it when report logic should match a specific variant observation field instead of relying on `alt`
- `label` optional
- `summary` optional
- `notes` optional

Variant-bound PGx sidecar include:

```yaml
findings:
  - schema: "bioscript:pgx-summary:1.0"
    id: "rs123_pgx_sidecar"
    include: "rs123-pgx.yaml"
    notes: "Detailed PGx findings are stored in the sidecar file."
```

Sidecar files use `schema: "bioscript:pgx-findings:1.0"` and contain dense PGx evidence for one variant. This keeps large ClinPGx/PharmGKB annotation tables out of the core variant identity file. A sidecar may contain both summary annotations and drug label annotations.

Summary annotations are variant/drug evidence interpretations. They use `schema: "bioscript:pgx-summary:1.0"` and normally include evidence levels, phenotype categories, and genotype-specific effects.

```yaml
findings:
  - schema: "bioscript:pgx-summary:1.0"
    id: "example_summary_annotation"
    authority_type: "evidence_summary"
    drugs:
      - name: "example drug"
    phenotype_categories:
      - "Toxicity"
    evidence_level: "3"
    evidence:
      source: "ClinPGx"
      kind: "summary_annotation"
      id: "1448427005"
      url: "https://www.clinpgx.org/variant/PA.../summaryAnnotation"
    effects:
      - id: "example_variant_bound_pgx_finding_alt_carrier"
        label: "C carrier"
        binding:
          source: "variant"
          variant: "rs123.yaml"
          allele: "C"
          operator: "dosage_in"
          values: [1, 2]
          description: "Applies when the participant carries one or two copies of C."
        text: "Reportable text for this allele dosage."
```

Drug label annotations are regulatory label statements. They use `schema: "bioscript:pgx-label:1.0"` and should carry regulatory/action fields instead of summary evidence levels.

```yaml
findings:
  - schema: "bioscript:pgx-label:1.0"
    id: "example_label_annotation"
    authority_type: "regulatory_label"
    genes:
      - "ABCG2"
    drugs:
      - name: "rosuvastatin"
        aliases:
          - "Crestor"
    regulatory_sources:
      - "FDA"
    pgx_action_level: "testing_recommended"
    prescribing_actions:
      - "dose_adjustment"
    evidence:
      source: "ClinPGx"
      kind: "label_annotation"
      id: "PA..."
      url: "https://www.clinpgx.org/variant/PA.../labelAnnotation"
    notes: "Regulatory drug label annotation."
```

Supported binding operators are:

- `equals` and `in` for matching literal analysis outputs or observation fields
- `dosage_equals` and `dosage_in` for variant allele dosage, where `allele` is the reference allele or one of the alternate alleles and dosage values are `0`, `1`, or `2`

Known PGx finding schemas are:

- `bioscript:pgx-summary:1.0` for ClinPGx/PharmGKB summary annotations
- `bioscript:pgx-label:1.0` for ClinPGx/PharmGKB drug label annotations
- `bioscript:pgx:1.0` is legacy; prefer one of the two specific schemas above

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
