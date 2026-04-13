# Assay Package Schema

This document defines the portable assay package contract for `bioscript:assay`.

The goal is to keep assay packages:

- portable across clients
- declarative rather than programmable
- expressive enough for third-party publication
- small enough to validate and render consistently

## Design Rules

1. `assay.yaml` is the source of truth for assay metadata.
2. Clients own rendering behavior, not assay-specific app code.
3. YAML may provide structured content, but it must not become a UI DSL.
4. Variant-level evidence belongs in variant YAML files.
5. Assay-level framing belongs in `assay.yaml`.

## Minimal Shape

```yaml
schema: "bioscript:assay"
version: "1.0"

assay_id: "herc2_eye_color"
label: "HERC2 Eye Color"
summary: "Classify the main HERC2 eye-color signal from one defining site."

metadata:
  category: "traits"
  tags:
    - "herc2"
    - "eye-color"
  disclaimer: "Research use only. Not medical advice."

package:
  assay_version: "0.1.0"
  source_of_truth: "script"

ui:
  template: "binary-call"
  version: "1.0"

compatibility:
  works_with:
    - "23andme"
    - "wgs"
    - "vcf"
  assemblies:
    - "grch37"
    - "grch38"
  notes: "Checks one major HERC2 eye-color signal."

privacy:
  mode: "local_only"
  uploads_data: false
  stores_results_locally: true
  external_urls: []

inputs:
  catalogue: "catalogue.yaml"

implementation:
  kind: "script"
  path: "herc2.py"
  runtime: "bioscript"

outputs:
  format: "tsv"
```

## Required Fields

- `schema`
- `version`
- `assay_id`
- `label`
- `summary`
- `metadata.category`
- `ui.template`
- `ui.version`
- `inputs.catalogue`
- `implementation.kind`
- `implementation.path`
- `implementation.runtime`
- `outputs.format`

## `metadata`

Portable assay identity and classification.

```yaml
metadata:
  category: "risk"
  tags:
    - "apol1"
    - "kidney"
  disclaimer: "Research use only. Not medical advice."
```

### `metadata.category`

Use the shared client-facing assay taxonomy:

- `traits`
- `ancestry`
- `pgx`
- `risk`

Do not use client-specific variants such as `health-risk`.

### `metadata.tags`

Simple search and grouping tags.

### `metadata.disclaimer`

Plain-text assay disclaimer shown by clients.

## `ui`

Declares the renderer family a client should use.

```yaml
ui:
  template: "risk-classifier"
  version: "1.0"
```

This is intentionally small:

- `template` selects a supported renderer family
- `version` declares the content contract that renderer expects

This is not a templating language. The client still owns layout, styling, and interaction behavior.

### Recommended Initial Templates

- `binary-call`
- `risk-classifier`
- `variant-panel`
- `drug-response`

Clients may support only a subset and should fail clearly on unsupported templates or versions.

## `interpretation`

Optional assay-level result content for portable third-party publishing.

This block is allowed because clients cannot be expected to ship assay-specific prose for every published assay. Keep it constrained to fixed content slots rather than arbitrary rendering instructions.

```yaml
interpretation:
  states:
    matched:
      headline: "Eye-color signal detected"
      body: "This assay detected the main HERC2 signal associated with lighter eye color."
      caveat: "Eye color is influenced by multiple genes."
    normal:
      headline: "No flagged signal found"
      body: "The checked site was present, but the lighter-eye-associated signal was not detected."
    missing:
      headline: "Not enough data"
      body: "This file did not include enough of the expected sites for a confident result."
    partial:
      headline: "Partial result"
      body: "Some expected rows were present, but the file did not cover the full assay."
```

Rules:

- only use these fixed state keys:
  - `matched`
  - `normal`
  - `missing`
  - `partial`
- each state may define:
  - `headline`
  - `body`
  - `caveat`
- no conditionals
- no formatting language
- no executable expressions

If omitted, clients may fall back to generic template-owned copy.

## `compatibility`

Portable input metadata, not a promise of complete support.

```yaml
compatibility:
  works_with:
    - "23andme"
    - "wgs"
    - "vcf"
  assemblies:
    - "grch37"
    - "grch38"
  notes: "Deletion-based calls may be absent from some chip exports."
```

### `works_with`

Human-readable source classes or file families such as:

- `23andme`
- `wgs`
- `vcf`
- `cram`

### `assemblies`

Known supported assemblies.

### `notes`

Plain-text compatibility caveats.

## `privacy`

Portable privacy contract for clients.

```yaml
privacy:
  mode: "local_only"
  uploads_data: false
  stores_results_locally: true
  external_urls: []
```

## Variant-Level Meaning

Do not duplicate variant semantics in `assay.yaml` when they belong in variant files.

Use variant YAML for:

- per-variant summary
- allele-specific findings
- source links

Use `assay.yaml` for:

- overall assay purpose
- category
- renderer type
- privacy
- compatibility
- assay-level interpretation framing

## Guidance For Clients

Clients should:

1. validate `schema`, `version`, `ui.template`, and `ui.version`
2. load assay metadata from `assay.yaml`
3. load assay members from `catalogue.yaml` and variant YAML files
4. render using the declared template family
5. use `interpretation` content when present, otherwise fall back to template defaults

Clients should not:

- ship assay-specific hardcoded metadata for package-authored assays
- treat YAML as a rendering DSL
- invent new category values that diverge from the shared set
