# Panel Schema

Use a panel when you want one manifest that points to a curated set of runnable variant records, assay manifests, and optional interpretation scripts derived from those records.

The Rust runner supports variant members directly. Test tooling can also run declared interpretation scripts and add their emitted fields to the generated report.

## Schema Identity

```yaml
schema: "bioscript:panel:1.0"
version: "1.0"
```

## Minimal Shape

```yaml
schema: "bioscript:panel:1.0"
version: "1.0"
name: "traits-common"
label: "Common Traits"
tags:
  - "type:trait"

members:
  - kind: "variant"
    path: "variants/rs671.yaml"
    version: "1.0"
  - kind: "assay"
    path: "../risk/APOL1/assay.yaml"
    version: "1.0"
  - kind: "variant"
    path: "variants/rs713598.yaml"
    version: "1.0"

analyses:
  - id: "taste_status"
    kind: "bioscript"
    path: "interpretations/taste.py"
    output_format: "tsv"
    label: "Taste status"
    derived_from:
      - "variants/rs713598.yaml"
    emits:
      - key: "taste_status"
        label: "Taste status"
        value_type: "string"
        format: "badge"
```

## Purpose

A panel is:

- a selection manifest
- a stable name for a bundle of variants
- something the Rust `bioscript` command can run directly
- a way to include smaller assay manifests in a broader bundle
- a place to declare interpretation chunks that derive custom report fields from member variants

It is not:

- a full remote package manager
- a replacement for richer assay manifests
- a place to hide variant metadata inside Python when YAML can describe it

## Members

Each member must currently be a local variant or assay:

```yaml
- kind: "variant"
  path: "variants/rs671.yaml"
  version: "1.0"
- kind: "assay"
  path: "../risk/APOL1/assay.yaml"
  version: "1.0"
```

Rules:

- `kind` is required
- exactly one of `path` or `download` is required
- current runner support is local `variant` and `assay` members
- `version` is recommended for local members
- `sha256` is optional for local members

## Analyses

Use `analyses` when a panel needs custom derived output that is not the same thing as a single variant observation. Examples include APOE epsilon genotype from rs429358/rs7412 or APOL1 G0/G1/G2 status from three sites. The older `interpretations` key is accepted for compatibility, but new manifests should use `analyses`.

```yaml
analyses:
  - id: "apoe_epsilon"
    kind: "bioscript"
    path: "variants/APOE/apoe.py"
    output_format: "tsv"
    label: "APOE epsilon genotype"
    derived_from:
      - "variants/APOE/rs429358.yaml"
      - "variants/APOE/rs7412.yaml"
    emits:
      - key: "apoe_status"
        label: "APOE status"
        value_type: "string"
        format: "badge"
    logic:
      source:
        name: "ClinPGx / PharmGKB"
        url: "https://www.clinpgx.org/variant/PA166155341/overview"
      description: >
        Optional human-readable description of the derivation logic implemented by the analysis script.
```

Rules:

- `id`, `kind`, `path`, and `derived_from` are required
- `kind` is currently `bioscript`
- `path` points to a BioScript-compatible Python file
- `output_format` is optional and defaults to `tsv`; use `json` or `jsonl` when the script writes structured JSON output
- `derived_from` lists the variant YAML files used by the interpretation
- `emits` is optional but recommended so report generators know which output columns to display and how to label them
- `logic` is optional; use `logic.description` and `logic.source.url` to document where the script's derivation rules came from
- keep variant identity, coordinates, alleles, findings, and provenance in YAML; keep cross-variant logic in the interpretation script

## Permissions And Downloads

Panels may declare remote downloads up front even if the current runner only executes local members.

```yaml
permissions:
  domains:
    - "https://example.org"

downloads:
  - id: "remote-rs671"
    url: "https://example.org/variants/rs671.yaml"
    sha256: "0123456789abcdef0123456789abcdef0123456789abcdef0123456789abcdef"
    version: "1.0"
```

Validation rules:

- `permissions.domains` entries must be origins only
- every `downloads[*].url` origin must also appear in `permissions.domains`
- `downloads[*].sha256` must be a 64-character lowercase hex digest

This keeps host approval focused on which remote origins may be contacted.

## Running Panels

Examples:

```bash
bioscript panel.yaml --input-file sample.txt --output-file output.tsv
bioscript panel.yaml --input-file sample.txt --filter tag=type:trait
bioscript panel.yaml --input-file sample.txt --filter name=rs671
```

Current filter keys:

- `kind`
- `name`
- `path`
- `tag`

## Future Relationship To Catalogues And Assays

- `panel` is the small runnable collection manifest
- `catalogue` can later become a larger published index over many panels and items
- `assay` can later become the richer multi-file runnable bundle with declared assets and overrides
