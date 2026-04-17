# Panel Schema

Use a panel when you want one manifest that points to a curated set of runnable variant records.

Right now the Rust runner supports variant members directly. Keep the shape simple.

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
  - kind: "variant"
    path: "variants/rs713598.yaml"
    version: "1.0"
```

## Purpose

A panel is:

- a selection manifest
- a stable name for a bundle of variants
- something the Rust `bioscript` command can run directly

It is not:

- a full remote package manager
- a replacement for richer assay manifests

## Members

Each member must currently be:

```yaml
- kind: "variant"
  path: "variants/rs671.yaml"
  version: "1.0"
```

Rules:

- `kind` is required
- exactly one of `path` or `download` is required
- current runner support is `variant` members only
- `version` is recommended for local members
- `sha256` is optional for local members

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
