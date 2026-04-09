# Variant Catalogue Schema

Use a catalogue when you want to manage many `bioscript:variant` files together as a panel.

This keeps each variant small and self-contained while giving you one manifest for:

- grouping
- version pinning
- caching
- batch validation
- panel execution

## Recommended Shape

```yaml
schema: "bioscript:catalogue"
version: "1.0"
catalogue_id: "apol1-panel"
label: "APOL1 panel"

variants:
  - id: "APOL1_G1_rs73885319"
    path: "variants/apol1/rs73885319.yaml"
    version: "1.0"
  - id: "APOL1_G1_rs60910145"
    path: "variants/apol1/rs60910145.yaml"
    version: "1.0"
  - id: "APOL1_G2"
    path: "variants/apol1/g2.yaml"
    version: "1.0"
```

## Why This Shape

- `path` is easy to load locally
- `id` lets you refer to a variant without reparsing the path
- `version` lets you cache and pin a known record revision

## Optional Extensions

Later you can add:

- `sha256`
- `source_url`
- `tags`
- `groups`
- `assemblies`

Example:

```yaml
schema: "bioscript:catalogue"
version: "1.0"
catalogue_id: "pgx-core"

variants:
  - id: "ABCG2_rs2231142"
    path: "variants/abcg2/rs2231142.yaml"
    version: "1.0"
    tags:
      - "pgx"
      - "statin"
```

## Recommended Workflow

1. Keep one variant per file.
2. Keep those files in a predictable directory layout.
3. Use a catalogue file to define the panel.
4. Validate both the catalogue and each referenced variant.
5. Cache by `id + version`.
