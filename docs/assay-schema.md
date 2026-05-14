# Assay Schema

Use an assay when a named test observes one or more variants and emits custom derived report fields.

An assay is different from a panel: a panel is a collection of mostly independent observations, while an assay has its own interpretation logic. APOL1 is an assay because it observes G1/G2 sites and reports one derived APOL1 status.

## Schema Identity

```yaml
schema: "bioscript:assay:1.0"
version: "1.0"
```

## Minimal Shape

```yaml
schema: "bioscript:assay:1.0"
version: "1.0"
name: "APOL1"
label: "APOL1 Risk Assay"
tags:
  - "type:risk"
  - "gene:APOL1"

members:
  - kind: "variant"
    path: "g1-site-1.yaml"
    version: "1.0"
  - kind: "variant"
    path: "g1-site-2.yaml"
    version: "1.0"
  - kind: "variant"
    path: "g2-site.yaml"
    version: "1.0"

analyses:
  - id: "apol1_status"
    kind: "bioscript"
    path: "apol1.py"
    output_format: "tsv"
    label: "APOL1 risk genotype"
    derived_from:
      - "g1-site-1.yaml"
      - "g1-site-2.yaml"
      - "g2-site.yaml"
    assets:
      - id: "variants"
        path: "variants.tsv"
      - id: "findings"
        path: "findings.tsv"
    emits:
      - key: "apol1_status"
        label: "APOL1 status"
        value_type: "string"
        format: "badge"
    logic:
      source:
        name: "Example derivation source"
        url: "https://example.org/assay-logic"
      description: >
        Optional human-readable description of the derivation logic implemented by the analysis script.
```

## Members

Assay members are currently local variant YAML files:

```yaml
- kind: "variant"
  path: "g1-site-1.yaml"
  version: "1.0"
```

Rules:

- `kind` is required and currently must be `variant`
- `path` is required
- `version` is recommended
- keep variant identity, coordinates, alleles, findings, and provenance in the variant YAML files

## Analyses

Use `analyses` for custom output derived from the member variants. The older `interpretations` key is accepted for compatibility, but new manifests should use `analyses`.

Rules:

- `id`, `kind`, `path`, and `derived_from` are required
- `kind` is currently `bioscript`
- `path` points to a BioScript-compatible Python file
- `output_format` is optional and defaults to `tsv`; use `json` or `jsonl` when the script writes structured JSON output
- `derived_from` lists the variant YAML files used by the interpretation
- `assets` is optional and lists local files the analysis script can read through the injected `asset_paths` dict
- `emits` is optional but recommended so report generators know which output columns to display and how to label them
- `logic` is optional; use `logic.description` and `logic.source.url` to document where the script's derivation rules came from
- Analysis rows may emit `notes` or `report_notes` as a reporting convention. HTML reports render those notes below the analysis table and omit them from the table columns; this avoids a manifest-level template language while still letting the script build human-readable text from computed values.

Analysis scripts receive these injected variables:

- `input_file`: input genotype/VCF/CRAM path
- `output_file`: path the script must write
- `participant_id`: current participant/sample identifier
- `observations_file`: TSV containing the variant observations already gathered from assay/panel members for this participant
- `asset_paths`: dict from each `assets[].id` to the resolved local path

This supports large catalogues where Rust/BioScript first gathers all variant observations, then Python joins those observations to attached TSV assets such as `findings.tsv`, `conditions.tsv`, or `rules.tsv` and emits derived classification rows.

## Findings

Use `findings` for evidence that binds either to a variant observation or an emitted analysis value. Keep the executable logic in `analyses`; keep PGx evidence and reporting semantics in YAML.

```yaml
findings:
  - schema: "bioscript:pgx-label:1.0"
    id: "clinpgx_PA166313401"
    label: "ClinPGx drug label annotation PA166313401"
    authority_type: "regulatory_label"
    binding:
      source: "analysis"
      analysis_id: "apoe_epsilon"
      key: "apoe_status"
      operator: "equals"
      value: "e4/e4"
    drugs:
      - name: "lecanemab"
        aliases:
          - "LEQEMBI"
    evidence:
      source: "ClinPGx"
      kind: "label_annotation"
      id: "PA166313401"
      url: "https://www.clinpgx.org/labelAnnotation/PA166313401"
    notes: "Drug label annotation applies when APOE status is e4/e4."
```

Binding rules:

- `source` is `analysis` or `variant`
- `analysis` bindings require `analysis_id`, `key`, and either `operator: equals` with `value` or `operator: in` with `values`
- `variant` bindings require `variant` or `path`, `key`, and either `equals`/`value` or `in`/`values`
- PGx label findings use `schema: "bioscript:pgx-label:1.0"` and should include `regulatory_sources`, `pgx_action_level` or `prescribing_actions` when known
- PGx summary findings use `schema: "bioscript:pgx-summary:1.0"` and should include `evidence_level`, `phenotype_categories`, and genotype-specific `effects` when known
- PGx findings should include `drugs` and should link to the exact ClinPGx/PharmGKB/ClinVar evidence page

## Inclusion In Panels

A larger panel may include an assay as a member:

```yaml
members:
  - kind: "assay"
    path: "../risk/APOL1/assay.yaml"
    version: "1.0"
```

When a panel includes an assay, the assay's variant observations can be expanded into the panel output, while report tooling can also run the assay's interpretation and include its emitted fields.
