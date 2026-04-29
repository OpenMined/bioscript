# Report Schema

This schema defines the normalized BioScript report object derived from one or more `bioscript:observation:1.0` rows.

The report answers: "What does this assay conclude for this participant?"

It is intentionally separate from observations:

- observations are row-based evidence for individual variants or sites
- reports are derived assay-level outputs for one participant and one assay
- rendered Markdown/HTML reports should be generated from this object plus the underlying observations

## Schema Identity

```yaml
schema: "bioscript:report:1.0"
version: "1.0"
```

## Minimal Example

```json
{
  "participant_id": "NA06985",
  "assay_id": "APOL1",
  "assay_version": "1.0",
  "report_status": "complete",
  "derived_from": ["G1_SITE_1", "G1_SITE_2", "G2_SITE"],
  "facets": "self_reported_ethnicity=African Caribbean;cohort=CariGenetics",
  "timing": {
    "started_at": "2026-04-29T10:00:00Z",
    "finished_at": "2026-04-29T10:00:00.120Z",
    "duration_ms": 120,
    "breakdown": {
      "load_ms": 30,
      "query_ms": 40,
      "compute_ms": 20,
      "write_ms": 30
    }
  },
  "fields": [
    {
      "key": "apol1_status",
      "label": "APOL1 status",
      "value": "G1/G0",
      "value_type": "string",
      "format": "badge"
    },
    {
      "key": "summary",
      "label": "Summary",
      "value": "APOL1 genotype pattern is G1/G0.",
      "value_type": "string",
      "format": "plain_text"
    }
  ],
  "metrics": {
    "n_variants": 2,
    "n_sites_tested": 3
  }
}
```

## Fields

Required identity fields:

- `participant_id`: participant/sample identifier.
- `assay_id`: stable assay or panel identifier.
- `assay_version`: assay or panel version.

Required report status:

- `report_status`: top-level completion state.

Provenance over observations:

- `derived_from`: list of observation `variant_key` values or stable observation IDs used to produce this report.
- `facets`: optional semicolon-separated key-value metadata copied from the run, participant, cohort, or observations.

Timing:

- `timing.started_at`: ISO 8601 timestamp.
- `timing.finished_at`: ISO 8601 timestamp.
- `timing.duration_ms`: total report-generation duration.
- `timing.breakdown.load_ms`: input loading time, or `null`.
- `timing.breakdown.query_ms`: observation lookup/query time, or `null`.
- `timing.breakdown.compute_ms`: assay-specific compute/interpretation time, or `null`.
- `timing.breakdown.write_ms`: output writing time, or `null`.

Report content:

- `fields`: ordered list of typed report fields.
- `metrics`: optional object for machine-readable aggregate values.
- `warnings`: optional list of warning strings.
- `notes`: optional list of note strings.

## Controlled Vocabularies

### `report_status`

Allowed values:

- `complete`
- `partial`
- `failed`

Use `complete` when all required observations were available and the assay produced its intended report fields.

Use `partial` when the report is usable but one or more observations, calls, or optional outputs were unavailable.

Use `failed` when the assay could not produce a trustworthy report object.

### `fields[*].value_type`

Recommended values:

- `string`
- `number`
- `boolean`
- `enum`
- `markdown`
- `json`

`value_type` describes the type of `value`, not just the renderer. For example, an APOL1 `G1/G0` status is `string` or `enum`; a count is `number`; a structured object is `json`.

### `fields[*].format`

Recommended values:

- `plain_text`
- `badge`
- `markdown`
- `table`
- `json`

`format` is a rendering hint. It should not change the meaning of `value`.

Examples:

```json
{
  "key": "apol1_status",
  "label": "APOL1 status",
  "value": "G1/G0",
  "value_type": "string",
  "format": "badge"
}
```

```json
{
  "key": "summary",
  "label": "Summary",
  "value": "APOL1 genotype pattern is G1/G0.",
  "value_type": "string",
  "format": "plain_text"
}
```

## Metrics

`metrics` is an optional object for aggregate values that are useful for filtering, quality checks, dashboards, or cohort summaries.

Examples:

```json
{
  "n_variants": 2,
  "n_sites_tested": 3,
  "n_sites_missing": 0
}
```

Metrics should be machine-oriented. Human prose belongs in `fields`, `warnings`, or `notes`.

## Facets

`facets` uses the same format as observation facets:

```text
key=value;key=value
```

Examples:

```text
self_reported_ethnicity=African Caribbean;cohort=CariGenetics
cohort=TestBatch1
```

Rules:

- keys should use lowercase snake_case
- values should not contain tabs or newlines
- semicolons in values should be percent-encoded as `%3B`
- equals signs in values should be percent-encoded as `%3D`

## Relationship To Observations

The report is derived from observation rows.

For APOL1:

- observations record `G1_SITE_1`, `G1_SITE_2`, and `G2_SITE`
- the report's `derived_from` lists those observation keys
- the report's `fields` include `apol1_status`
- the report's `metrics` can count tested sites, called sites, missing sites, or variant sites

This gives downstream systems two clean query levels:

- query observations for variant evidence and cohort-level allele counts
- query reports for assay-level participant results

## Storage Guidance

Recommended storage layout:

- `observations.tsv` or `observations` table for `bioscript:observation:1.0`
- `reports.jsonl` or `reports` table for `bioscript:report:1.0`
- `report_fields` table when field-level SQL querying is needed

A relational split can look like:

```text
reports(
  participant_id,
  assay_id,
  assay_version,
  report_status,
  derived_from,
  facets,
  started_at,
  finished_at,
  duration_ms,
  metrics_json,
  warnings_json,
  notes_json
)

report_fields(
  participant_id,
  assay_id,
  assay_version,
  key,
  label,
  value_json,
  value_type,
  format
)
```

The JSON form remains the canonical exchange format.
