"""HTML reporting helpers for the minimal VNtyper BioScript port."""

from __future__ import annotations

from html import escape


def render_html_report(report: dict) -> str:
    metadata = report.get("metadata", {})
    coverage = report.get("coverage", {})
    kestrel_rows = report.get("kestrel_variants", [])
    pipeline_log = report.get("pipeline_log", [])
    return "\n".join(
        [
            "<!doctype html>",
            '<html lang="en">',
            "<head>",
            '<meta charset="utf-8">',
            "<title>VNtyper BioScript Report</title>",
            _style(),
            "</head>",
            "<body>",
            "<main>",
            "<h1>VNtyper BioScript Report</h1>",
            _section("Screening Summary", f"<p>{_trusted_breaks(report.get('screening_summary', ''))}</p>"),
            _section("Run Metadata", _definition_list(metadata)),
            _section("VNTR Coverage QC", _definition_list(coverage)),
            _section("Kestrel Identified Variants", _variant_table(kestrel_rows)),
            _section("Pipeline Log", _log_block(pipeline_log)),
            "</main>",
            "</body>",
            "</html>",
        ]
    )


def write_html_report(path: str, report: dict) -> None:
    with open(path, "w", encoding="utf-8") as handle:
        handle.write(render_html_report(report))


def _section(title: str, body: str) -> str:
    return f"<section><h2>{escape(title)}</h2>{body}</section>"


def _definition_list(values: dict) -> str:
    if not values:
        return "<p>Not available</p>"
    rows = []
    for key, value in values.items():
        rows.append(f"<dt>{escape(str(key))}</dt><dd>{escape(_display_value(value))}</dd>")
    return "<dl>" + "".join(rows) + "</dl>"


def _variant_table(rows: list[dict]) -> str:
    columns = [
        "Motif",
        "Variant",
        "POS",
        "REF",
        "ALT",
        "Motif_sequence",
        "Estimated_Depth_AlternateVariant",
        "Estimated_Depth_Variant_ActiveRegion",
        "Depth_Score",
        "Confidence",
        "Flag",
    ]
    if not rows:
        return "<p>No Kestrel variants reported.</p>"
    header = "".join(f"<th>{escape(column)}</th>" for column in columns)
    body_rows = []
    for row in rows:
        cells = "".join(f"<td>{escape(_display_value(row.get(column, '')))}</td>" for column in columns)
        body_rows.append(f"<tr>{cells}</tr>")
    return f"<table><thead><tr>{header}</tr></thead><tbody>{''.join(body_rows)}</tbody></table>"


def _log_block(lines: list[str]) -> str:
    if not lines:
        return "<p>No pipeline log entries recorded.</p>"
    return "<pre>" + escape("\n".join(str(line) for line in lines)) + "</pre>"


def _trusted_breaks(value: str) -> str:
    return escape(str(value)).replace("&lt;br&gt;", "<br>")


def _display_value(value) -> str:
    if value is None:
        return "Not available"
    if isinstance(value, list):
        return ", ".join(str(item) for item in value) if value else "None"
    return str(value)


def _style() -> str:
    return """<style>
body{font-family:Arial,sans-serif;margin:0;background:#f6f7f9;color:#222}
main{max-width:1120px;margin:0 auto;padding:32px}
section{margin:20px 0;padding:18px;background:#fff;border:1px solid #ddd}
h1{font-size:28px;margin:0 0 24px}
h2{font-size:18px;margin:0 0 12px}
dl{display:grid;grid-template-columns:220px 1fr;gap:8px 16px;margin:0}
dt{font-weight:700}
dd{margin:0}
table{border-collapse:collapse;width:100%;font-size:13px}
th,td{border:1px solid #ddd;padding:6px;text-align:left}
th{background:#eef1f5}
pre{white-space:pre-wrap;background:#111;color:#f7f7f7;padding:12px;overflow:auto}
</style>"""
