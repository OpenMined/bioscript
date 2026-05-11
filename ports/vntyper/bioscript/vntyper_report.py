"""HTML reporting helpers for the minimal VNtyper BioScript port."""

from __future__ import annotations

from html import escape


def render_html_report(report: dict) -> str:
    metadata = report.get("metadata", {})
    coverage = report.get("coverage", {})
    kestrel_rows = report.get("kestrel_variants", [])
    pipeline_log = report.get("pipeline_log", [])
    igv = report.get("igv", {})
    return "\n".join(
        [
            "<!doctype html>",
            '<html lang="en">',
            "<head>",
            '<meta charset="utf-8">',
            "<title>VNtyper BioScript Report</title>",
            _style(),
            _script(),
            "</head>",
            "<body>",
            "<main>",
            "<h1>VNtyper BioScript Report</h1>",
            _section("Screening Summary", f"<p>{_trusted_breaks(report.get('screening_summary', ''))}</p>"),
            _section("Run Metadata", _definition_list(metadata)),
            _details_section("VNTR Coverage QC", _definition_list(coverage), open_by_default=True),
            _section("Kestrel Identified Variants", _variant_table(kestrel_rows)),
            _section("IGV Visualization", _igv_section(igv, kestrel_rows)),
            _details_section("Pipeline Log", _log_block(pipeline_log), open_by_default=False),
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


def _details_section(title: str, body: str, open_by_default: bool = False) -> str:
    open_attr = " open" if open_by_default else ""
    return f"<section><details{open_attr}><summary>{escape(title)}</summary>{body}</details></section>"


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
    controls = (
        '<div class="table-tools">'
        '<label>Search <input id="variant-search" type="search" oninput="filterVariants()" '
        'placeholder="Filter variants"></label>'
        '<label><input id="show-flagged" type="checkbox" checked onchange="filterVariants()"> '
        "Show flagged rows</label>"
        "</div>"
    )
    header = "".join(
        f'<th><button type="button" onclick="sortVariants({idx})">{escape(column)}</button></th>'
        for idx, column in enumerate(columns)
    )
    body_rows = []
    for row in rows:
        flagged = row.get("Flag", "Not flagged") != "Not flagged"
        cells = "".join(_variant_cell(column, row.get(column, "")) for column in columns)
        body_rows.append(f'<tr data-flagged="{str(flagged).lower()}">{cells}</tr>')
    table = f'<table id="kestrel-table"><thead><tr>{header}</tr></thead><tbody>{"".join(body_rows)}</tbody></table>'
    return controls + table


def _variant_cell(column: str, value) -> str:
    content = escape(_display_value(value))
    if column == "Confidence":
        css = "confidence " + _confidence_class(str(value))
        return f'<td><span class="{css}">{content}</span></td>'
    if column == "Flag":
        flagged = str(value) not in ("", "None", "Not flagged", "Not applicable")
        icon = "!" if flagged else "-"
        title = "Flagged variant" if flagged else "Not flagged"
        return f'<td><span class="flag" title="{title}">{icon}</span> {content}</td>'
    return f"<td>{content}</td>"


def _confidence_class(value: str) -> str:
    normalized = value.lower().replace("*", "star").replace("_", "-")
    if "high-precision" in normalized:
        return "confidence-high"
    if "low-precision" in normalized:
        return "confidence-low"
    if "negative" in normalized:
        return "confidence-negative"
    return "confidence-other"


def _log_block(lines: list[str]) -> str:
    if not lines:
        return "<p>No pipeline log entries recorded.</p>"
    return "<pre>" + escape("\n".join(str(line) for line in lines)) + "</pre>"


def _igv_section(igv: dict, variants: list[dict]) -> str:
    if not igv:
        return "<p>IGV visualization is not configured for this report.</p>"
    required = ["reference", "bam", "vcf"]
    missing = [key for key in required if not igv.get(key)]
    if missing:
        return f"<p>IGV visualization is missing: {escape(', '.join(missing))}</p>"
    selector = _igv_variant_selector(variants)
    config = {
        "reference": igv["reference"],
        "bam": igv["bam"],
        "bai": igv.get("bai"),
        "vcf": igv["vcf"],
        "locus": igv.get("locus"),
    }
    attrs = " ".join(f'data-{key}="{escape(_display_value(value))}"' for key, value in config.items() if value)
    return (
        selector
        + f'<div id="igv-viewer" {attrs}></div>'
        + '<script src="https://cdn.jsdelivr.net/npm/igv@2.15.13/dist/igv.min.js"></script>'
        + _igv_script()
    )


def _igv_variant_selector(variants: list[dict]) -> str:
    if not variants:
        return "<p>No variants available for IGV selection.</p>"
    rows = []
    for row in variants:
        label = f"{row.get('CHROM', 'MUC1')}:{row.get('POS', '')} {row.get('REF', '')}>{row.get('ALT', '')}"
        locus = f"{row.get('CHROM', 'MUC1')}:{row.get('POS', '')}"
        rows.append(
            '<tr>'
            f"<td>{escape(label)}</td>"
            f'<td><button type="button" data-locus="{escape(locus)}" onclick="jumpIgv(this.dataset.locus)">View</button></td>'
            "</tr>"
        )
    return (
        '<table class="variant-selector"><thead><tr><th>Variant</th><th>IGV</th></tr></thead>'
        f"<tbody>{''.join(rows)}</tbody></table>"
    )


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
th button{border:0;background:transparent;font:inherit;font-weight:700;cursor:pointer}
.table-tools{display:flex;gap:18px;align-items:center;margin:0 0 12px}
.table-tools input[type=search]{padding:6px;border:1px solid #bbb}
.confidence{font-weight:700}
.confidence-high{color:#116329}
.confidence-low{color:#8a5a00}
.confidence-negative{color:#6b7280}
.flag{display:inline-block;min-width:16px;text-align:center;font-weight:700}
details summary{cursor:pointer;font-weight:700;font-size:18px;margin:0 0 12px}
pre{white-space:pre-wrap;background:#111;color:#f7f7f7;padding:12px;overflow:auto}
</style>"""


def _script() -> str:
    return """<script>
function filterVariants(){
  const q=document.getElementById('variant-search').value.toLowerCase();
  const showFlagged=document.getElementById('show-flagged').checked;
  document.querySelectorAll('#kestrel-table tbody tr').forEach(row=>{
    const flagged=row.dataset.flagged==='true';
    const match=row.innerText.toLowerCase().includes(q);
    row.style.display=(match && (showFlagged || !flagged)) ? '' : 'none';
  });
}
function sortVariants(index){
  const tbody=document.querySelector('#kestrel-table tbody');
  const rows=Array.from(tbody.querySelectorAll('tr'));
  const numeric=value=>value!=='' && !Number.isNaN(Number(value));
  rows.sort((a,b)=>{
    const av=a.children[index].innerText.trim();
    const bv=b.children[index].innerText.trim();
    if(numeric(av) && numeric(bv)){return Number(av)-Number(bv);}
    return av.localeCompare(bv);
  });
  rows.forEach(row=>tbody.appendChild(row));
}
</script>"""


def _igv_script() -> str:
    return """<script>
let bioscriptIgvBrowser=null;
function initBioScriptIgv(){
  const el=document.getElementById('igv-viewer');
  if(!el || !window.igv){return;}
  const options={
    genome:{fastaURL:el.dataset.reference},
    locus:el.dataset.locus || undefined,
    tracks:[
      {name:'BAM',type:'alignment',format:'bam',url:el.dataset.bam,indexURL:el.dataset.bai},
      {name:'Kestrel VCF',type:'variant',format:'vcf',url:el.dataset.vcf}
    ]
  };
  igv.createBrowser(el, options).then(browser=>{bioscriptIgvBrowser=browser;});
}
function jumpIgv(locus){
  if(bioscriptIgvBrowser){bioscriptIgvBrowser.search(locus);}
}
window.addEventListener('load', initBioScriptIgv);
</script>"""
