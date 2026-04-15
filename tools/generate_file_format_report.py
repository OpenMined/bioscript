#!/usr/bin/env python3
from __future__ import annotations

import argparse
import html
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path


INSPECT_EXTENSIONS = (
    ".zip",
    ".vcf",
    ".vcf.gz",
    ".cram",
    ".bam",
    ".fa",
    ".fasta",
    ".fa.gz",
    ".fna",
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--bioscript", required=True, type=Path)
    parser.add_argument("--root", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    return parser.parse_args()


def iter_targets(test_data_root: Path) -> list[Path]:
    targets: list[Path] = []
    for path in sorted(test_data_root.rglob("*")):
        if not path.is_file():
            continue
        lower = path.name.lower()
        if lower.endswith(INSPECT_EXTENSIONS):
            targets.append(path)
    return targets


def inspect_file(bioscript: Path, root: Path, path: Path) -> dict[str, str]:
    proc = subprocess.run(
        [str(bioscript), "inspect", str(path)],
        cwd=root,
        text=True,
        capture_output=True,
        check=False,
    )
    result = {
        "path": str(path.relative_to(root)),
        "status": "ok" if proc.returncode == 0 else "error",
        "stderr": proc.stderr.strip(),
    }
    if proc.returncode != 0:
        return result
    for line in proc.stdout.splitlines():
        if "\t" not in line:
            continue
        key, value = line.split("\t", 1)
        result[key] = value
    return result


def cell(value: str) -> str:
    return html.escape(value or "")


def render_html(rows: list[dict[str, str]]) -> str:
    generated_at = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M:%S UTC")
    total = len(rows)
    failures = sum(1 for row in rows if row["status"] != "ok")
    headers = [
        "path",
        "container",
        "kind",
        "confidence",
        "vendor",
        "platform_version",
        "assembly",
        "phased",
        "selected_entry",
        "has_index",
        "index_path",
        "reference_matches",
        "source_confidence",
        "duration_ms",
        "source_evidence",
        "evidence",
        "warnings",
        "stderr",
    ]
    lines = [
        "<!doctype html>",
        "<html lang='en'>",
        "<head>",
        "<meta charset='utf-8'>",
        "<title>Bioscript File Format Report</title>",
        "<style>",
        "body { font-family: ui-sans-serif, system-ui, sans-serif; margin: 24px; color: #1f2937; }",
        "h1 { margin-bottom: 8px; }",
        ".meta { margin-bottom: 20px; color: #4b5563; }",
        "table { border-collapse: collapse; width: 100%; font-size: 14px; }",
        "th, td { border: 1px solid #d1d5db; padding: 8px; text-align: left; vertical-align: top; }",
        "th { position: sticky; top: 0; background: #f3f4f6; }",
        "tr:nth-child(even) { background: #f9fafb; }",
        ".ok { background: #ecfdf5; }",
        ".error { background: #fef2f2; }",
        "code { font-family: ui-monospace, SFMono-Regular, monospace; font-size: 12px; }",
        "</style>",
        "</head>",
        "<body>",
        "<h1>Bioscript File Format Report</h1>",
        f"<div class='meta'>Generated {cell(generated_at)}. Files inspected: {total}. Failures: {failures}.</div>",
        "<table>",
        "<thead><tr>" + "".join(f"<th>{cell(header)}</th>" for header in headers) + "</tr></thead>",
        "<tbody>",
    ]
    for row in rows:
        css = "ok" if row["status"] == "ok" else "error"
        lines.append(f"<tr class='{css}'>")
        for header in headers:
            value = row.get(header, "")
            if header == "path":
                lines.append(f"<td><code>{cell(value)}</code></td>")
            else:
                lines.append(f"<td>{cell(value)}</td>")
        lines.append("</tr>")
    lines.extend(["</tbody>", "</table>", "</body>", "</html>"])
    return "\n".join(lines)


def main() -> int:
    args = parse_args()
    test_data_root = args.root / "test-data"
    if not args.bioscript.exists():
        print(f"bioscript binary not found: {args.bioscript}", file=sys.stderr)
        return 1
    if not test_data_root.exists():
        print(f"test-data directory not found: {test_data_root}", file=sys.stderr)
        return 1

    rows = [inspect_file(args.bioscript, args.root, path) for path in iter_targets(test_data_root)]
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(render_html(rows), encoding="utf-8")
    print(args.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
