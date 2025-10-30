from __future__ import annotations

import csv
import sys
from collections.abc import Iterator
from typing import Iterable

from .types import VariantRow

REQUIRED = ("rsid", "chromosome", "position", "genotype")
COMMENT_PREFIXES = ("#", "//")


def _normalize_name(name: str) -> str:
    """
    Normalize header/column names for comparison by removing whitespace,
    hyphens, and underscores and lowercasing the result.
    """
    return "".join(ch for ch in name.strip().lower() if ch not in " _-")


FIELD_ALIASES = {
    "rsid": {"rsid", "name", "snp", "marker", "id"},
    "chromosome": {"chromosome", "chr", "chrom"},
    "position": {"position", "pos", "coordinate", "basepairposition", "basepair"},
    "genotype": {
        "genotype",
        "gt",
        "result",
        "results",
        "result1",
        "call",
        "calls",
        "yourcode",
        "code",
        "genotypevalue",
        "variation",
    },
    "allele1": {"allele1", "allelea", "allele_a", "allele1top"},
    "allele2": {"allele2", "alleleb", "allele_b", "allele2top"},
    "gs": {"gs"},
    "baf": {"baf"},
    "lrr": {"lrr"},
}


def _float_or_none(x: str | None) -> float | None:
    if x is None:
        return None
    x = x.strip()
    if not x:
        return None
    try:
        return float(x)
    except ValueError:
        return None


def _int_or_none(x: str | None) -> int | None:
    if x is None:
        return None
    x = x.strip()
    if not x:
        return None
    try:
        return int(x)
    except ValueError:
        return None


def _detect_delimiter(lines: Iterable[str]) -> str:
    """
    Detect the delimiter used in the genotype file by inspecting the first non-comment line.
    """
    for raw in lines:
        stripped = raw.strip()
        if not stripped or stripped.startswith(COMMENT_PREFIXES):
            continue
        if "\t" in raw:
            return "\t"
        if "," in raw:
            return ","
        if len(stripped.split()) > 1:
            return " "
    return "\t"


def _parse_fields(text: str, delimiter: str) -> list[str]:
    """Parse a single line into fields while respecting quoting."""
    reader = csv.reader([text], delimiter=delimiter, skipinitialspace=True)
    fields = [field.strip() for field in next(reader)]
    if len(fields) == 1 and delimiter in {"\t", " "}:
        split_fields = text.strip().split()
        if len(split_fields) > 1:
            return split_fields
    return fields


def _looks_like_header(fields: list[str]) -> bool:
    if not fields:
        return False
    first = _normalize_name(fields[0])
    return first in FIELD_ALIASES["rsid"]


def _default_header(field_count: int) -> list[str]:
    base = list(REQUIRED)
    if field_count <= len(base):
        return base[:field_count]
    extras = [f"extra_{i}" for i in range(field_count - len(base))]
    return base + extras


def _strip_inline_comment(value: str) -> str:
    """
    Remove inline comments starting with '#' or '//' from a field value.
    """
    for marker in ("#", "//"):
        idx = value.find(marker)
        if idx != -1:
            return value[:idx].strip()
    return value.strip()


def _normalize_genotype_value(value: str) -> str:
    """
    Normalize genotype strings, accepting formats like 'A/G' or 'TTATAA/-'.
    Treat '#N/A' and empty values as '--'.
    Return 'ID' when a deletion marker ('-') appears alongside another allele.
    """
    cleaned = (value or "").strip()
    if not cleaned or cleaned.upper() in {"NA", "N/A", "#N/A", "NONE"}:
        return "--"

    cleaned = cleaned.replace(" ", "").upper()
    if "/" in cleaned:
        parts = [part.strip() for part in cleaned.split("/") if part is not None]
        if any(part in {"", "-"} for part in parts):
            return "ID"
        cleaned = "".join(parts)

    return cleaned.replace("/", "")


def _extract_header_and_data(lines: list[str]) -> tuple[list[str], list[tuple[str, list[str]]]]:
    """
    Returns (header_fields, data_entries) where each entry is (raw_line, parsed_fields).
    Supports commented headers and infers a default header when none is present.
    """
    normalized_lines = [ln.rstrip("\n\r") for ln in lines]
    delimiter = _detect_delimiter(normalized_lines)

    header_fields: list[str] | None = None
    comment_header: list[str] | None = None
    data_entries: list[tuple[str, list[str]]] = []

    for line in normalized_lines:
        stripped = line.strip()
        if not stripped:
            continue

        if stripped.startswith(COMMENT_PREFIXES):
            candidate = stripped.lstrip("#/").strip()
            if not candidate:
                continue
            fields = _parse_fields(candidate, delimiter)
            if _looks_like_header(fields):
                comment_header = [field.strip() for field in fields]
            continue

        fields = _parse_fields(line, delimiter)

        if header_fields is None:
            if _looks_like_header(fields):
                header_fields = [field.strip() for field in fields]
                continue
            if comment_header is not None:
                header_fields = comment_header
                comment_header = None
            else:
                header_fields = _default_header(len(fields))

        data_entries.append((line, fields))

    if header_fields is None:
        raise ValueError("No header found or inferred in genotype file.")

    return header_fields, data_entries


def load_variants_tsv(path: str) -> Iterator[VariantRow]:
    """
    Load a genotype file (TSV/CSV) and yield VariantRow objects.
    Supports:
      - plain TSV files with rsid/chromosome/position/genotype columns
      - CSV files (e.g. FamilyTreeDNA, MyHeritage, deCODEme)
      - Commented headers (e.g. 23andMe style)
      - Split allele columns (allele1/allele2) merged into genotype
      - Optional columns: gs, baf, lrr

    Args:
        path: Path to genotype file
    """
    with open(path, encoding="utf-8-sig", newline="") as f:
        all_lines = f.readlines()

    header, data_entries = _extract_header_and_data(all_lines)
    normalized_header = [_normalize_name(h) for h in header]
    column_keys = dict(zip(header, normalized_header))

    def _lookup(row_map: dict[str, str], key: str) -> str | None:
        aliases = FIELD_ALIASES.get(key, set())
        for candidate in aliases:
            value = row_map.get(candidate)
            if value is not None and value != "":
                return value
        return None

    available = set(normalized_header)
    required_missing = []
    for required in ("rsid", "chromosome", "position"):
        if available.isdisjoint(FIELD_ALIASES[required]):
            required_missing.append(required)
    genotype_available = not available.isdisjoint(FIELD_ALIASES["genotype"])
    allele_split_available = not available.isdisjoint(
        FIELD_ALIASES["allele1"]
    ) and not available.isdisjoint(FIELD_ALIASES["allele2"])
    if not genotype_available and not allele_split_available:
        required_missing.append("genotype (or allele1/allele2)")
    if required_missing:
        raise ValueError(f"Missing required column(s): {', '.join(required_missing)}")

    for raw_line, values in data_entries:
        row_map: dict[str, str] = {}
        for idx, value in enumerate(values):
            if idx >= len(header):
                continue
            key = column_keys[header[idx]]
            row_map[key] = _strip_inline_comment(value)

        rsid = (_lookup(row_map, "rsid") or "").strip()
        chrom = (_lookup(row_map, "chromosome") or "").strip()
        pos = _int_or_none(_lookup(row_map, "position"))
        genotype = _lookup(row_map, "genotype")

        if genotype is None:
            allele1 = _lookup(row_map, "allele1") or ""
            allele2 = _lookup(row_map, "allele2") or ""
            genotype = (allele1 + allele2).strip()

        genotype = _normalize_genotype_value(genotype)

        if not rsid or not chrom or pos is None:
            print(
                f"[bioscript] Skipping row with missing required fields: {row_map}",
                file=sys.stderr,
            )
            continue

        if not genotype:
            genotype = "--"

        gs = _float_or_none(_lookup(row_map, "gs"))
        baf = _float_or_none(_lookup(row_map, "baf"))
        lrr = _float_or_none(_lookup(row_map, "lrr"))

        yield VariantRow(
            rsid=rsid,
            chromosome=chrom,
            position=pos,
            genotype=genotype,
            assembly=None,  # Always None - user must set when constructing VariantRow
            gs=gs,
            baf=baf,
            lrr=lrr,
            raw_line=raw_line,
        )
