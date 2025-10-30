from __future__ import annotations

from pathlib import Path

import pytest

from bioscript import load_variants_tsv

REPO_ROOT = Path(__file__).resolve().parents[3]
GENOTYPE_DIR = REPO_ROOT / "cli/tests/data/genotype_files"

TARGET_GENOTYPES = {
    "rs73885319": "AG",
    "rs60910145": "TG",
    "rs71785313": "ID",
    "rs143830837": "ID",
}

SUPPORTED_FILES = [
    "01_23andme_grch36_standard_4col.txt",
    "02_23andme_grch37_standard_4col.txt",
    "03_23andme_standard_4col.txt",
    "04_ancestrydna_grch37_v1_split_alleles_5col.txt",
    "05_ancestrydna_split_alleles_5col.txt",
    "06_decodeme_grch36_variation_6col.csv",
    "07_dynamicdna_grch38_extended_7col.txt",
    "08_familytreedna_grch37_standard_4col.csv",
    "09_genesforgood_standard_4col.txt",
    "10_livingdna_grch37_standard_4col.txt",
    "11_myheritage_grch37_standard_4col.csv",
    "12_unknown_no_header.txt",
    "13_unknown_standard_4col.txt",
]


@pytest.mark.parametrize("file_name", SUPPORTED_FILES)
def test_load_variants_handles_supported_formats(file_name: str) -> None:
    path = GENOTYPE_DIR / file_name
    rows = list(load_variants_tsv(str(path)))
    by_rsid = {row.rsid: row for row in rows}

    for rsid, expected in TARGET_GENOTYPES.items():
        assert rsid in by_rsid, f"{file_name}: missing {rsid}"
        assert by_rsid[rsid].genotype == expected, (
            f"{file_name}: expected {rsid}={expected}, got {by_rsid[rsid].genotype}"
        )


def test_dynamicdna_optional_fields_preserved() -> None:
    path = GENOTYPE_DIR / "07_dynamicdna_grch38_extended_7col.txt"
    rows = list(load_variants_tsv(str(path)))
    target = next(row for row in rows if row.rsid == "rs73885319")

    assert pytest.approx(target.gs, rel=1e-4) == 0.7421
    assert pytest.approx(target.baf, rel=1e-4) == 0.5217
    assert pytest.approx(target.lrr, rel=1e-4) == -0.0843
