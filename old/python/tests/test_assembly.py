"""Tests for genome assembly version support."""

import tempfile
from pathlib import Path

import pytest

from bioscript import GenotypeGenerator, GRCh, load_variants_tsv
from bioscript.data import create_test_variants


def test_assembly_enum_parsing():
    """Test GRCh.parse() with various formats."""
    # Test case-insensitive strings
    assert GRCh.parse("GRCh38") == GRCh.GRCH38
    assert GRCh.parse("grch38") == GRCh.GRCH38
    assert GRCh.parse("GRCH38") == GRCh.GRCH38

    # Test short versions
    assert GRCh.parse("38") == GRCh.GRCH38
    assert GRCh.parse("37") == GRCh.GRCH37
    assert GRCh.parse("36") == GRCh.GRCH36

    # Test with spaces
    assert GRCh.parse(" GRCh38 ") == GRCh.GRCH38

    # Test GRCh enum passthrough
    assert GRCh.parse(GRCh.GRCH38) == GRCh.GRCH38

    # Test None
    assert GRCh.parse(None) is None

    # Test invalid
    with pytest.raises(ValueError):
        GRCh.parse("GRCh99")

    with pytest.raises(ValueError):
        GRCh.parse("invalid")


def test_assembly_in_genotype_generator():
    """Test GenotypeGenerator with assembly parameter."""
    gen = GenotypeGenerator(
        [
            {"rsid": "rs73885319", "chromosome": "22", "position": 36265860},
            {"rsid": "rs60910145", "chromosome": "22", "position": 36265988},
        ],
        assembly="GRCh38",
    )

    variants = list(gen(["AA", "TT"]))
    assert len(variants) == 2
    assert variants[0].assembly == GRCh.GRCH38
    assert variants[1].assembly == GRCh.GRCH38


def test_assembly_in_create_test_variants():
    """Test create_test_variants with assembly parameter."""
    variants = [
        {"rsid": "rs123", "chromosome": "1", "position": 1000, "genotype": "AA"},
        {"rsid": "rs456", "chromosome": "2", "position": 2000, "genotype": "TT"},
    ]

    variant_rows = list(create_test_variants(variants, assembly="GRCh37"))
    assert len(variant_rows) == 2
    assert variant_rows[0].assembly == GRCh.GRCH37
    assert variant_rows[1].assembly == GRCh.GRCH37


def test_assembly_from_file_always_none():
    """Test that assembly is always None when loading from TSV files."""
    with tempfile.NamedTemporaryFile(mode="w", suffix=".txt", delete=False) as f:
        f.write("# rsid\tchromosome\tposition\tgenotype\n")
        f.write("rs123\t1\t1000\tAA\n")
        f.write("rs456\t2\t2000\tTT\n")
        temp_path = f.name

    try:
        variants = list(load_variants_tsv(temp_path))
        assert len(variants) == 2
        # assembly is always None from files - user must set when constructing
        assert variants[0].assembly is None
        assert variants[1].assembly is None
    finally:
        Path(temp_path).unlink()


def test_assembly_programmatic_creation():
    """Test that assembly is set correctly when creating variants programmatically."""
    # Test with create_test_variants
    variants = [
        {"rsid": "rs123", "chromosome": "1", "position": 1000, "genotype": "AA"},
        {"rsid": "rs456", "chromosome": "2", "position": 2000, "genotype": "TT"},
    ]

    # assembly parameter is used when creating variants in-memory
    variant_rows = list(create_test_variants(variants, assembly="GRCh38"))
    assert len(variant_rows) == 2
    assert variant_rows[0].assembly == GRCh.GRCH38
    assert variant_rows[1].assembly == GRCh.GRCH38

    # Files don't store assembly - always None when loaded back
    with tempfile.NamedTemporaryFile(mode="w", suffix=".txt", delete=False) as f:
        temp_path = f.name

    try:
        output_path = create_test_variants(variants, output_file=temp_path, assembly="GRCh38")
        loaded_variants = list(load_variants_tsv(str(output_path)))
        assert len(loaded_variants) == 2
        # assembly is NOT persisted to file - always None when loaded
        assert loaded_variants[0].assembly is None
        assert loaded_variants[1].assembly is None
    finally:
        Path(temp_path).unlink()
