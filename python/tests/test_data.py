"""Tests for data module."""

import tempfile
from pathlib import Path

from bioscript import GenotypeGenerator, create_test_variants, load_variants_tsv


def test_create_test_variants_in_memory():
    """Test creating variants in-memory."""
    variants = list(
        create_test_variants(
            [
                {"rsid": "rs123", "chromosome": "1", "position": 1000, "genotype": "AA"},
                {"rsid": "rs456", "chromosome": "2", "position": 2000, "genotype": "AT"},
            ]
        )
    )

    assert len(variants) == 2
    assert variants[0].rsid == "rs123"
    assert variants[0].chromosome == "1"
    assert variants[0].position == 1000
    assert variants[0].genotype == "AA"
    assert variants[1].rsid == "rs456"


def test_create_test_variants_with_optional_fields():
    """Test creating variants with optional fields."""
    variants = list(
        create_test_variants(
            [
                {
                    "rsid": "rs123",
                    "chromosome": "1",
                    "position": 1000,
                    "genotype": "AA",
                    "gs": 0.95,
                    "baf": 0.5,
                    "lrr": 0.1,
                },
            ]
        )
    )

    assert variants[0].gs == 0.95
    assert variants[0].baf == 0.5
    assert variants[0].lrr == 0.1


def test_create_test_variants_to_file():
    """Test creating variants and writing to file."""
    with tempfile.NamedTemporaryFile(mode="w", suffix=".txt", delete=False) as f:
        temp_path = f.name

    try:
        result_path = create_test_variants(
            [
                {"rsid": "rs123", "chromosome": "1", "position": 1000, "genotype": "AA"},
                {"rsid": "rs456", "chromosome": "2", "position": 2000, "genotype": "AT"},
            ],
            output_file=temp_path,
        )

        assert isinstance(result_path, Path)
        assert result_path.exists()

        # Read back using load_variants_tsv
        variants = list(load_variants_tsv(str(result_path)))
        assert len(variants) == 2
        assert variants[0].rsid == "rs123"
        assert variants[1].rsid == "rs456"
    finally:
        Path(temp_path).unlink(missing_ok=True)


def test_create_test_variants_missing_fields():
    """Test that missing required fields raise an error."""
    import pytest

    with pytest.raises(ValueError, match="missing required fields"):
        list(
            create_test_variants(
                [
                    {"rsid": "rs123", "chromosome": "1"},  # Missing position and genotype
                ]
            )
        )


def test_genotype_generator():
    """Test GenotypeGenerator basic functionality."""
    gen = GenotypeGenerator(
        [
            {"rsid": "rs123", "chromosome": "1", "position": 1000},
            {"rsid": "rs456", "chromosome": "2", "position": 2000},
        ]
    )

    # Generate first set of variants
    variants1 = list(gen(["AA", "TT"]))
    assert len(variants1) == 2
    assert variants1[0].rsid == "rs123"
    assert variants1[0].genotype == "AA"
    assert variants1[1].rsid == "rs456"
    assert variants1[1].genotype == "TT"

    # Generate second set with different genotypes
    variants2 = list(gen(["AT", "TC"]))
    assert len(variants2) == 2
    assert variants2[0].genotype == "AT"
    assert variants2[1].genotype == "TC"

    # Same positions
    assert variants2[0].rsid == "rs123"
    assert variants2[1].rsid == "rs456"


def test_genotype_generator_length_mismatch():
    """Test that mismatched genotype count raises error."""
    import pytest

    gen = GenotypeGenerator(
        [
            {"rsid": "rs123", "chromosome": "1", "position": 1000},
        ]
    )

    with pytest.raises(ValueError, match="must match"):
        list(gen(["AA", "TT"]))  # Too many genotypes


def test_genotype_generator_with_file():
    """Test GenotypeGenerator writing to file."""
    with tempfile.NamedTemporaryFile(mode="w", suffix=".txt", delete=False) as f:
        temp_path = f.name

    try:
        gen = GenotypeGenerator(
            [
                {"rsid": "rs123", "chromosome": "1", "position": 1000},
                {"rsid": "rs456", "chromosome": "2", "position": 2000},
            ]
        )

        result_path = gen(["AA", "TT"], output_file=temp_path)
        assert isinstance(result_path, Path)

        # Read back
        variants = list(load_variants_tsv(str(result_path)))
        assert len(variants) == 2
        assert variants[0].genotype == "AA"
    finally:
        Path(temp_path).unlink(missing_ok=True)


def test_genotype_generator_invalid_template():
    """Test that templates with genotype field raise error."""
    import pytest

    with pytest.raises(ValueError, match="should not include 'genotype'"):
        GenotypeGenerator(
            [
                {"rsid": "rs123", "chromosome": "1", "position": 1000, "genotype": "AA"},
            ]
        )
