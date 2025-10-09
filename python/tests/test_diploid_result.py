"""Tests for DiploidResult sorting functionality."""

from bioscript.classifier import DiploidResult, GenotypeEnum


class SampleGenotypes(GenotypeEnum):
    HIGH = "HIGH"
    MEDIUM = "MEDIUM"
    LOW = "LOW"


def test_diploid_result_sorted_respects_enum_order():
    """Test that sorted() returns genotypes in enum definition order."""
    # LOW comes after HIGH in enum, so LOW/HIGH should become HIGH/LOW
    result = DiploidResult(SampleGenotypes.LOW, SampleGenotypes.HIGH)
    sorted_result = result.sorted()

    assert sorted_result.genotype1 == SampleGenotypes.HIGH
    assert sorted_result.genotype2 == SampleGenotypes.LOW
    assert str(sorted_result) == "HIGH/LOW"


def test_diploid_result_sorted_already_sorted():
    """Test that already sorted results remain unchanged."""
    result = DiploidResult(SampleGenotypes.HIGH, SampleGenotypes.MEDIUM)
    sorted_result = result.sorted()

    assert sorted_result.genotype1 == SampleGenotypes.HIGH
    assert sorted_result.genotype2 == SampleGenotypes.MEDIUM
    assert str(sorted_result) == "HIGH/MEDIUM"


def test_diploid_result_sorted_same_genotype():
    """Test sorting with identical genotypes."""
    result = DiploidResult(SampleGenotypes.MEDIUM, SampleGenotypes.MEDIUM)
    sorted_result = result.sorted()

    assert sorted_result.genotype1 == SampleGenotypes.MEDIUM
    assert sorted_result.genotype2 == SampleGenotypes.MEDIUM
    assert str(sorted_result) == "MEDIUM/MEDIUM"


def test_diploid_result_sorted_reverse_order():
    """Test sorting with genotypes in reverse order."""
    result = DiploidResult(SampleGenotypes.LOW, SampleGenotypes.HIGH)
    sorted_result = result.sorted()

    assert sorted_result.genotype1 == SampleGenotypes.HIGH
    assert sorted_result.genotype2 == SampleGenotypes.LOW


def test_diploid_result_sorted_non_enum():
    """Test that non-enum values are returned as-is."""
    result = DiploidResult("A", "B")
    sorted_result = result.sorted()

    # Non-enums should be returned unchanged
    assert sorted_result.genotype1 == "A"
    assert sorted_result.genotype2 == "B"
