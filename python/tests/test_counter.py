"""Tests for allele counting utilities."""

from bioscript import AlleleCounter, GenotypeGenerator
from bioscript.types import Alleles, MatchList, Nucleotide, VariantCall


def test_allele_counter_heterozygous():
    """Test counting heterozygous genotype."""
    # A>G variant
    call = VariantCall(rsid="rs123", ref=Alleles.A, alt=Alleles.G)
    counter = AlleleCounter(call)

    # Create AG genotype
    gen = GenotypeGenerator([{"rsid": "rs123", "chromosome": "1", "position": 1000}])
    variants = gen(["AG"])

    matches = MatchList(variant_calls=[call])
    matches.match_rows(variants)

    result = counter.count(matches)

    assert result.ref_count == 1
    assert result.alt_count == 1
    assert result.has_variant()
    assert result.is_heterozygous()
    assert not result.is_homozygous_variant()
    assert not result.is_homozygous_reference()
    assert result.count(Nucleotide.A) == 1
    assert result.count(Nucleotide.G) == 1


def test_allele_counter_homozygous_variant():
    """Test counting homozygous variant genotype."""
    call = VariantCall(rsid="rs123", ref=Alleles.A, alt=Alleles.G)
    counter = AlleleCounter(call)

    gen = GenotypeGenerator([{"rsid": "rs123", "chromosome": "1", "position": 1000}])
    variants = gen(["GG"])

    matches = MatchList(variant_calls=[call])
    matches.match_rows(variants)

    result = counter.count(matches)

    assert result.ref_count == 0
    assert result.alt_count == 2
    assert result.has_variant()
    assert result.is_homozygous_variant()
    assert not result.is_heterozygous()
    assert result.count(Nucleotide.G) == 2


def test_allele_counter_homozygous_reference():
    """Test counting homozygous reference genotype."""
    call = VariantCall(rsid="rs123", ref=Alleles.A, alt=Alleles.G)
    counter = AlleleCounter(call)

    gen = GenotypeGenerator([{"rsid": "rs123", "chromosome": "1", "position": 1000}])
    variants = gen(["AA"])

    matches = MatchList(variant_calls=[call])
    matches.match_rows(variants)

    result = counter.count(matches)

    assert result.ref_count == 2
    assert result.alt_count == 0
    assert not result.has_variant()
    assert result.is_homozygous_reference()
    assert not result.is_heterozygous()
    assert result.count(Nucleotide.A) == 2


def test_allele_counter_no_match():
    """Test counter when variant not found."""
    call = VariantCall(rsid="rs123", ref=Alleles.A, alt=Alleles.G)
    counter = AlleleCounter(call)

    # Different rsID
    gen = GenotypeGenerator([{"rsid": "rs999", "chromosome": "1", "position": 1000}])
    variants = gen(["AG"])

    matches = MatchList(variant_calls=[call])
    matches.match_rows(variants)

    result = counter.count(matches)

    assert result.ref_count == 0
    assert result.alt_count == 0
    assert not result.has_variant()
    assert result.genotype is None


def test_allele_counter_multiple_alt_alleles():
    """Test with multiple possible alternate alleles."""
    # A can be anything but A
    call = VariantCall(rsid="rs123", ref=Alleles.A, alt=Alleles.NOT_A)
    counter = AlleleCounter(call)

    gen = GenotypeGenerator([{"rsid": "rs123", "chromosome": "1", "position": 1000}])
    variants = gen(["TC"])  # Both are NOT_A

    matches = MatchList(variant_calls=[call])
    matches.match_rows(variants)

    result = counter.count(matches)

    assert result.ref_count == 0
    assert result.alt_count == 2
    assert result.is_homozygous_variant()
    assert result.count(Nucleotide.T) == 1
    assert result.count(Nucleotide.C) == 1


def test_allele_counter_repr():
    """Test string representations."""
    call = VariantCall(rsid="rs123", ref=Alleles.A, alt=Alleles.G)
    counter = AlleleCounter(call)

    assert "rs123" in repr(counter)

    gen = GenotypeGenerator([{"rsid": "rs123", "chromosome": "1", "position": 1000}])
    variants = gen(["AG"])

    matches = MatchList(variant_calls=[call])
    matches.match_rows(variants)

    result = counter.count(matches)

    assert "AG" in repr(result)
    assert "ref=1" in repr(result)
    assert "alt=1" in repr(result)
