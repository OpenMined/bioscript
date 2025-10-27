"""Tests for VariantMatch helper properties and MatchList lookup ergonomics."""

from bioscript import GenotypeGenerator, load_variants_tsv
from bioscript.types import Alleles, MatchList, Nucleotide, VariantCall


def _build_matches(call: VariantCall, genotype: str) -> MatchList:
    rsid_value = sorted(call.rsid.aliases)[0] if hasattr(call.rsid, "aliases") else call.rsid

    templates = [
        {"rsid": rsid_value, "chromosome": "1", "position": 1000},
    ]
    gen = GenotypeGenerator(templates)
    variants = gen([genotype])

    matches = MatchList(variant_calls=[call])
    matches.match_rows(variants)
    return matches


def test_variant_match_counts_heterozygous():
    call = VariantCall(rsid="rs123", ref=Alleles.A, alt=Alleles.G)

    matches = _build_matches(call, "AG")
    match = matches.get(call)

    assert match is not None
    assert match.ref_count == 1
    assert match.alt_count == 1
    assert match.has_variant
    assert match.is_heterozygous
    assert not match.is_homozygous_variant
    assert not match.is_homozygous_reference
    assert match.count(Nucleotide.A) == 1
    assert match.count(Nucleotide.G) == 1
    assert match.genotype_sorted == "AG"
    assert match.raw_line is None


def test_variant_match_counts_homozygous_variant():
    call = VariantCall(rsid="rs123", ref=Alleles.A, alt=Alleles.G)
    matches = _build_matches(call, "GG")

    match = matches.rs123

    assert match.ref_count == 0
    assert match.alt_count == 2
    assert match.has_variant
    assert match.is_homozygous_variant
    assert not match.is_heterozygous
    assert match.count(Nucleotide.G) == 2


def test_variant_match_counts_homozygous_reference():
    call = VariantCall(rsid="rs123", ref=Alleles.A, alt=Alleles.G)
    matches = _build_matches(call, "AA")

    match = matches.get("rs123")
    assert match.ref_count == 2
    assert match.alt_count == 0
    assert not match.has_variant
    assert match.is_homozygous_reference
    assert not match.is_heterozygous
    assert match.count(Nucleotide.A) == 2


def test_match_lookup_missing_returns_none():
    call = VariantCall(rsid="rs123", ref=Alleles.A, alt=Alleles.G)
    matches = _build_matches(call, "AA")

    assert matches.get("rs999") is None


def test_variant_match_multiple_alt_alleles():
    call = VariantCall(rsid="rs123", ref=Alleles.A, alt=Alleles.NOT_A)
    matches = _build_matches(call, "TC")

    match = matches.get(call)
    assert match.ref_count == 0
    assert match.alt_count == 2
    assert match.is_homozygous_variant
    assert match.count(Nucleotide.T) == 1
    assert match.count(Nucleotide.C) == 1


def test_match_lookup_uses_coordinate_aliases():
    call = VariantCall(
        rsid="rs123",
        ref=Alleles.A,
        alt=Alleles.G,
        chromosome="1",
        position=1000,
    )
    matches = _build_matches(call, "AG")

    assert matches.get("chr1:1000").genotype_sorted == "AG"
    assert matches.get("1:1000").genotype_sorted == "AG"


def test_variant_match_raw_line(tmp_path):
    call = VariantCall(rsid="rs123", ref=Alleles.A, alt=Alleles.G)

    path = tmp_path / "variants.tsv"
    path.write_text("# rsid\tchromosome\tposition\tgenotype\nrs123\t1\t1000\tAG\n")

    variants = list(load_variants_tsv(str(path)))
    matches = MatchList(variant_calls=[call]).match_rows(variants)

    match = matches.get(call)
    assert match is not None
    assert match.genotype_sorted == "AG"
    assert match.raw_line == "rs123\t1\t1000\tAG"
