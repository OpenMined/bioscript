"""Tests for complex variant matching (insertions, deletions, indels, MNVs)."""

from bioscript.types import (
    DiploidSNP,
    InDel,
    MatchType,
    Nucleotide,
    VariantCall,
    VariantType,
)


def test_variant_type_auto_detection_snv():
    """Test auto-detection of SNV (single nucleotide variant)."""
    vc = VariantCall(rsid="rs12345", ref="A", alt="G", chromosome="17", position=43094464)
    assert vc.variant_type == VariantType.SNV


def test_variant_type_auto_detection_insertion():
    """Test auto-detection of insertion (ref=1, alt>1)."""
    vc = VariantCall(rsid="rs12345", ref="A", alt="AG", chromosome="17", position=43094464)
    assert vc.variant_type == VariantType.INSERTION


def test_variant_type_auto_detection_deletion():
    """Test auto-detection of deletion (ref>1, alt=1)."""
    vc = VariantCall(rsid="rs12345", ref="AG", alt="A", chromosome="17", position=43094464)
    assert vc.variant_type == VariantType.DELETION


def test_variant_type_auto_detection_indel():
    """Test auto-detection of indel (same length, >1)."""
    vc = VariantCall(rsid="rs397509270", ref="AC", alt="CT", chromosome="17", position=43094464)
    assert vc.variant_type == VariantType.INDEL


def test_variant_type_auto_detection_mnv():
    """Test auto-detection of MNV (different lengths, both >1)."""
    vc = VariantCall(rsid="rs12345", ref="CTG", alt="TTTA", chromosome="17", position=43094464)
    assert vc.variant_type == VariantType.MNV


def test_insertion_matching_i_alleles():
    """Test insertion matching: I alleles are pathogenic."""
    vc = VariantCall(rsid="rs12345", ref="A", alt="AG", chromosome="17", position=43094464)

    # II genotype = homozygous insertion = VARIANT_CALL
    snp_ii = DiploidSNP(InDel.I, InDel.I)
    assert vc._match_insertion(snp_ii) == MatchType.VARIANT_CALL

    # DD genotype = no insertion = REFERENCE_CALL
    snp_dd = DiploidSNP(InDel.D, InDel.D)
    assert vc._match_insertion(snp_dd) == MatchType.REFERENCE_CALL

    # ID genotype = heterozygous insertion = VARIANT_CALL
    snp_id = DiploidSNP(InDel.I, InDel.D)
    assert vc._match_insertion(snp_id) == MatchType.VARIANT_CALL


def test_deletion_matching_d_alleles():
    """Test deletion matching: D alleles are pathogenic."""
    vc = VariantCall(rsid="rs12345", ref="AG", alt="A", chromosome="17", position=43094464)

    # DD genotype = homozygous deletion = VARIANT_CALL
    snp_dd = DiploidSNP(InDel.D, InDel.D)
    assert vc._match_deletion(snp_dd) == MatchType.VARIANT_CALL

    # II genotype = no deletion = REFERENCE_CALL
    snp_ii = DiploidSNP(InDel.I, InDel.I)
    assert vc._match_deletion(snp_ii) == MatchType.REFERENCE_CALL

    # ID genotype = heterozygous deletion = VARIANT_CALL
    snp_id = DiploidSNP(InDel.I, InDel.D)
    assert vc._match_deletion(snp_id) == MatchType.VARIANT_CALL


def test_complex_indel_matching_rs397509270():
    """Test complex indel matching: rs397509270 ref=AC alt=CT."""
    vc = VariantCall(rsid="rs397509270", ref="AC", alt="CT", chromosome="17", position=43094464)

    # AA genotype: A is in ref (AC), not exclusively in alt (CT) = REFERENCE_CALL
    snp_aa = DiploidSNP(Nucleotide.A, Nucleotide.A)
    assert vc._match_complex(snp_aa) == MatchType.REFERENCE_CALL

    # TT genotype: T is in alt (CT), not in ref (AC) = VARIANT_CALL
    snp_tt = DiploidSNP(Nucleotide.T, Nucleotide.T)
    assert vc._match_complex(snp_tt) == MatchType.VARIANT_CALL

    # CC genotype: C is in both ref and alt = NO_CALL (ambiguous)
    snp_cc = DiploidSNP(Nucleotide.C, Nucleotide.C)
    assert vc._match_complex(snp_cc) == MatchType.NO_CALL


def test_complex_mnv_matching_rs397508883():
    """Test complex indel matching: rs397508883 ref=CTGC alt=TTTA (same length=4)."""
    vc = VariantCall(rsid="rs397508883", ref="CTGC", alt="TTTA", chromosome="17", position=43094464)

    # AA genotype: A is in alt (TTTA), not in ref (CTGC) = VARIANT_CALL
    snp_aa = DiploidSNP(Nucleotide.A, Nucleotide.A)
    assert vc._match_complex(snp_aa) == MatchType.VARIANT_CALL

    # GG genotype: G is in ref (CTGC), not in alt (TTTA) = REFERENCE_CALL
    snp_gg = DiploidSNP(Nucleotide.G, Nucleotide.G)
    assert vc._match_complex(snp_gg) == MatchType.REFERENCE_CALL

    # TT genotype: T is in both ref and alt = NO_CALL (ambiguous)
    snp_tt = DiploidSNP(Nucleotide.T, Nucleotide.T)
    assert vc._match_complex(snp_tt) == MatchType.NO_CALL


def test_complex_matching_with_missing_alleles():
    """Test complex matching handles missing alleles correctly."""
    vc = VariantCall(rsid="rs12345", ref="AC", alt="CT", chromosome="17", position=43094464)

    # Missing genotype = NO_CALL
    snp_missing = DiploidSNP(Nucleotide.MISSING, Nucleotide.MISSING)
    assert vc._match_complex(snp_missing) == MatchType.NO_CALL


def test_complex_matching_heterozygous():
    """Test complex matching with heterozygous genotypes."""
    vc = VariantCall(rsid="rs12345", ref="AC", alt="GT", chromosome="17", position=43094464)

    # AG genotype: A in ref, G in alt = NO_CALL (mixed)
    snp_ag = DiploidSNP(Nucleotide.A, Nucleotide.G)
    assert vc._match_complex(snp_ag) == MatchType.NO_CALL

    # GG genotype: G only in alt = VARIANT_CALL
    snp_gg = DiploidSNP(Nucleotide.G, Nucleotide.G)
    assert vc._match_complex(snp_gg) == MatchType.VARIANT_CALL


def test_determine_match_type_uses_variant_type():
    """Test that _determine_match_type dispatches to correct method based on variant_type."""
    # SNV
    vc_snv = VariantCall(rsid="rs12345", ref="A", alt="G", chromosome="17", position=43094464)
    snp_gg = DiploidSNP(Nucleotide.G, Nucleotide.G)
    assert vc_snv._determine_match_type(snp_gg) == MatchType.VARIANT_CALL

    # Insertion
    vc_ins = VariantCall(rsid="rs12345", ref="A", alt="AG", chromosome="17", position=43094464)
    snp_ii = DiploidSNP(InDel.I, InDel.I)
    assert vc_ins._determine_match_type(snp_ii) == MatchType.VARIANT_CALL

    # Deletion
    vc_del = VariantCall(rsid="rs12345", ref="AG", alt="A", chromosome="17", position=43094464)
    snp_dd = DiploidSNP(InDel.D, InDel.D)
    assert vc_del._determine_match_type(snp_dd) == MatchType.VARIANT_CALL

    # Complex indel
    vc_indel = VariantCall(
        rsid="rs397509270", ref="AC", alt="CT", chromosome="17", position=43094464
    )
    snp_tt = DiploidSNP(Nucleotide.T, Nucleotide.T)
    assert vc_indel._determine_match_type(snp_tt) == MatchType.VARIANT_CALL


def test_ref_and_alt_string_storage():
    """Test that ref/alt strings are stored directly and variant type is detected."""
    vc = VariantCall(rsid="rs397508883", ref="CTGC", alt="TTTA", chromosome="17", position=43094464)

    # Strings should be stored as-is
    assert vc.ref == "CTGC"
    assert vc.alt == "TTTA"

    # Variant type should be auto-detected (same length = INDEL)
    assert vc.variant_type == VariantType.INDEL


def test_snv_heterozygous_false_positive_bug():
    """
    Test that heterozygous genotypes don't falsely match as VARIANT_CALL
    when neither allele matches the alt allele.
    """
    # Variant: ref=A, alt=G
    vc = VariantCall(rsid="rs12345", ref="A", alt="G", chromosome="17", position=43094464)

    # CT genotype: neither C nor T is the alt allele G
    # Should be NO_CALL, not VARIANT_CALL
    snp_ct = DiploidSNP(Nucleotide.C, Nucleotide.T)
    assert vc._match_snv(snp_ct) == MatchType.NO_CALL

    # AG genotype: one allele (G) matches alt
    # Should be VARIANT_CALL
    snp_ag = DiploidSNP(Nucleotide.A, Nucleotide.G)
    assert vc._match_snv(snp_ag) == MatchType.VARIANT_CALL

    # GG genotype: both alleles match alt
    # Should be VARIANT_CALL
    snp_gg = DiploidSNP(Nucleotide.G, Nucleotide.G)
    assert vc._match_snv(snp_gg) == MatchType.VARIANT_CALL

    # AA genotype: homozygous ref
    # Should be REFERENCE_CALL
    snp_aa = DiploidSNP(Nucleotide.A, Nucleotide.A)
    assert vc._match_snv(snp_aa) == MatchType.REFERENCE_CALL
