"""Tests for APOL1 genotype classification.

APOL1 (Apolipoprotein L1) genetic variants are associated with kidney disease risk.
This test suite validates the classification logic for G0, G1, and G2 genotypes.

Reference SNPs:
- rs73885319: A>G (G1 variant, position 1)
- rs60910145: T>C (G1 variant, position 2)
- rs71785313: Insertion/Deletion (G2 variant)

Classification Rules (count-based, no phase information):
- G2 requires deletion variant at rs71785313 (D allele)
- G1 requires variants at BOTH rs73885319 AND rs60910145
- Counts are simply allele counts (0-4 total for G1)
"""

from bioscript import AlleleCounter, GenotypeGenerator
from bioscript.classifier import Classifier, DiploidResult, GenotypeEnum
from bioscript.types import Alleles, MatchList, VariantCall

# Define APOL1 variant calls
# rs73885319: A>G at chr22:36265860 (GRCh38)
rs73885319 = VariantCall(rsid="rs73885319", ref=Alleles.A, alt=Alleles.NOT_A)

# rs60910145: T>C at chr22:36265988 (GRCh38)
rs60910145 = VariantCall(rsid="rs60910145", ref=Alleles.T, alt=Alleles.NOT_T)

# rs71785313: INDEL at chr22:36266000 (GRCh38)
# Has multiple rsID aliases
rs71785313 = VariantCall(
    rsid=["rs71785313", "rs1317778148", "rs143830837"], ref=Alleles.I, alt=Alleles.D
)


# Define APOL1 genotype categories
class APOL1Genotypes(GenotypeEnum):
    G2 = "G2"
    G1 = "G1"
    G0 = "G0"


MISSING = "G-"


class APOL1Classifier(Classifier):
    """
    Classify APOL1 genotypes based on simple allele counting.

    Without phase information, we use a count-based approach:
    - Count D alleles at rs71785313 (0, 1, or 2)
    - Count variant alleles at BOTH G1 positions (0-4 total)
    - G1 only counts if variants present at BOTH sites
    """

    def classify(self, matches) -> DiploidResult:
        # Create counters for each position
        g2_counter = AlleleCounter(rs71785313)
        g1_site1_counter = AlleleCounter(rs73885319)
        g1_site2_counter = AlleleCounter(rs60910145)

        # Count variant alleles
        g2_result = g2_counter.count(matches)
        site1_result = g1_site1_counter.count(matches)
        site2_result = g1_site2_counter.count(matches)

        # Check if we have any APOL1 data
        has_data = (
            g2_result.genotype is not None
            or site1_result.genotype is not None
            or site2_result.genotype is not None
        )
        if not has_data:
            return DiploidResult(MISSING, MISSING)

        d_count = g2_result.alt_count  # D alleles (0, 1, or 2)

        # G1 requires variants at BOTH positions
        site1_variants = site1_result.alt_count  # 0, 1, or 2
        site2_variants = site2_result.alt_count  # 0, 1, or 2

        # Only count as G1 if both sites have at least one variant
        has_g1 = site1_variants > 0 and site2_variants > 0
        g1_total = site1_variants + site2_variants if has_g1 else 0  # 0, 2, 3, or 4

        # Simple count-based classification
        if d_count == 2:  # Homozygous deletion
            return DiploidResult(APOL1Genotypes.G2, APOL1Genotypes.G2)
        elif d_count == 1:  # Heterozygous deletion
            if g1_total >= 2:  # At least one G1 copy
                return DiploidResult(APOL1Genotypes.G2, APOL1Genotypes.G1)
            else:
                return DiploidResult(APOL1Genotypes.G2, APOL1Genotypes.G0)
        else:  # No deletion
            if g1_total == 4:  # Both sites homozygous variant
                return DiploidResult(APOL1Genotypes.G1, APOL1Genotypes.G1)
            elif g1_total >= 2:  # At least one G1 copy
                return DiploidResult(APOL1Genotypes.G1, APOL1Genotypes.G0)
            else:
                return DiploidResult(APOL1Genotypes.G0, APOL1Genotypes.G0)


# Test fixture: genotype generator for APOL1 SNPs
apol1_genotypes = GenotypeGenerator(
    [
        {"rsid": "rs73885319", "chromosome": "22", "position": 36265860},
        {"rsid": "rs60910145", "chromosome": "22", "position": 36265988},
        {"rsid": "rs71785313", "chromosome": "22", "position": 36266000},
    ]
)


def test_apol1_g0_g0_healthy():
    """Test G0/G0 (healthy, no risk variants)."""
    classifier = APOL1Classifier()

    # Wild-type genotypes: AA, TT, II (no variants)
    variants = apol1_genotypes(["AA", "TT", "II"])

    calls = MatchList(variant_calls=[rs73885319, rs60910145, rs71785313])
    matches = calls.match_rows(variants)
    result = classifier.classify(matches)

    assert str(result) == "G0/G0"
    assert result.genotype1 == APOL1Genotypes.G0
    assert result.genotype2 == APOL1Genotypes.G0


def test_apol1_g1_g0_heterozygous():
    """Test G1/G0 (heterozygous for G1 variants)."""
    classifier = APOL1Classifier()

    # Heterozygous for both G1 SNPs, no deletion
    variants = apol1_genotypes(["AG", "TC", "II"])

    calls = MatchList(variant_calls=[rs73885319, rs60910145, rs71785313])
    matches = calls.match_rows(variants)
    result = classifier.classify(matches)

    assert str(result) == "G1/G0"
    assert result.genotype1 == APOL1Genotypes.G1
    assert result.genotype2 == APOL1Genotypes.G0


def test_apol1_g1_g1_homozygous():
    """Test G1/G1 (homozygous for G1 variants)."""
    classifier = APOL1Classifier()

    # Homozygous for both G1 SNPs
    variants = apol1_genotypes(["GG", "CC", "II"])

    calls = MatchList(variant_calls=[rs73885319, rs60910145, rs71785313])
    matches = calls.match_rows(variants)
    result = classifier.classify(matches)

    assert str(result) == "G1/G1"
    assert result.genotype1 == APOL1Genotypes.G1
    assert result.genotype2 == APOL1Genotypes.G1


def test_apol1_g2_g0_heterozygous_deletion():
    """Test G2/G0 (heterozygous deletion, no G1)."""
    classifier = APOL1Classifier()

    # Wild-type for G1 SNPs, heterozygous deletion
    variants = apol1_genotypes(["AA", "TT", "ID"])

    calls = MatchList(variant_calls=[rs73885319, rs60910145, rs71785313])
    matches = calls.match_rows(variants)
    result = classifier.classify(matches)

    assert str(result) == "G2/G0"
    assert result.genotype1 == APOL1Genotypes.G2
    assert result.genotype2 == APOL1Genotypes.G0


def test_apol1_g2_g2_homozygous_deletion():
    """Test G2/G2 (homozygous deletion)."""
    classifier = APOL1Classifier()

    # Homozygous deletion
    variants = apol1_genotypes(["AA", "TT", "DD"])

    calls = MatchList(variant_calls=[rs73885319, rs60910145, rs71785313])
    matches = calls.match_rows(variants)
    result = classifier.classify(matches)

    assert str(result) == "G2/G2"
    assert result.genotype1 == APOL1Genotypes.G2
    assert result.genotype2 == APOL1Genotypes.G2


def test_apol1_g2_g1_compound_heterozygous():
    """Test G2/G1 (compound heterozygous - highest risk)."""
    classifier = APOL1Classifier()

    # Heterozygous deletion + homozygous G1 variants
    variants = apol1_genotypes(["GG", "CC", "ID"])

    calls = MatchList(variant_calls=[rs73885319, rs60910145, rs71785313])
    matches = calls.match_rows(variants)
    result = classifier.classify(matches)

    assert str(result) == "G2/G1"
    assert result.genotype1 == APOL1Genotypes.G2
    assert result.genotype2 == APOL1Genotypes.G1


def test_apol1_g2_g1_het_deletion_het_g1():
    """Test heterozygous deletion with heterozygous G1 variants.

    Without phase information, we just count alleles:
    - 1 D allele at rs71785313
    - 1 variant at rs73885319 + 1 variant at rs60910145 = 2 total
    - Since both G1 sites have variants (total >= 2), this is G2/G1
    """
    classifier = APOL1Classifier()

    # Het G1 variants + het deletion
    variants = apol1_genotypes(["AG", "TC", "ID"])

    calls = MatchList(variant_calls=[rs73885319, rs60910145, rs71785313])
    matches = calls.match_rows(variants)
    result = classifier.classify(matches)

    # Simple count: 1 D allele + 2 G1 alleles (1 at each site) = G2/G1
    assert str(result) == "G2/G1"
    assert result.genotype1 == APOL1Genotypes.G2
    assert result.genotype2 == APOL1Genotypes.G1


def test_apol1_missing_no_data():
    """Test that missing APOL1 data returns G-/G-."""
    classifier = APOL1Classifier()

    # Create variants for different SNPs (not APOL1)
    other_gen = GenotypeGenerator(
        [
            {"rsid": "rs999999", "chromosome": "1", "position": 1000},
        ]
    )
    variants = other_gen(["AA"])

    calls = MatchList(variant_calls=[rs73885319, rs60910145, rs71785313])
    matches = calls.match_rows(variants)
    result = classifier.classify(matches)

    assert str(result) == "G-/G-"


def test_apol1_partial_g1_not_classified():
    """Test that partial G1 (only one SNP variant) doesn't count as G1.

    G1 requires BOTH rs73885319 AND rs60910145 variants in CIS.
    Having only one should result in G0/G0.
    """
    classifier = APOL1Classifier()

    # Only first G1 SNP is variant, second is wild-type
    variants = apol1_genotypes(["AG", "TT", "II"])

    calls = MatchList(variant_calls=[rs73885319, rs60910145, rs71785313])
    matches = calls.match_rows(variants)
    result = classifier.classify(matches)

    # Should be G0/G0 because G1 requires both SNPs
    assert str(result) == "G0/G0"
    assert result.genotype1 == APOL1Genotypes.G0
    assert result.genotype2 == APOL1Genotypes.G0
