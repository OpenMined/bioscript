"""APOL1 genotype classifier for kidney disease risk assessment.

This classifier identifies G0, G1, and G2 APOL1 risk variants based on:
- rs73885319: A>G (G1 variant, position 1)
- rs60910145: T>C (G1 variant, position 2)
- rs71785313: Insertion/Deletion (G2 variant)

Export convention for bioscript CLI:
- variant_calls: List of VariantCall objects to match
- classifier: Classifier instance to run
- name: Column name for output (optional, defaults to filename)
"""

from bioscript import AlleleCounter
from bioscript.classifier import DiploidResult, GenotypeClassifier, GenotypeEnum
from bioscript.types import Alleles, VariantCall

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


class APOL1Classifier(GenotypeClassifier):
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


# ==============================================================================
# CLI Export Convention
# ==============================================================================
# Export __bioscript__ dictionary for the bioscript CLI

__bioscript__ = {
    "variant_calls": [rs73885319, rs60910145, rs71785313],
    "classifier": APOL1Classifier(),
    "name": "APOL1",  # Optional, defaults to script filename
}
