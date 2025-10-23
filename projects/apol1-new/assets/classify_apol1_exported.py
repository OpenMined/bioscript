from bioscript import AlleleCounter
from bioscript.classifier import DiploidResult, GenotypeClassifier, GenotypeEnum
from bioscript.types import Alleles, VariantCall

# Define APOL1 variant calls
rs73885319 = VariantCall(rsid="rs73885319", ref=Alleles.A, alt=Alleles.NOT_A)
rs60910145 = VariantCall(rsid="rs60910145", ref=Alleles.T, alt=Alleles.NOT_T)
rs71785313 = VariantCall(
    rsid=["rs71785313", "rs1317778148", "rs143830837"], ref=Alleles.I, alt=Alleles.D
)

class APOL1Genotypes(GenotypeEnum):
    G2 = "G2"
    G1 = "G1"
    G0 = "G0"

MISSING = "G-"

class APOL1Classifier(GenotypeClassifier):
    def classify(self, matches) -> DiploidResult:
        g2_counter = AlleleCounter(rs71785313)
        g1_site1_counter = AlleleCounter(rs73885319)
        g1_site2_counter = AlleleCounter(rs60910145)

        g2_result = g2_counter.count(matches)
        site1_result = g1_site1_counter.count(matches)
        site2_result = g1_site2_counter.count(matches)

        has_data = (
            g2_result.genotype is not None
            or site1_result.genotype is not None
            or site2_result.genotype is not None
        )
        if not has_data:
            return DiploidResult(MISSING, MISSING)

        d_count = g2_result.alt_count
        site1_variants = site1_result.alt_count
        site2_variants = site2_result.alt_count

        has_g1 = site1_variants > 0 and site2_variants > 0
        g1_total = site1_variants + site2_variants if has_g1 else 0

        if d_count == 2:
            return DiploidResult(APOL1Genotypes.G2, APOL1Genotypes.G2)
        elif d_count == 1:
            if g1_total >= 2:
                return DiploidResult(APOL1Genotypes.G2, APOL1Genotypes.G1)
            else:
                return DiploidResult(APOL1Genotypes.G2, APOL1Genotypes.G0)
        else:
            if g1_total == 4:
                return DiploidResult(APOL1Genotypes.G1, APOL1Genotypes.G1)
            elif g1_total >= 2:
                return DiploidResult(APOL1Genotypes.G1, APOL1Genotypes.G0)
            else:
                return DiploidResult(APOL1Genotypes.G0, APOL1Genotypes.G0)

__bioscript__ = {
    "variant_calls": [rs73885319, rs60910145, rs71785313],
    "classifier": APOL1Classifier(),
    "name": "APOL1",
}

from bioscript import VariantFixture
from bioscript.types import MatchList

fixture = VariantFixture(
    [
        {"rsid": "rs73885319", "chromosome": "22", "position": 36265860},
        {"rsid": "rs60910145", "chromosome": "22", "position": 36265988},
        {"rsid": "rs71785313", "chromosome": "22", "position": 36266000},
    ],
    assembly="GRCh38",
)

def test_g0_homozygous():
    variants = fixture(["AA", "TT", "II"])
    matches = MatchList([rs73885319, rs60910145, rs71785313]).match_rows(variants)
    classifier = APOL1Classifier()
    result = classifier(matches)
    assert result == "G0/G0"

def test_g1_homozygous():
    variants = fixture(["GG", "CC", "II"])
    matches = MatchList([rs73885319, rs60910145, rs71785313]).match_rows(variants)
    classifier = APOL1Classifier()
    result = classifier(matches)
    assert result == "G1/G1"
