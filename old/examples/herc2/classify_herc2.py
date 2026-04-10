from bioscript import write_tsv
from bioscript.classifier import GenotypeClassifier
from bioscript.types import Alleles, MatchList, MatchType, Nucleotide, VariantCall

# def rs12913832 = 38/37 ['15:28120472-28120472', '15:28365618-28365618']
# https://www.ncbi.nlm.nih.gov/snp/?term=rs12913832
# https://www.ncbi.nlm.nih.gov/snp/rs12913832
rs12913832 = VariantCall(rsid=["rs12913832", "rs60078917"], ref=Alleles.A, alt=Alleles.NOT_A, gene="HERC2")
# (A;A) yields brown eye color ~80% of the time.
# (A;G) also tends toward brown.
# (G;G) gives blue eye color ~99% of the time.

def _format_allele_label(allele):
    if isinstance(allele, Alleles):
        return ",".join(sorted(a.value for a in allele))
    if isinstance(allele, str):
        return allele
    if hasattr(allele, "__iter__"):
        return ",".join(sorted(str(a) for a in allele))
    return str(allele)

class HERC2Classifier(GenotypeClassifier):
    def classify(self, matches):
        match = matches.get(rs12913832)
        
        # Determine eye color result
        eye_color_map = {
            "AA": "Brown",
            "AG": "Brown",
            "GG": "Blue",
        }
        
        if not match or match.has_missing:
            result = "No call"
            genotype_sorted = None
            match_type = MatchType.NO_CALL.value
            ref_count = 0
            alt_count = 0
        else:
            genotype_sorted = match.genotype_sorted
            result = eye_color_map.get(genotype_sorted, "Unknown")
            match_type = match.match_type.value
            ref_count = match.ref_count
            alt_count = match.alt_count
        
        # Extract properties from match (source_row has the actual data from the file)
        if match and match.source_row:
            rsid = match.source_row.rsid
            chromosome = match.source_row.chromosome
            position = match.source_row.position
        else:
            rsid = "rs12913832"
            chromosome = None
            position = None
        
        # Create report row using match properties
        ref_label = _format_allele_label(rs12913832.ref)
        alt_label = _format_allele_label(rs12913832.alt)
        report_row = {
            "participant_id": self.participant_id,
            "filename": self.filename,
            "rsid": rsid,
            "chromosome": chromosome,
            "position": position,
            "ref": ref_label,
            "alt": alt_label,
            "genotype": genotype_sorted,
            "match_type": match_type,
            "eye_color": result,
            "ref_count": ref_count,
            "alt_count": alt_count,
        }
        
        # Write to TSV file (as a list with one row)
        write_tsv(f"{self.output_basename}.tsv", [report_row])
        
        # Return list for testing (consistent with BRCA)
        return [report_row]

__bioscript__ = {
    "variant_calls": [rs12913832],
    "classifier": HERC2Classifier,
    "name": "HERC2",
}

from bioscript import VariantFixture

# Use the regular VariantFixture which now includes raw_line functionality
fixture = VariantFixture(
    [
        {"rsid": "rs12913832", "chromosome": "15", "position": 28120472}
    ],
    assembly="GRCh38",
)

import os

def classify_fixture(genotype):
    variants = fixture([genotype])
    matches = MatchList([rs12913832]).match_rows(variants)
    classifier = HERC2Classifier(participant_id="TEST_ID", name="HERC2", filename="test.txt")
    return classifier(matches)

def test_brown_homozygous():
    result = classify_fixture("AA")
    assert len(result) == 1
    assert result[0]["eye_color"] == "Brown"
    assert result[0]["genotype"] == "AA"
    assert result[0]["match_type"] == MatchType.REFERENCE_CALL.value
    assert result[0]["participant_id"] == "TEST_ID"
    assert result[0]["filename"] == "test.txt"
    assert result[0]["rsid"] == "rs12913832"
    assert result[0]["chromosome"] == "15"
    assert result[0]["position"] == 28120472
    assert result[0]["ref_count"] == 2
    assert result[0]["alt_count"] == 0
    # Cleanup
    os.remove("result_HERC2_TEST_ID.tsv")

def test_brown_heterozygous_unsorted():
    result = classify_fixture("GA")
    assert len(result) == 1
    assert result[0]["eye_color"] == "Brown"
    assert result[0]["genotype"] == "AG"
    assert result[0]["match_type"] == MatchType.VARIANT_CALL.value
    assert result[0]["participant_id"] == "TEST_ID"
    assert result[0]["filename"] == "test.txt"
    assert result[0]["ref_count"] == 1
    assert result[0]["alt_count"] == 1
    # Cleanup
    os.remove("result_HERC2_TEST_ID.tsv")

def test_blue_homozygous():
    result = classify_fixture("GG")
    assert len(result) == 1
    assert result[0]["eye_color"] == "Blue"
    assert result[0]["genotype"] == "GG"
    assert result[0]["match_type"] == MatchType.VARIANT_CALL.value
    assert result[0]["participant_id"] == "TEST_ID"
    assert result[0]["filename"] == "test.txt"
    assert result[0]["ref_count"] == 0
    assert result[0]["alt_count"] == 2
    # Cleanup
    os.remove("result_HERC2_TEST_ID.tsv")

def test_no_call():
    result = classify_fixture("--")
    assert len(result) == 1
    assert result[0]["eye_color"] == "No call"
    assert result[0]["genotype"] is None
    assert result[0]["match_type"] == MatchType.NO_CALL.value
    assert result[0]["participant_id"] == "TEST_ID"
    assert result[0]["filename"] == "test.txt"
    assert result[0]["ref_count"] == 0
    assert result[0]["alt_count"] == 0
    # Cleanup
    os.remove("result_HERC2_TEST_ID.tsv")

def test_unexpected_allele_c():
    'Test handling of unexpected C allele (not in reference map).'
    result = classify_fixture("AC")
    assert len(result) == 1
    assert result[0]["eye_color"] == "Unknown"
    assert result[0]["genotype"] == "AC"
    assert result[0]["match_type"] == MatchType.VARIANT_CALL.value
    assert result[0]["participant_id"] == "TEST_ID"
    assert result[0]["filename"] == "test.txt"
    assert result[0]["ref_count"] == 1
    assert result[0]["alt_count"] == 1
    # Cleanup
    os.remove("result_HERC2_TEST_ID.tsv")
