import pandas as pd
from bioscript import optional_int, optional_str, write_tsv
from bioscript.classifier import GenotypeClassifier
from bioscript.types import VariantCall

def generate_variant_calls(df: pd.DataFrame) -> list[VariantCall]:
    """Generate VariantCall objects from ClinVar DataFrame."""
    vcs: list[VariantCall] = []
    for _, row in df.iterrows():
        vcs.append(
            VariantCall(
                rsid=optional_str(row["rsid"]),
                ref=optional_str(row["ref"]),
                alt=optional_str(row["alt"]),
                chromosome=optional_str(row["chromosome"]),
                position=optional_int(row["position"]),
                gene=optional_str(row.get("gene"), upper=True),
            )
        )
    return vcs

def get_vcs() -> list[VariantCall]:
    """Load BRCA1 and BRCA2 variant calls from ClinVar TSV files."""
    dfs = [pd.read_csv(f, sep="\t") for f in ["brca1_clinvar.tsv", "brca2_clinvar.tsv"]]
    df = pd.concat(dfs, ignore_index=True)
    print(f"Loaded {len(df)} variants from BRCA1 and BRCA2")
    return generate_variant_calls(df)

class BRCAClassifier(GenotypeClassifier):
    def classify(self, matches):
        """Classify BRCA variants and write results to TSV files."""
        if not matches.all_matches:
            print("No variant matches were found.", flush=True)

        # Get categorized matches as report rows
        ref_rows, var_rows, no_rows = matches.categorize_report_rows(
            self.participant_id, self.filename
        )

        if self.debug:
            write_tsv(f"{self.output_basename}_ref.tsv", ref_rows)
            write_tsv(f"{self.output_basename}_no.tsv", no_rows)

        write_tsv(f"{self.output_basename}.tsv", var_rows)
        
        # Return variant rows for testing
        return var_rows

__bioscript__ = {
    "variant_calls": get_vcs,
    "classifier": BRCAClassifier,
    "name": "BRCA",
}

from bioscript import VariantFixture
from bioscript.types import MatchList
import os

# Create test fixtures for BRCA1 and BRCA2 variants
fixture = VariantFixture(
    [
        {"rsid": "rs80357336", "chromosome": "17", "position": 43045711},
        {"rsid": "rs886040303", "chromosome": "17", "position": 43045728},
        {"rsid": "rs397509295", "chromosome": "17", "position": 43045729},
        {"rsid": "rs80358650", "chromosome": "13", "position": 32316463},
        {"rsid": "rs397507571", "chromosome": "13", "position": 32316470},
        {"rsid": "rs80358622", "chromosome": "13", "position": 32316497},
    ],
    assembly="GRCh38",
)

def test_brca1_heterozygous_variants():
    """Test detection of heterozygous BRCA1 variants."""
    # Create test data with heterozygous variants (one alt allele)
    variants = fixture(["GC", "GA", "GT", "GG", "GG", "GG"])
    
    # Create mini variant call list for testing
    test_vcs = [
        VariantCall(rsid="rs80357336", ref="G", alt="C", chromosome="17", position=43045711, gene="BRCA1"),
        VariantCall(rsid="rs886040303", ref="G", alt="A", chromosome="17", position=43045728, gene="BRCA1"),
        VariantCall(rsid="rs397509295", ref="G", alt="T", chromosome="17", position=43045729, gene="BRCA1"),
    ]
    
    matches = MatchList(variant_calls=test_vcs).match_rows(variants)
    classifier = BRCAClassifier(participant_id="TEST_HET", name="BRCA", filename="test.txt")
    result = classifier(matches)
    
    assert len(result) == 3, f"Expected 3 variant rows, got {len(result)}"
    assert all(row["gene"] == "BRCA1" for row in result), "All variants should be BRCA1"
    assert all(row["match_type"] == "VARIANT_CALL" for row in result), "All should be variant calls"
    
    # Cleanup output file
    os.remove("result_BRCA_TEST_HET.tsv")

def test_brca2_homozygous_variant():
    """Test detection of homozygous BRCA2 variant."""
    # Create test data with one homozygous variant (two alt alleles)
    variants = fixture(["GG", "GG", "GG", "AA", "GG", "GG"])
    
    test_vcs = [
        VariantCall(rsid="rs80358650", ref="G", alt="A", chromosome="13", position=32316463, gene="BRCA2"),
    ]

    matches = MatchList(variant_calls=test_vcs).match_rows(variants)
    classifier = BRCAClassifier(participant_id="TEST_HOM", name="BRCA", filename="test.txt")
    result = classifier(matches)
    
    assert len(result) == 1, f"Expected 1 variant row, got {len(result)}"
    assert result[0]["gene"] == "BRCA2", "Variant should be BRCA2"
    assert result[0]["genotype"] == "AA", "Should be homozygous AA"
    
    # Cleanup output file
    os.remove("result_BRCA_TEST_HOM.tsv")

def test_no_variants():
    """Test classifier with no matching variants."""
    # All reference genotypes
    variants = fixture(["GG", "GG", "GG", "GG", "GG", "GG"])
    
    test_vcs = [
        VariantCall(rsid="rs80357336", ref="G", alt="C", chromosome="17", position=43045711, gene="BRCA1"),
    ]
    
    matches = MatchList(variant_calls=test_vcs).match_rows(variants)
    classifier = BRCAClassifier(participant_id="TEST_REF", name="BRCA", filename="test.txt")
    result = classifier(matches)
    
    assert len(result) == 0, f"Expected 0 variant rows, got {len(result)}"
    
    # Cleanup output file
    os.remove("result_BRCA_TEST_REF.tsv")
