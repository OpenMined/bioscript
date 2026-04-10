import pandas as pd
from bioscript import optional_int, optional_str, write_tsv
from bioscript.classifier import GenotypeClassifier
from bioscript.types import VariantCall
from bioscript import assets_dir

ASSETS_DIR = assets_dir()
CLINVAR_TSV = 'thalassemia_clinvar.tsv'
RESULT_HEADERS = [
    'participant_id',
    'filename',
    'gene',
    'rsid',
    'chromosome',
    'position',
    'genotype',
    'ref',
    'alt',
    'variant_type',
    'match_type',
    'ref_count',
    'alt_count',
]

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
    """Load thalassemia-associated variant calls from a ClinVar TSV file."""
    df = pd.read_csv(ASSETS_DIR / CLINVAR_TSV, sep='	')
    print(f'Loaded {len(df)} variants from {CLINVAR_TSV}')
    return generate_variant_calls(df)

class ThalassemiaClassifier(GenotypeClassifier):
    def classify(self, matches):
        """Classify thalassemia-associated variants and write results to TSV files."""
        if not matches.all_matches:
            print('No variant matches were found.', flush=True)

        # Get categorized matches as report rows
        ref_rows, var_rows, no_rows = matches.categorize_report_rows(
            self.participant_id, self.filename
        )

        if self.debug:
            write_tsv(f'{self.output_basename}_ref.tsv', ref_rows)
            write_tsv(f'{self.output_basename}_no.tsv', no_rows)

        write_tsv(f'{self.output_basename}.tsv', var_rows, headers=RESULT_HEADERS)

        # Return variant rows for testing
        return var_rows

__bioscript__ = {
    'variant_calls': get_vcs,
    'classifier': ThalassemiaClassifier,
    'name': 'THALASSEMIA',
}

from bioscript import VariantFixture
from bioscript.types import MatchList
import os

# Create test fixtures for thalassemia-associated HBB variants (subset from thalassemia_clinvar.tsv)
fixture = VariantFixture(
    [
        {'rsid': 'rs33985472', 'chromosome': '11', 'position': 5225485},
        {'rsid': 'rs63751128', 'chromosome': '11', 'position': 5225487},
        {'rsid': 'rs33978907', 'chromosome': '11', 'position': 5225488},
        {'rsid': 'rs34809925', 'chromosome': '11', 'position': 5225592},
        {'rsid': 'rs35117167', 'chromosome': '11', 'position': 5225605},
        {'rsid': 'rs33971634', 'chromosome': '11', 'position': 5225660},
    ],
    assembly='GRCh38',
)

def test_thalassemia_heterozygous_variants():
    """Test detection of heterozygous thalassemia-associated variants."""
    variants = fixture(['TC', 'TC', 'AG', 'GG', 'TT', 'GG'])

    # Create mini variant call list for testing
    test_vcs = [
        VariantCall(rsid='rs33985472', ref='T', alt='C', chromosome='11', position=5225485, gene='HBB'),
        VariantCall(rsid='rs63751128', ref='T', alt='C', chromosome='11', position=5225487, gene='HBB'),
        VariantCall(rsid='rs33978907', ref='A', alt='G', chromosome='11', position=5225488, gene='HBB'),
    ]

    matches = MatchList(variant_calls=test_vcs).match_rows(variants)
    classifier = ThalassemiaClassifier(participant_id='TEST_HET', name='THALASSEMIA', filename='test.txt')
    result = classifier(matches)

    assert len(result) == 3, f'Expected 3 variant rows, got {len(result)}'
    assert all(row['gene'] == 'HBB' for row in result), 'All variants should be HBB'
    assert all(row['match_type'] == 'VARIANT_CALL' for row in result), 'All should be variant calls'

    # Cleanup output file
    os.remove('result_THALASSEMIA_TEST_HET.tsv')

def test_thalassemia_homozygous_variant():
    """Test detection of a homozygous thalassemia-associated variant."""
    variants = fixture(['TT', 'TT', 'AA', 'CC', 'TT', 'GG'])

    test_vcs = [
        VariantCall(rsid='rs34809925', ref='G', alt='C', chromosome='11', position=5225592, gene='HBB'),
    ]

    matches = MatchList(variant_calls=test_vcs).match_rows(variants)
    classifier = ThalassemiaClassifier(participant_id='TEST_HOM', name='THALASSEMIA', filename='test.txt')
    result = classifier(matches)

    assert len(result) == 1, f'Expected 1 variant row, got {len(result)}'
    assert result[0]['gene'] == 'HBB', 'Variant should be HBB'
    assert result[0]['genotype'] == 'CC', 'Should be homozygous CC'

    # Cleanup output file
    os.remove('result_THALASSEMIA_TEST_HOM.tsv')

def test_no_variants():
    """Test classifier with no matching variants."""
    variants = fixture(['TT', 'TT', 'AA', 'GG', 'TT', 'GG'])

    test_vcs = [
        VariantCall(rsid='rs33985472', ref='T', alt='C', chromosome='11', position=5225485, gene='HBB'),
    ]

    matches = MatchList(variant_calls=test_vcs).match_rows(variants)
    classifier = ThalassemiaClassifier(participant_id='TEST_REF', name='THALASSEMIA', filename='test.txt')
    result = classifier(matches)

    assert len(result) == 0, f'Expected 0 variant rows, got {len(result)}'

    # Cleanup output file
    os.remove('result_THALASSEMIA_TEST_REF.tsv')
