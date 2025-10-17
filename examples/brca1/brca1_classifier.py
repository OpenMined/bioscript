from bioscript import AlleleCounter
from bioscript.classifier import DiploidResult, GenotypeClassifier, GenotypeEnum
from bioscript.types import Alleles, VariantCall

import pandas as pd

def filter_snvs(df: pd.DataFrame) -> pd.DataFrame:
    """
    Return only rows where clnvc == 'single_nucleotide_variant' (case-insensitive).
    """
    mask = df["clnvc"].str.lower() == "single_nucleotide_variant"
    return df[mask].reset_index(drop=True)

def generate_variant_calls(df: pd.DataFrame) -> list[str]:
    vcs = []
    for _, row in df.iterrows():
        rsid = str(row["rsid"]).strip()
        ref = str(row["ref"]).strip().upper()
        alt = str(row["alt"]).strip().upper()

        # Build readable variant call
        vc = VariantCall(rsid=rsid, ref=Alleles.from_letter(ref), alt=Alleles.from_not_letter(ref))
        vcs.append(vc)

    return vcs

def get_vcs():
    # Path to your TSV file
    tsv_path = "brca1_clinvar.tsv"
    
    # Load the TSV file
    df = pd.read_csv(
        tsv_path,
        sep="\t",
        dtype={
            "rsid": "string",
            "gene": "string",
            "chromosome": "string",
            "position": "Int64",
            "ref": "string",
            "alt": "string",
            "clnrevstat": "string",
            "clnsig": "string",
            "clnvc": "string",
        }
    )

    # Example usage:
    df_snvs = filter_snvs(df)
    vcs = generate_variant_calls(df_snvs)
    return vcs

class BRCA1Classifier(GenotypeClassifier):
    def classify(self, matches) -> DiploidResult:
        print(type(matches))
        print(len(matches.all_matches))
        print(len(matches.variant_matches))

__bioscript__ = {
    "variant_calls": get_vcs(),
    "classifier": BRCA1Classifier(),
    "name": "BRCA1",
}
