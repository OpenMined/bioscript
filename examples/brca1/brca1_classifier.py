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

def generate_variant_calls(df: pd.DataFrame) -> list[VariantCall]:
    vcs: list[VariantCall] = []
    for _, row in df.iterrows():
        rsid_value: str | None = None
        if not pd.isna(row["rsid"]):
            candidate = str(row["rsid"]).strip()
            invalid_tokens = {".", "<NA>", "NA", "N/A"}
            if candidate and candidate.upper() not in invalid_tokens:
                rsid_value = candidate

        ref = str(row["ref"]).strip().upper()
        alt = str(row["alt"]).strip().upper()
        chrom = None
        if not pd.isna(row["chromosome"]):
            chrom_candidate = str(row["chromosome"]).strip()
            if chrom_candidate:
                chrom = chrom_candidate
        pos = None
        if not pd.isna(row["position"]):
            pos = int(row["position"])

        vc = VariantCall(
            rsid=rsid_value,
            ref=Alleles.from_letter(ref),
            alt=Alleles.from_not_letter(ref),
            chromosome=chrom,
            position=pos,
        )
        vc.ref_label = ref
        vc.alt_label = alt

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
        if not matches.all_matches:
            print("No variant matches were found.", flush=True)

        buckets = {
            "REFERENCE_CALL": [],
            "VARIANT_CALL": [],
            "NO_CALL": [],
        }

        def format_genotype(snp):
            from bioscript.types import Nucleotide

            symbols = []
            for allele in snp:
                if isinstance(allele, Nucleotide) and allele == Nucleotide.MISSING:
                    symbols.append('-')
                else:
                    symbols.append(allele.value)
            return ''.join(symbols)

        for match in matches.all_matches:
            call = match.variant_call
            chrom = call.chromosome or "?"
            pos = call.position if call.position is not None else "?"
            if call.rsid:
                if hasattr(call.rsid, "aliases"):
                    rsid_repr = ", ".join(sorted(call.rsid.aliases))
                else:
                    rsid_repr = str(call.rsid)
            else:
                rsid_repr = f"chr{chrom}:{pos}"

            ref = getattr(
                call,
                "ref_label",
                "/".join('-' if n.value == '.' else n.value for n in sorted(call.ref, key=lambda x: x.value)),
            )
            alt = getattr(
                call,
                "alt_label",
                "/".join('-' if n.value == '.' else n.value for n in sorted(call.alt, key=lambda x: x.value)),
            )
            genotype = format_genotype(match.snp)
            line = (
                f"{match.match_type.name:<14} {rsid_repr} chr{chrom}:{pos} "
                f"ref={ref} alt={alt} genotype={genotype}"
            )
            buckets.setdefault(match.match_type.name, []).append(line)

        def emit(header: str, lines: list[str]):
            if not lines:
                return
            print(header, flush=True)
            for idx, line in enumerate(lines, 1):
                print(f"  {idx:>3} {line}", flush=True)

        emit("Reference matches", buckets["REFERENCE_CALL"])
        emit("Variant matches", buckets["VARIANT_CALL"])
        emit("Diagnostics (no call)", buckets["NO_CALL"])

        return DiploidResult("DEBUG", "DEBUG")

__bioscript__ = {
    "variant_calls": get_vcs(),
    "classifier": BRCA1Classifier(),
    "name": "BRCA1",
}
