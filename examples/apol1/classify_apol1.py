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

from bioscript import write_tsv
from bioscript.types import MatchType

def _format_allele_label(allele):
    if isinstance(allele, Alleles):
        return ",".join(sorted(a.value for a in allele))
    if isinstance(allele, str):
        return allele
    if hasattr(allele, "__iter__"):
        return ",".join(sorted(str(a) for a in allele))
    return str(allele)

class APOL1Classifier(GenotypeClassifier):
    def classify(self, matches) -> list[dict[str, object]]:
        g2_match = matches.get(rs71785313)
        site1_match = matches.get(rs73885319)
        site2_match = matches.get(rs60910145)

        variant_matches = [
            ("rs71785313", rs71785313, g2_match),
            ("rs73885319", rs73885319, site1_match),
            ("rs60910145", rs60910145, site2_match),
        ]

        if not any(match is not None for _, _, match in variant_matches):
            diploid_result = DiploidResult(MISSING, MISSING)
        else:
            d_count = g2_match.alt_count if g2_match else 0
            site1_variants = site1_match.alt_count if site1_match else 0
            site2_variants = site2_match.alt_count if site2_match else 0

            has_g1 = site1_variants > 0 and site2_variants > 0
            g1_total = site1_variants + site2_variants if has_g1 else 0

            if d_count == 2:
                diploid_result = DiploidResult(APOL1Genotypes.G2, APOL1Genotypes.G2)
            elif d_count == 1:
                if g1_total >= 2:
                    diploid_result = DiploidResult(APOL1Genotypes.G2, APOL1Genotypes.G1)
                else:
                    diploid_result = DiploidResult(APOL1Genotypes.G2, APOL1Genotypes.G0)
            else:
                if g1_total == 4:
                    diploid_result = DiploidResult(APOL1Genotypes.G1, APOL1Genotypes.G1)
                elif g1_total >= 2:
                    diploid_result = DiploidResult(APOL1Genotypes.G1, APOL1Genotypes.G0)
                else:
                    diploid_result = DiploidResult(APOL1Genotypes.G0, APOL1Genotypes.G0)

        apol1_status = str(diploid_result.sorted())

        report_rows = []
        for fallback_rsid, variant_call, match in variant_matches:
            if match and match.source_row:
                rsid = match.source_row.rsid
                chromosome = match.source_row.chromosome
                position = match.source_row.position
            else:
                aliases = getattr(getattr(variant_call, "rsid", None), "aliases", None)
                rsid = sorted(aliases)[0] if aliases else fallback_rsid
                chromosome = getattr(variant_call, "chromosome", None)
                position = getattr(variant_call, "position", None)

            if match:
                genotype = match.genotype_sorted
                raw_match_type = match.match_type.value
                if match.has_missing or raw_match_type == MatchType.NO_CALL.value:
                    match_type = MatchType.NO_CALL.value
                    ref_count = 0
                    alt_count = 0
                else:
                    match_type = raw_match_type
                    ref_count = match.ref_count
                    alt_count = match.alt_count
            else:
                genotype = None
                match_type = MatchType.NO_CALL.value
                ref_count = 0
                alt_count = 0

            ref_label = _format_allele_label(variant_call.ref)
            alt_label = _format_allele_label(variant_call.alt)

            report_rows.append(
                {
                    "participant_id": self.participant_id,
                    "filename": self.filename,
                    "rsid": rsid,
                    "chromosome": chromosome,
                    "position": position,
                    "ref": ref_label,
                    "alt": alt_label,
                    "genotype": genotype,
                    "match_type": match_type,
                    "ref_count": ref_count,
                    "alt_count": alt_count,
                    "apol1_status": apol1_status,
                }
            )

        write_tsv(f"{self.output_basename}.tsv", report_rows)
        return report_rows

__bioscript__ = {
    "variant_calls": [rs73885319, rs60910145, rs71785313],
    "classifier": APOL1Classifier,
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

import os
from bioscript.types import MatchType

OUTPUT_FILE = "result_APOL1_TEST_ID.tsv"

def cleanup_output():
    if os.path.exists(OUTPUT_FILE):
        os.remove(OUTPUT_FILE)

def classify_fixture(genotypes):
    cleanup_output()
    variants = fixture(genotypes)
    matches = MatchList([rs73885319, rs60910145, rs71785313]).match_rows(variants)
    classifier = APOL1Classifier(participant_id="TEST_ID", name="APOL1", filename="test.txt")
    result = classifier(matches)
    assert isinstance(result, list)
    assert len(result) == 3
    rows = {row["rsid"]: row for row in result}
    assert set(rows.keys()) == {"rs73885319", "rs60910145", "rs71785313"}
    for row in rows.values():
        assert row["participant_id"] == "TEST_ID"
        assert row["filename"] == "test.txt"
        assert "ref" in row
        assert "alt" in row
        assert "ref_count" in row
        assert "alt_count" in row
    return rows

def test_g0_homozygous():
    rows = classify_fixture(["AA", "TT", "II"])
    assert rows["rs71785313"]["apol1_status"] == "G0/G0"
    assert rows["rs71785313"]["genotype"] == "II"
    assert rows["rs71785313"]["match_type"] == MatchType.NO_CALL.value
    assert rows["rs71785313"]["ref_count"] == 0
    assert rows["rs71785313"]["alt_count"] == 0

    assert rows["rs73885319"]["apol1_status"] == "G0/G0"
    assert rows["rs73885319"]["genotype"] == "AA"
    assert rows["rs73885319"]["match_type"] == MatchType.REFERENCE_CALL.value
    assert rows["rs73885319"]["ref_count"] == 2
    assert rows["rs73885319"]["alt_count"] == 0

    assert rows["rs60910145"]["apol1_status"] == "G0/G0"
    assert rows["rs60910145"]["genotype"] == "TT"
    assert rows["rs60910145"]["match_type"] == MatchType.REFERENCE_CALL.value
    assert rows["rs60910145"]["ref_count"] == 2
    assert rows["rs60910145"]["alt_count"] == 0

    cleanup_output()

def test_all_apol1_status_combinations():
    cases = [
        ("G0/G0", ["AA", "TT", "II"], {"rs73885319": "AA", "rs60910145": "TT", "rs71785313": "II"}),
        ("G1/G0", ["AG", "TC", "II"], {"rs73885319": "AG", "rs60910145": "CT", "rs71785313": "II"}),
        ("G1/G1", ["GG", "CC", "II"], {"rs73885319": "GG", "rs60910145": "CC", "rs71785313": "II"}),
        ("G2/G0", ["AA", "TT", "ID"], {"rs73885319": "AA", "rs60910145": "TT", "rs71785313": "DI"}),
        ("G2/G1", ["AG", "TC", "ID"], {"rs73885319": "AG", "rs60910145": "CT", "rs71785313": "DI"}),
        ("G2/G2", ["AA", "TT", "DD"], {"rs73885319": "AA", "rs60910145": "TT", "rs71785313": "DD"}),
    ]

    for expected_status, genotypes, expected_genotypes in cases:
        rows = classify_fixture(genotypes)

        assert set(rows.keys()) == {"rs73885319", "rs60910145", "rs71785313"}
        assert all(row["apol1_status"] == expected_status for row in rows.values())

        for rsid, expected_genotype in expected_genotypes.items():
            assert rows[rsid]["genotype"] == expected_genotype

        # rs71785313 is present in every test case and currently emitted as No call by the matcher.
        assert rows["rs71785313"]["match_type"] == MatchType.NO_CALL.value
        assert rows["rs71785313"]["ref_count"] == 0
        assert rows["rs71785313"]["alt_count"] == 0

        # rs73885319 and rs60910145 must still be callable and typed correctly.
        assert rows["rs73885319"]["match_type"] in {
            MatchType.REFERENCE_CALL.value,
            MatchType.VARIANT_CALL.value,
        }
        assert rows["rs60910145"]["match_type"] in {
            MatchType.REFERENCE_CALL.value,
            MatchType.VARIANT_CALL.value,
        }

    cleanup_output()

def test_g1_homozygous():
    rows = classify_fixture(["GG", "CC", "II"])
    assert all(row["apol1_status"] == "G1/G1" for row in rows.values())

    assert rows["rs73885319"]["genotype"] == "GG"
    assert rows["rs73885319"]["match_type"] == MatchType.VARIANT_CALL.value
    assert rows["rs73885319"]["ref_count"] == 0
    assert rows["rs73885319"]["alt_count"] == 2

    assert rows["rs60910145"]["genotype"] == "CC"
    assert rows["rs60910145"]["match_type"] == MatchType.VARIANT_CALL.value
    assert rows["rs60910145"]["ref_count"] == 0
    assert rows["rs60910145"]["alt_count"] == 2

    assert rows["rs71785313"]["genotype"] == "II"
    assert rows["rs71785313"]["match_type"] == MatchType.NO_CALL.value
    assert rows["rs71785313"]["ref_count"] == 0
    assert rows["rs71785313"]["alt_count"] == 0

    cleanup_output()
