import importlib.util
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]
REGIONS_PATH = ROOT / "ports" / "vntyper" / "bioscript" / "vntyper_regions.py"


spec = importlib.util.spec_from_file_location("vntyper_regions", REGIONS_PATH)
vntyper_regions = importlib.util.module_from_spec(spec)
spec.loader.exec_module(vntyper_regions)


class VntyperRegionTests(unittest.TestCase):
    def test_reference_assembly_aliases_match_upstream_coordinates(self):
        self.assertEqual(vntyper_regions.get_coordinate_system("hg19"), "GRCh37")
        self.assertEqual(vntyper_regions.get_coordinate_system("hg38"), "GRCh38")
        self.assertEqual(
            vntyper_regions.get_coordinates("hg19", "bam_region_coords"),
            "155158000-155163000",
        )
        self.assertEqual(
            vntyper_regions.get_coordinates("hg38", "vntr_region_coords"),
            "155188000-155192500",
        )

    def test_region_strings_follow_reference_source_naming(self):
        self.assertEqual(
            vntyper_regions.region_string("hg19", "bam_region_coords"),
            "chr1:155158000-155163000",
        )
        self.assertEqual(
            vntyper_regions.region_string("hg19_ncbi", "bam_region_coords"),
            "NC_000001.10:155158000-155163000",
        )
        self.assertEqual(
            vntyper_regions.region_string("hg38_ensembl", "vntr_region_coords"),
            "1:155188000-155192500",
        )

    def test_detect_naming_convention_matches_upstream_patterns(self):
        self.assertEqual(vntyper_regions.detect_naming_convention(["chr1", "chr2", "chrX"]), "ucsc")
        self.assertEqual(vntyper_regions.detect_naming_convention(["1", "2", "X"]), "ensembl")
        self.assertEqual(
            vntyper_regions.detect_naming_convention(["NC_000001.10", "NC_000002.11"]),
            "ncbi",
        )
        self.assertEqual(vntyper_regions.detect_naming_convention([]), "unknown")

    def test_rejects_unknown_assembly_and_invalid_coordinates(self):
        with self.assertRaises(ValueError):
            vntyper_regions.normalize_assembly_name("mm10")
        with self.assertRaises(ValueError):
            vntyper_regions.build_region_string("chr1", "10-1")
        with self.assertRaises(ValueError):
            vntyper_regions.build_region_string("bad_chr", "1-10")


if __name__ == "__main__":
    unittest.main()
