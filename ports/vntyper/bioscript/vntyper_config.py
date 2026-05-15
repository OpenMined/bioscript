"""Explicit VNtyper data and configuration used by the BioScript port."""

from __future__ import annotations

DEFAULT_KESTREL_JAR = "ports/vntyper/kestrel/kestrel.jar"
DEFAULT_MUC1_REFERENCE = "ports/vntyper/vntyper/reference/All_Pairwise_and_Self_Merged_MUC1_motifs_filtered.fa"

COORDINATE_SYSTEMS = {
    "GRCh37": {
        "chromosome": 1,
        "bam_region_coords": "155158000-155163000",
        "vntr_region_coords": "155160500-155162000",
    },
    "GRCh38": {
        "chromosome": 1,
        "bam_region_coords": "155184000-155194000",
        "vntr_region_coords": "155188000-155192500",
    },
}

ASSEMBLY_METADATA = {
    "hg19": {"coordinate_system": "GRCh37", "reference_source": "ucsc"},
    "hg38": {"coordinate_system": "GRCh38", "reference_source": "ucsc"},
    "GRCh37": {"coordinate_system": "GRCh37", "reference_source": "ncbi"},
    "GRCh38": {"coordinate_system": "GRCh38", "reference_source": "ncbi"},
    "hg19_ncbi": {"coordinate_system": "GRCh37", "reference_source": "ncbi"},
    "hg38_ncbi": {"coordinate_system": "GRCh38", "reference_source": "ncbi"},
    "hg19_ensembl": {"coordinate_system": "GRCh37", "reference_source": "ensembl"},
    "hg38_ensembl": {"coordinate_system": "GRCh38", "reference_source": "ensembl"},
}

ASSEMBLY_ALIASES = {name: name for name in ASSEMBLY_METADATA}

KNOWN_NCBI_ACCESSIONS = {
    "GRCh37": "NC_000001.10",
    "GRCh38": "NC_000001.11",
}

NATIVE_KESTREL_MAX_HAPLOTYPES = 2
NATIVE_KESTREL_MAX_SAVED_STATES = 2
NATIVE_KESTREL_MAX_BASES = 120
NATIVE_KESTREL_MIN_KMER_COUNT = 5

OPTIONAL_VALIDATION_DEFAULTS = {
    "advntr_enabled": False,
    "advntr_result_when_disabled": "none",
}

REPORT_SCHEMA_KEYS = [
    "sample_name",
    "version",
    "report_date",
    "metadata",
    "input_files",
    "coverage",
    "fastp",
    "algorithm_results",
    "screening_summary",
    "kestrel_variants",
    "advntr_variants",
    "cross_match_summary",
    "pipeline_log",
]

DEFAULT_KESTREL_CONFIG = {
    "confidence_assignment": {
        "depth_score_thresholds": {
            "low": 0.00469,
            "high": 0.00515,
        },
        "alt_depth_thresholds": {
            "low": 20,
            "mid_low": 21,
            "mid_high": 100,
        },
        "var_active_region_threshold": 200,
        "confidence_levels": {
            "low_precision": "Low_Precision",
            "high_precision": "High_Precision",
            "high_precision_star": "High_Precision*",
        },
    },
    "alt_filtering": {
        "gg_alt_value": "GG",
        "gg_depth_score_threshold": 0.00469,
        "exclude_alts": [],
    },
    "motif_filtering": {
        "use_uniform_filtering": False,
        "position_threshold": 60,
        "exclude_motifs_right": ["8", "9", "7", "6p", "6"],
        "alt_for_motif_right_gg": "GG",
        "motifs_for_alt_gg": [],
        "exclude_alts_combined": ["CCGCC", "CGGCG", "CGGCC"],
        "exclude_motifs_combined": ["6", "6p", "7"],
    },
    "flagging_rules": {
        "False_Positive_4bp_Insertion": "(REF == 'C') and (ALT == 'CGGCA')",
        "Low_Depth_Conserved_Motifs": "(Depth_Score < 0.4) and (Motif in ['1', '2', '3', '4', '6', '7', '8', '9'])",
    },
    "duplicate_flagging": {
        "enabled": False,
        "flag_name": "Potential_Duplicate",
        "group_by": ["REF", "ALT"],
        "sort_by": [
            {"column": "Depth_Score", "ascending": False},
            {"column": "Motifs", "ascending": True},
            {"column": "POS", "ascending": True},
        ],
    },
}

DEFAULT_REPORT_CONFIG = {
    "mean_vntr_coverage_threshold": 100,
    "algorithm_logic": {
        "kestrel": {
            "rules": [
                {
                    "conditions": {
                        "Confidence": {"operator": "in", "value": ["High_Precision", "High_Precision*"]},
                        "Flag": {"operator": "==", "value": "Not flagged"},
                    },
                    "result": "High_Precision",
                },
                {
                    "conditions": {
                        "Confidence": {"operator": "in", "value": ["Low_Precision"]},
                        "Flag": {"operator": "==", "value": "Not flagged"},
                    },
                    "result": "Low_Precision",
                },
                {
                    "conditions": {
                        "Confidence": {"operator": "in", "value": ["High_Precision", "High_Precision*"]},
                        "Flag": {"operator": "!=", "value": "Not flagged"},
                    },
                    "result": "High_Precision_flagged",
                },
                {
                    "conditions": {
                        "Confidence": {"operator": "in", "value": ["Low_Precision"]},
                        "Flag": {"operator": "!=", "value": "Not flagged"},
                    },
                    "result": "Low_Precision_flagged",
                },
            ],
            "default": "negative",
        },
        "advntr": {
            "rules": [
                {
                    "conditions": {
                        "VID": {"operator": "!=", "value": "Negative"},
                        "Flag": {"operator": "==", "value": "Not flagged"},
                    },
                    "result": "positive",
                },
                {
                    "conditions": {
                        "Flag": {"operator": "not in", "value": ["Not flagged", "Not applicable", "None"]},
                    },
                    "result": "positive flagged",
                },
            ],
            "default": "negative",
        },
    },
    "screening_summary_default": "The screening was negative (no valid Kestrel or adVNTR data).",
    "screening_summary_rules": [
        {
            "conditions": {
                "kestrel_result": "High_Precision",
                "advntr_result": "none",
                "quality_metrics_pass": True,
            },
            "message": "Kestrel detected a high-precision pathogenic variant.<br>Note: adVNTR genotyping was not performed.<br>It is recommended to perform adVNTR and validate the result using orthogonal methods (e.g., SNaPshot, long-read sequencing).",
        },
        {
            "conditions": {
                "kestrel_result": "High_Precision",
                "advntr_result": "none",
                "quality_metrics_pass": False,
            },
            "message": "Kestrel detected a high-precision pathogenic variant with quality metrics below threshold, and adVNTR genotyping was not performed.<br>Further validation using alternative methods (e.g., SNaPshot, long-read sequencing) is strongly recommended.",
        },
        {
            "conditions": {
                "kestrel_result": "High_Precision_flagged",
                "advntr_result": "none",
                "quality_metrics_pass": True,
            },
            "message": "Kestrel detected a high-precision pathogenic variant with a flagged result.<br>Note: adVNTR genotyping was not performed.<br>It is recommended to perform adVNTR and validate the finding using orthogonal methods (e.g., SNaPshot, long-read sequencing).",
        },
        {
            "conditions": {
                "kestrel_result": "Low_Precision",
                "advntr_result": "none",
                "quality_metrics_pass": True,
            },
            "message": "Kestrel detected a pathogenic variant with low precision.<br>Note: adVNTR genotyping was not performed.<br>It is recommended to perform adVNTR and validate the result using alternative methods (e.g., SNaPshot, long-read sequencing).",
        },
        {
            "conditions": {
                "kestrel_result": "negative",
                "advntr_result": "none",
                "quality_metrics_pass": True,
            },
            "message": "No variant detected.<br>Note: adVNTR genotyping was not performed.",
        },
    ],
}
