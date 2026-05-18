"""Minimal VNtyper region/reference helpers for the BioScript port."""

from __future__ import annotations

import re

try:
    from . import vntyper_config
except ImportError:
    import vntyper_config

COORDINATE_SYSTEMS = vntyper_config.COORDINATE_SYSTEMS
ASSEMBLY_METADATA = vntyper_config.ASSEMBLY_METADATA
ASSEMBLY_ALIASES = vntyper_config.ASSEMBLY_ALIASES
KNOWN_NCBI_ACCESSIONS = vntyper_config.KNOWN_NCBI_ACCESSIONS


def normalize_assembly_name(user_input: str) -> str:
    if user_input not in ASSEMBLY_ALIASES:
        supported = ", ".join(sorted(ASSEMBLY_ALIASES))
        raise ValueError(f"Unknown assembly '{user_input}'. Supported assemblies: {supported}")
    return ASSEMBLY_ALIASES[user_input]


def get_coordinate_system(assembly_name: str) -> str:
    canonical = normalize_assembly_name(assembly_name)
    return ASSEMBLY_METADATA[canonical]["coordinate_system"]


def get_reference_source(assembly_name: str) -> str:
    canonical = normalize_assembly_name(assembly_name)
    return ASSEMBLY_METADATA[canonical]["reference_source"]


def get_coordinates(assembly_name: str, region_type: str) -> str:
    coordinate_system = get_coordinate_system(assembly_name)
    coordinates = COORDINATE_SYSTEMS[coordinate_system].get(region_type)
    if coordinates is None:
        raise ValueError(f"Unknown region type '{region_type}' for assembly '{assembly_name}'")
    return coordinates


def detect_naming_convention(contig_names: list[str]) -> str:
    if not contig_names:
        return "unknown"

    counts = {"ucsc": 0, "ncbi": 0, "ensembl": 0}
    for name in contig_names:
        if re.match(r"^chr[0-9XYM]+$", name, re.IGNORECASE):
            counts["ucsc"] += 1
        elif re.match(r"^NC_\d{6}\.\d+$", name):
            counts["ncbi"] += 1
        elif re.match(r"^([0-9]+|X|Y|MT?)$", name, re.IGNORECASE):
            counts["ensembl"] += 1

    total = len(contig_names)
    for convention, count in counts.items():
        if count / total >= 0.5:
            return convention
    return "unknown"


def chromosome_name(chromosome_number: int, assembly_name: str, convention: str | None = None) -> str:
    coordinate_system = get_coordinate_system(assembly_name)
    source = convention or get_reference_source(assembly_name)
    if source == "ucsc":
        return f"chr{chromosome_number}"
    if source == "ensembl":
        return str(chromosome_number)
    if source == "ncbi" and chromosome_number == 1:
        return KNOWN_NCBI_ACCESSIONS[coordinate_system]
    raise ValueError(f"Unsupported chromosome source '{source}' for chromosome {chromosome_number}")


def validate_chromosome_name(name: str) -> bool:
    if not name:
        return False
    patterns = [
        r"^chr[0-9]+$",
        r"^chr[XYM]$",
        r"^[0-9]+$",
        r"^[XYMT]+$",
        r"^NC_\d{6}\.\d+$",
    ]
    return any(re.match(pattern, name, re.IGNORECASE) for pattern in patterns)


def build_region_string(chromosome: str, coordinates: str) -> str:
    if not validate_chromosome_name(chromosome):
        raise ValueError(f"Invalid chromosome name: '{chromosome}'")
    if "-" not in coordinates:
        raise ValueError(f"Invalid coordinate format: '{coordinates}'")
    start, end = coordinates.split("-", maxsplit=1)
    start_i = int(start)
    end_i = int(end)
    if end_i < start_i:
        raise ValueError(f"Invalid coordinate range: '{coordinates}'")
    return f"{chromosome}:{start_i}-{end_i}"


def region_string(assembly_name: str, region_type: str, convention: str | None = None) -> str:
    coordinate_system = get_coordinate_system(assembly_name)
    chromosome_number = COORDINATE_SYSTEMS[coordinate_system]["chromosome"]
    chromosome = chromosome_name(chromosome_number, assembly_name, convention=convention)
    coordinates = get_coordinates(assembly_name, region_type)
    return build_region_string(chromosome, coordinates)
