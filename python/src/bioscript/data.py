"""
BioScript Data Utilities

Fetch and manage test genomic datasets for bioinformatics workflows.
"""

from __future__ import annotations

import os
import zipfile
from collections.abc import Iterator
from pathlib import Path
from urllib.request import urlretrieve

from .types import VariantRow

# Sample data registry
SAMPLES = {
    "23andme_v4": {
        "url": "https://github.com/OpenMined/biovault-data/raw/main/snp/23andme_genome_v4_Full.zip",
        "description": "23andMe Genome v4 Full - Complete SNP genotyping dataset",
        "filename": "23andme_genome_v4_Full.zip",
    }
}


def fetch_sample(sample_name: str, output_dir: str | None = None, force: bool = False) -> Path:
    """
    Download and extract a sample genomic dataset.

    Args:
        sample_name: Name of the sample dataset (e.g., "23andme_v4")
        output_dir: Directory to download/extract to (default: current directory)
        force: Re-download even if file already exists

    Returns:
        Path to the main data file (e.g., genome_*.txt)

    Raises:
        ValueError: If sample_name is not recognized

    Example:
        >>> from bioscript.data import fetch_sample
        >>> data_file = fetch_sample("23andme_v4")
        >>> print(f"Data file: {data_file}")
        >>> # Use with pandas
        >>> import pandas as pd
        >>> df = pd.read_csv(data_file, sep='\\t', comment='#')
    """
    if sample_name not in SAMPLES:
        available = ", ".join(SAMPLES.keys())
        raise ValueError(f"Unknown sample '{sample_name}'. Available samples: {available}")

    sample = SAMPLES[sample_name]

    # Determine output directory
    if output_dir is None:
        output_dir = os.getcwd()
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    # Download file
    zip_path = output_path / sample["filename"]

    if not zip_path.exists() or force:
        print(f"📥 Downloading {sample_name}...")
        print(f"   URL: {sample['url']}")
        print(f"   Destination: {zip_path}")

        urlretrieve(sample["url"], zip_path)
        print(f"✅ Downloaded: {zip_path.name}")
    else:
        print(f"✅ Using cached file: {zip_path.name}")

    # Extract zip
    extract_dir = output_path / sample_name

    if not extract_dir.exists() or force:
        print(f"📦 Extracting {zip_path.name}...")

        with zipfile.ZipFile(zip_path, "r") as zip_ref:
            zip_ref.extractall(extract_dir)

        # List extracted files
        extracted_files = list(extract_dir.rglob("*"))
        file_count = len([f for f in extracted_files if f.is_file()])
        print(f"✅ Extracted {file_count} file(s) to: {extract_dir}")

        # Show extracted files
        for file in sorted(extracted_files):
            if file.is_file():
                size_mb = file.stat().st_size / (1024 * 1024)
                print(f"   - {file.name} ({size_mb:.2f} MB)")
    else:
        print(f"✅ Using cached extraction: {extract_dir}")

    # Find and return the main data file (first .txt file)
    data_files = list(extract_dir.glob("*.txt"))
    if data_files:
        return data_files[0]

    # Fallback: return any file in the directory
    all_files = [f for f in extract_dir.rglob("*") if f.is_file()]
    if all_files:
        return all_files[0]

    # If no files found, return the directory
    return extract_dir


def list_samples() -> None:
    """
    List all available sample datasets.

    Example:
        >>> from bioscript.data import list_samples
        >>> list_samples()
    """
    print("📊 Available Sample Datasets:\n")

    for name, info in SAMPLES.items():
        print(f"  {name}")
        print(f"    Description: {info['description']}")
        print(f"    URL: {info['url']}")
        print()


def get_sample_info(sample_name: str) -> dict:
    """
    Get metadata about a sample dataset.

    Args:
        sample_name: Name of the sample dataset

    Returns:
        Dictionary with sample metadata

    Raises:
        ValueError: If sample_name is not recognized
    """
    if sample_name not in SAMPLES:
        available = ", ".join(SAMPLES.keys())
        raise ValueError(f"Unknown sample '{sample_name}'. Available samples: {available}")

    return SAMPLES[sample_name].copy()


def create_test_variants(
    variants: list[dict[str, str | int | float]],
    output_file: str | Path | None = None,
) -> Iterator[VariantRow] | Path:
    """
    Create test variant data for testing purposes.

    Args:
        variants: List of variant dictionaries with keys:
            - rsid (required): SNP identifier (e.g., "rs123")
            - chromosome (required): Chromosome number/name (e.g., "1", "X")
            - position (required): Position on chromosome (int)
            - genotype (required): Genotype string (e.g., "AA", "AT")
            - gs (optional): GenCall score (float)
            - baf (optional): B Allele Frequency (float)
            - lrr (optional): Log R Ratio (float)
        output_file: If provided, write to this file path and return the path.
                    If None, return an iterator of VariantRow objects directly.

    Returns:
        If output_file is None: Iterator of VariantRow objects
        If output_file is provided: Path to the created file

    Examples:
        >>> # In-memory variant rows (for testing)
        >>> from bioscript.data import create_test_variants
        >>> variants = create_test_variants([
        ...     {"rsid": "rs123", "chromosome": "1", "position": 1000, "genotype": "AA"},
        ...     {"rsid": "rs456", "chromosome": "2", "position": 2000, "genotype": "AT"},
        ... ])
        >>> for variant in variants:
        ...     print(variant.rsid, variant.genotype)

        >>> # Write to file and use with load_variants_tsv
        >>> from bioscript import load_variants_tsv
        >>> test_file = create_test_variants([
        ...     {"rsid": "rs123", "chromosome": "1", "position": 1000, "genotype": "AA"},
        ... ], output_file="test_data.txt")
        >>> variants = load_variants_tsv(test_file)

        >>> # Use temporary file
        >>> import tempfile
        >>> with tempfile.NamedTemporaryFile(mode='w', suffix='.txt', delete=False) as f:
        ...     test_file = create_test_variants(variants_data, output_file=f.name)
    """
    # Validate required fields
    required_fields = {"rsid", "chromosome", "position", "genotype"}
    for i, variant in enumerate(variants):
        missing = required_fields - set(variant.keys())
        if missing:
            raise ValueError(f"Variant at index {i} missing required fields: {', '.join(missing)}")

    if output_file is None:
        # Return iterator of VariantRow objects directly
        def variant_iterator():
            for v in variants:
                yield VariantRow(
                    rsid=str(v["rsid"]),
                    chromosome=str(v["chromosome"]),
                    position=int(v["position"]),
                    genotype=str(v["genotype"]),
                    gs=float(v["gs"]) if "gs" in v else None,
                    baf=float(v["baf"]) if "baf" in v else None,
                    lrr=float(v["lrr"]) if "lrr" in v else None,
                )

        return variant_iterator()

    # Write to file
    output_path = Path(output_file)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    with open(output_path, "w") as f:
        # Write header
        header = ["rsid", "chromosome", "position", "genotype"]
        optional_fields = []
        if any("gs" in v for v in variants):
            optional_fields.append("gs")
        if any("baf" in v for v in variants):
            optional_fields.append("baf")
        if any("lrr" in v for v in variants):
            optional_fields.append("lrr")
        header.extend(optional_fields)

        f.write("# " + "\t".join(header) + "\n")

        # Write data rows
        for v in variants:
            row = [
                str(v["rsid"]),
                str(v["chromosome"]),
                str(v["position"]),
                str(v["genotype"]),
            ]
            for field in optional_fields:
                row.append(str(v.get(field, "")))
            f.write("\t".join(row) + "\n")

    return output_path


class GenotypeGenerator:
    """
    Helper for creating test variants with fixed positions but varying genotypes.

    This is useful for testing different genotype combinations on the same variant positions.

    Args:
        variant_templates: List of variant dictionaries WITHOUT genotype field.
                          Must include: rsid, chromosome, position
                          Optional: gs, baf, lrr

    Examples:
        >>> from bioscript import GenotypeGenerator
        >>> gen = GenotypeGenerator([
        ...     {"rsid": "rs123", "chromosome": "1", "position": 1000},
        ...     {"rsid": "rs456", "chromosome": "2", "position": 2000},
        ... ])
        >>> # Generate variants with different genotypes
        >>> variants1 = gen(["AA", "TT"])  # First scenario
        >>> variants2 = gen(["AT", "TC"])  # Second scenario
        >>> for v in variants1:
        ...     print(v.rsid, v.genotype)

        >>> # Can also write to file
        >>> test_file = gen(["AA", "TT"], output_file="test.txt")

        >>> # Access original templates
        >>> gen.templates
    """

    def __init__(self, variant_templates: list[dict[str, str | int | float]]):
        """Initialize generator with variant templates (without genotypes)."""
        # Validate templates
        required_fields = {"rsid", "chromosome", "position"}
        for i, template in enumerate(variant_templates):
            missing = required_fields - set(template.keys())
            if missing:
                raise ValueError(
                    f"Template at index {i} missing required fields: {', '.join(missing)}"
                )
            if "genotype" in template:
                raise ValueError(f"Template at index {i} should not include 'genotype' field")

        self.templates = variant_templates

    def __call__(
        self, genotypes: list[str], output_file: str | Path | None = None
    ) -> Iterator[VariantRow] | Path:
        """
        Generate variants with the given genotypes.

        Args:
            genotypes: List of genotype strings, must match length of templates
            output_file: If provided, write to file and return path.
                        If None, return iterator of VariantRow objects.

        Returns:
            If output_file is None: Iterator of VariantRow objects
            If output_file is provided: Path to the created file

        Raises:
            ValueError: If genotypes list length doesn't match templates length
        """
        if len(genotypes) != len(self.templates):
            raise ValueError(
                f"Number of genotypes ({len(genotypes)}) must match "
                f"number of templates ({len(self.templates)})"
            )

        # Combine templates with genotypes
        variants = []
        for template, genotype in zip(self.templates, genotypes):
            variant = template.copy()
            variant["genotype"] = genotype
            variants.append(variant)

        return create_test_variants(variants, output_file=output_file)

    def __len__(self) -> int:
        """Return number of variant templates."""
        return len(self.templates)

    def __repr__(self) -> str:
        """Return string representation."""
        return f"GenotypeGenerator({len(self.templates)} variants)"
