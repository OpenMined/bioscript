"""
BioScript Data Utilities

Fetch and manage test genomic datasets for bioinformatics workflows.
"""

import os
import zipfile
from pathlib import Path
from typing import Optional
from urllib.request import urlretrieve


# Sample data registry
SAMPLES = {
    "23andme_v4": {
        "url": "https://github.com/OpenMined/biovault-data/raw/main/snp/23andme_genome_v4_Full.zip",
        "description": "23andMe Genome v4 Full - Complete SNP genotyping dataset",
        "filename": "23andme_genome_v4_Full.zip",
    }
}


def fetch_sample(
    sample_name: str,
    output_dir: Optional[str] = None,
    force: bool = False
) -> Path:
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
        raise ValueError(
            f"Unknown sample '{sample_name}'. Available samples: {available}"
        )

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

        with zipfile.ZipFile(zip_path, 'r') as zip_ref:
            zip_ref.extractall(extract_dir)

        # List extracted files
        extracted_files = list(extract_dir.rglob('*'))
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
    all_files = [f for f in extract_dir.rglob('*') if f.is_file()]
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
        raise ValueError(
            f"Unknown sample '{sample_name}'. Available samples: {available}"
        )

    return SAMPLES[sample_name].copy()
