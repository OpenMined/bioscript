"""BioScript - A library for analyzing biological scripts and genetic data."""

# BioVault integration
from .biovault import (
    BioVaultProject,
    create_bioscript_project,
    export_notebook_as_project,
    export_workflow,
    load_project,
    new_project,
)
from .classifier import DiploidResult, GenotypeClassifier, GenotypeEnum
from .data import GenotypeGenerator, create_test_variants
from .reader import load_variants_tsv
from .testing import VariantFixture, discover_tests, export_from_notebook, run_tests
from .types import GRCh, MatchType, Nucleotide, VariantCall
from .utils import optional_int, optional_str
from .writer import write_csv, write_tsv

__version__ = "0.1.1"

__all__ = [
    "DiploidResult",
    "GRCh",
    "GenotypeClassifier",
    "GenotypeEnum",
    "GenotypeGenerator",
    "MatchType",
    "Nucleotide",
    "VariantFixture",
    "VariantCall",
    "create_test_variants",
    "discover_tests",
    "export_from_notebook",
    "load_variants_tsv",
    "optional_int",
    "optional_str",
    "run_tests",
    "write_csv",
    "write_tsv",
    # BioVault integration
    "BioVaultProject",
    "create_bioscript_project",
    "export_notebook_as_project",
    "export_workflow",
    "load_project",
    "new_project",
]
