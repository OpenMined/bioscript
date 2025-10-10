"""BioScript - A library for analyzing biological scripts and genetic data."""

from .classifier import DiploidResult, GenotypeClassifier, GenotypeEnum
from .counter import AlleleCount, AlleleCounter
from .data import GenotypeGenerator, create_test_variants
from .reader import load_variants_tsv
from .testing import VariantFixture, discover_tests, export_from_notebook, run_tests
from .types import GRCh, MatchType, Nucleotide, VariantCall

__version__ = "0.1.0"

__all__ = [
    "AlleleCount",
    "AlleleCounter",
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
    "run_tests",
]
