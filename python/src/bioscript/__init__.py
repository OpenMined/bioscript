"""BioScript - A library for analyzing biological scripts and genetic data."""

from .classifier import Classifier, DiploidResult, GenotypeEnum
from .counter import AlleleCount, AlleleCounter
from .data import GenotypeGenerator, create_test_variants
from .reader import load_variants_tsv
from .types import MatchType, Nucleotide, VariantCall

__version__ = "0.1.0"

__all__ = [
    "AlleleCount",
    "AlleleCounter",
    "Classifier",
    "DiploidResult",
    "GenotypeEnum",
    "GenotypeGenerator",
    "MatchType",
    "Nucleotide",
    "VariantCall",
    "create_test_variants",
    "load_variants_tsv",
]
