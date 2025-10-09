"""Basic tests for bioscript package."""

import bioscript


def test_version():
    """Test that version is defined."""
    assert hasattr(bioscript, "__version__")
    assert isinstance(bioscript.__version__, str)


def test_imports():
    """Test that main classes can be imported."""
    assert hasattr(bioscript, "Classifier")
    assert hasattr(bioscript, "DiploidResult")
    assert hasattr(bioscript, "GenotypeEnum")
    assert hasattr(bioscript, "VariantCall")
    assert hasattr(bioscript, "load_variants_tsv")
