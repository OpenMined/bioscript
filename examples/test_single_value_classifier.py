"""Test classifier that returns a single value instead of a dict."""

from bioscript.classifier import GenotypeClassifier
from bioscript.types import Alleles, VariantCall

rs123 = VariantCall(rsid="rs123", ref=Alleles.A, alt=Alleles.T)

class SimpleClassifier(GenotypeClassifier):
    def classify(self, matches):
        match = matches.get(rs123)
        if not match:
            return "No match"
        if match.genotype_sorted == "AA":
            return "Type A"
        elif match.genotype_sorted == "AT":
            return "Type B"
        elif match.genotype_sorted == "TT":
            return "Type C"
        else:
            return "Unknown"

__bioscript__ = {
    "variant_calls": [rs123],
    "classifier": SimpleClassifier(),
    "name": "SIMPLE",
}