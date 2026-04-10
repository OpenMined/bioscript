from __future__ import annotations

from enum import Enum


class GenotypeEnum(Enum):
    pass


class GenotypeClassifier:
    """
    Base classifier for genotype-based classification.

    Designed for diploid genotype data, returns a DiploidResult
    containing two genotype values (one per chromosome copy).

    Standard interface for bioscript CLI:
    - __call__(matches) -> str: Returns classification result as string
    """

    def __init__(self, *_args, **kwargs):
        """Initialize classifier with optional parameters.

        Common kwargs:
            participant_id: Optional participant/sample ID
            debug: Optional debug flag
            name: Optional classifier name for output files
            filename: Optional input filename
        """
        self.participant_id = kwargs.get("participant_id")
        self.debug = kwargs.get("debug", False)
        self.name = kwargs.get("name", "result")
        self.filename = kwargs.get("filename")

        # Generate base output filename (without extension)
        participant_part = f"_{self.participant_id}" if self.participant_id else ""
        self.output_basename = f"result_{self.name}{participant_part}"

    def classify(self, matches):
        """
        Classify genotypes based on variant matches.

        Args:
            matches: MatchList containing variant matches

        Returns:
            Can return different types depending on classifier needs:
            - DiploidResult: For single-value genotype classification (e.g., APOL1)
            - dict: For multiple summary values (creates namespaced columns)
            - list: For detailed variant data (writes to files, CLI shows count)
            - None: For file-based output only (no CLI output)
            - str/other: Converted to string for output

        Examples:
            # Single value classification
            return DiploidResult(G0, G1)  # → "G0/G1"

            # Multiple summary values
            return {"variant_count": 2, "status": "complete"}  # → BRCA_variant_count, BRCA_status

            # Detailed variant data (for testing and file output)
            write_tsv(f"{self.output_basename}.tsv", variant_data)
            return variant_data  # → CLI shows BRCA_count=2

            # File-based output only (no CLI column)
            write_tsv(f"{self.output_basename}.tsv", variant_data)
            return None  # No CLI output
        """
        raise NotImplementedError("Subclasses must implement classify()")

    def __call__(self, matches):
        """
        Standard callable interface for bioscript CLI.

        Args:
            matches: MatchList containing variant matches

        Returns:
            Classification result as string, dict, list, or None
        """
        result = self.classify(matches)
        if result is None:
            return None
        if isinstance(result, (dict, list)):
            return result
        if hasattr(result, "sorted"):
            result = result.sorted()
        return str(result)


class DiploidResult:
    def __init__(self, genotype1: GenotypeEnum, genotype2: GenotypeEnum):
        self.genotype1 = genotype1
        self.genotype2 = genotype2

    def sorted(self):
        """
        Returns a new DiploidResult with genotypes sorted by enum definition order.

        The genotypes are sorted according to their position in the enum class.
        For example, if the enum is defined as G2, G1, G0, then G0/G2 becomes G2/G0.

        Returns:
            New DiploidResult with sorted genotypes
        """
        # Get the enum class from the first genotype
        if not hasattr(self.genotype1, "__class__") or not issubclass(
            self.genotype1.__class__, Enum
        ):
            # Not an enum, return as-is
            return DiploidResult(self.genotype1, self.genotype2)

        enum_class = self.genotype1.__class__

        # Get all enum members in definition order
        enum_members = list(enum_class)

        # Create a mapping from enum value to its index
        enum_order = {member: i for i, member in enumerate(enum_members)}

        # Sort the genotypes
        sorted_genotypes = sorted(
            [self.genotype1, self.genotype2],
            key=lambda g: enum_order.get(g, float("inf")),
        )

        return DiploidResult(sorted_genotypes[0], sorted_genotypes[1])

    def __str__(self):
        """
        Provides a pretty-printed string representation of the DiploidResult.
        Displays the two genotypes in the format VAL1/VAL2.
        If both genotypes are the same, only display one.
        """
        g1 = self.genotype1.value if hasattr(self.genotype1, "value") else self.genotype1
        g2 = self.genotype2.value if hasattr(self.genotype2, "value") else self.genotype2
        return f"{g1}/{g2}"

    def __repr__(self):
        return self.__str__()
