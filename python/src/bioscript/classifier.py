from __future__ import annotations

from enum import Enum


class GenotypeEnum(Enum):
    pass


class Classifier:
    def classify(self, matches) -> tuple[GenotypeEnum, GenotypeEnum]:
        raise NotImplementedError("Subclasses should implement this method.")


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
            [self.genotype1, self.genotype2], key=lambda g: enum_order.get(g, float("inf"))
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
