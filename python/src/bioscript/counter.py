"""Allele counting utilities for genotype classification."""

from __future__ import annotations

from .types import VariantCall


class AlleleCount:
    """
    Result of counting alleles at a genomic position.

    Provides convenient access to allele counts and genotype information.

    Attributes:
        match: The VariantMatch object (None if no match found)
        variant_call: The VariantCall used for matching
        ref_count: Number of reference alleles (0-2)
        alt_count: Number of alternate/variant alleles (0-2)
        genotype: The DiploidSNP genotype (None if no match)
    """

    def __init__(self, match=None, variant_call=None):
        """Initialize allele count from a match."""
        self.match = match
        self.variant_call = variant_call
        self.ref_count = 0
        self.alt_count = 0
        self.genotype = None

        if match and variant_call:
            self.genotype = match.snp

            # Count reference and alternate alleles
            for allele in match.snp:
                if allele in variant_call.ref:
                    self.ref_count += 1
                if allele in variant_call.alt:
                    self.alt_count += 1

    def count(self, allele):
        """
        Count occurrences of a specific allele (0-2).

        Args:
            allele: Nucleotide or InDel to count

        Returns:
            Number of times the allele appears (0, 1, or 2)

        Example:
            >>> from bioscript.types import Nucleotide
            >>> result.count(Nucleotide.A)  # How many A alleles?
            1
        """
        if self.genotype:
            return self.genotype.count(allele)
        return 0

    def has_variant(self) -> bool:
        """Check if any variant (alt) alleles are present."""
        return self.alt_count > 0

    def is_homozygous_variant(self) -> bool:
        """Check if homozygous for variant allele (e.g., GG when ref is A)."""
        return self.alt_count == 2

    def is_heterozygous(self) -> bool:
        """Check if heterozygous (one ref, one alt)."""
        return self.ref_count == 1 and self.alt_count == 1

    def is_homozygous_reference(self) -> bool:
        """Check if homozygous for reference allele."""
        return self.ref_count == 2

    def __repr__(self) -> str:
        """String representation."""
        if self.genotype:
            return (
                f"AlleleCount(genotype={''.join(a.value for a in self.genotype)}, "
                f"ref={self.ref_count}, alt={self.alt_count})"
            )
        return "AlleleCount(no_match)"


class AlleleCounter:
    """
    Count alleles at a genomic position across matches.

    This helper extracts and counts alleles for a specific variant call,
    making it easy to classify genotypes without dealing with match iteration.

    Args:
        variant_call: The VariantCall to count alleles for

    Example:
        >>> from bioscript import AlleleCounter
        >>> from bioscript.types import VariantCall, Alleles
        >>> call = VariantCall(rsid="rs123", ref=Alleles.A, alt=Alleles.G)
        >>> counter = AlleleCounter(call)
        >>> result = counter.count(matches)
        >>> result.alt_count  # How many G alleles?
        1
        >>> result.is_heterozygous()
        True
    """

    def __init__(self, variant_call: VariantCall):
        """Initialize counter for a specific variant call."""
        self.variant_call = variant_call

    def count(self, matches) -> AlleleCount:
        """
        Count alleles in the match list.

        Args:
            matches: MatchList object containing variant matches

        Returns:
            AlleleCount with genotype and count information
        """
        # Find the first match for this rsID (any match type)
        for match in matches.all_matches:
            if match.variant_call.rsid.matches(self.variant_call.rsid):
                return AlleleCount(match=match, variant_call=self.variant_call)

        # No match found
        return AlleleCount()

    def __repr__(self) -> str:
        """String representation."""
        rsid = (
            self.variant_call.rsid.aliases
            if hasattr(self.variant_call.rsid, "aliases")
            else self.variant_call.rsid
        )
        return f"AlleleCounter(rsid={rsid})"
