"""
Position-aware variant matching for bioscript.

This module provides an alternative matching strategy that groups variants
by position and provides clearer NO_CALL semantics.
"""

from collections import defaultdict
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Set

from .types import MatchType, VariantCall, VariantRow


@dataclass
class PositionGroup:
    """Group of variants at the same genomic position."""

    chromosome: str
    position: int
    variants: List[VariantCall] = field(default_factory=list)
    rsids: Set[str] = field(default_factory=set)

    def add_variant(self, variant: VariantCall):
        """Add a variant to this position group."""
        self.variants.append(variant)
        if variant.rsid:
            if isinstance(variant.rsid, list):
                self.rsids.update(variant.rsid)
            else:
                self.rsids.add(str(variant.rsid))

    def describe(self) -> str:
        """Get a human-readable description of variants at this position."""
        descriptions = []
        for v in self.variants:
            ref_str = getattr(v, "ref_label", str(v.ref))
            alt_str = getattr(v, "alt_label", str(v.alt))
            rsid_str = str(v.rsid) if v.rsid else "no_rsid"
            descriptions.append(f"{rsid_str}:{ref_str}→{alt_str}")
        return f"chr{self.chromosome}:{self.position} [{', '.join(descriptions)}]"


@dataclass
class PositionMatch:
    """Result of matching a genotype against all variants at a position."""

    position_group: PositionGroup
    genotype: str
    matched_variant: Optional[VariantCall] = None
    match_type: Optional[MatchType] = None
    all_variants_tested: List[VariantCall] = field(default_factory=list)

    @property
    def is_no_call(self) -> bool:
        """True if this is a NO_CALL at the position."""
        return self.match_type == MatchType.NO_CALL or self.matched_variant is None

    def describe(self) -> str:
        """Get a human-readable description of the match."""
        if self.matched_variant:
            ref_str = getattr(self.matched_variant, "ref_label", str(self.matched_variant.ref))
            alt_str = getattr(self.matched_variant, "alt_label", str(self.matched_variant.alt))
            return (
                f"{self.match_type.name}: {self.matched_variant.rsid or 'unknown'} "
                f"at chr{self.position_group.chromosome}:{self.position_group.position} "
                f"{ref_str}→{alt_str} genotype={self.genotype}"
            )
        else:
            # No match - show all variants that were tested
            variants_desc = []
            for v in self.all_variants_tested:
                ref_str = getattr(v, "ref_label", str(v.ref))
                alt_str = getattr(v, "alt_label", str(v.alt))
                variants_desc.append(f"{ref_str}→{alt_str}")

            return (
                f"NO_CALL at chr{self.position_group.chromosome}:{self.position_group.position} "
                f"genotype={self.genotype} doesn't match any of: {', '.join(variants_desc)}"
            )


class PositionAwareMatchList:
    """
    Alternative to MatchList that groups variants by position before matching.
    This provides clearer NO_CALL semantics.
    """

    def __init__(self, variant_calls: List[VariantCall]):
        """Initialize with variant calls to match against."""
        self.variant_calls = variant_calls
        self.position_groups: Dict[tuple, PositionGroup] = defaultdict(lambda: PositionGroup("", 0))

        # Group variants by position
        for variant in variant_calls:
            if variant.chromosome and variant.position:
                key = (variant.chromosome, variant.position)
                if not self.position_groups[key].chromosome:
                    self.position_groups[key].chromosome = variant.chromosome
                    self.position_groups[key].position = variant.position
                self.position_groups[key].add_variant(variant)

    def match_rows(self, variant_rows: List[VariantRow]) -> List[PositionMatch]:
        """
        Match rows against variant calls, grouped by position.

        Returns a list of PositionMatch objects, one per unique position
        in the input that matches our variant calls.
        """
        results = []
        seen_positions = set()

        for row in variant_rows:
            key = (str(row.chromosome), int(row.position))

            # Skip if we've already processed this position
            if key in seen_positions:
                continue
            seen_positions.add(key)

            # Check if we have variants at this position
            if key not in self.position_groups:
                continue

            position_group = self.position_groups[key]

            # Try to match against each variant at this position
            best_match = None
            best_match_type = None
            best_variant = None
            tested_variants = []

            for variant in position_group.variants:
                match = variant.filter_variant_row(row)
                tested_variants.append(variant)

                if match:
                    # Prioritize actual matches over NO_CALLs
                    if match.match_type in (MatchType.REFERENCE_CALL, MatchType.VARIANT_CALL):
                        best_match = match
                        best_match_type = match.match_type
                        best_variant = variant
                        break  # Found a real match, stop looking
                    elif best_match is None:
                        # Keep NO_CALL as fallback
                        best_match = match
                        best_match_type = match.match_type
                        best_variant = variant

            # Create position match result
            position_match = PositionMatch(
                position_group=position_group,
                genotype=row.genotype,
                matched_variant=best_variant,
                match_type=best_match_type,
                all_variants_tested=tested_variants,
            )
            results.append(position_match)

        return results

    def summarize(self, matches: List[PositionMatch]) -> Dict[str, List[PositionMatch]]:
        """
        Summarize matches by category.

        Returns a dict with keys:
        - 'reference': Positions matching reference
        - 'variant': Positions with variant calls
        - 'no_call': Positions with no matching genotype
        - 'multiple_variants': Positions with multiple possible variants
        """
        summary = {"reference": [], "variant": [], "no_call": [], "multiple_variants": []}

        for match in matches:
            if len(match.position_group.variants) > 1:
                summary["multiple_variants"].append(match)

            if match.match_type == MatchType.REFERENCE_CALL:
                summary["reference"].append(match)
            elif match.match_type == MatchType.VARIANT_CALL:
                summary["variant"].append(match)
            elif match.is_no_call:
                summary["no_call"].append(match)

        return summary


def demo_position_aware_matching():
    """Demo showing how position-aware matching provides clearer results."""

    # Example: Multiple variants at position 17:43045711
    variants = [
        VariantCall(rsid="rs80357336", chromosome="17", position=43045711),  # G→C
        VariantCall(rsid="rs80357336", chromosome="17", position=43045711),  # G→T
        VariantCall(rsid="rs80357629", chromosome="17", position=43045711),  # G→GT
    ]

    # Input row with DD genotype
    rows = [VariantRow(rsid="rs80357629", chromosome="17", position=43045711, genotype="DD")]

    # Match with position awareness
    matcher = PositionAwareMatchList(variants)
    matches = matcher.match_rows(rows)

    # Show results
    for match in matches:
        print(match.describe())

    # Summarize
    summary = matcher.summarize(matches)
    print("\nSummary:")
    print(f"  Positions with multiple variants: {len(summary['multiple_variants'])}")
    print(f"  NO_CALLs: {len(summary['no_call'])}")
    for match in summary["no_call"]:
        print(f"    - {match.position_group.describe()}")
