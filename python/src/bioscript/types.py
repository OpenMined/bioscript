from __future__ import annotations

from collections.abc import Iterable
from dataclasses import dataclass, field
from enum import Enum
from typing import Any


class GRCh(str, Enum):
    """Genome Reference Consortium human genome builds."""

    GRCH36 = "GRCh36"
    GRCH37 = "GRCh37"
    GRCH38 = "GRCh38"

    @classmethod
    def parse(cls, value: str | GRCh | None) -> GRCh | None:
        """
        Parse GRCh version from string (case-insensitive) or enum.

        Args:
            value: String like "GRCh38", "grch38", "38" or GRCh enum

        Returns:
            GRCh enum or None if value is None

        Raises:
            ValueError: If string doesn't match a valid GRCh version

        Examples:
            >>> GRCh.parse("GRCh38")
            <GRCh.GRCH38: 'GRCh38'>
            >>> GRCh.parse("grch38")
            <GRCh.GRCH38: 'GRCh38'>
            >>> GRCh.parse("38")
            <GRCh.GRCH38: 'GRCh38'>
            >>> GRCh.parse(GRCh.GRCH38)
            <GRCh.GRCH38: 'GRCh38'>
        """
        if value is None:
            return None

        if isinstance(value, cls):
            return value

        if not isinstance(value, str):
            raise ValueError(f"Expected str or GRCh, got {type(value)}")

        # Normalize to uppercase
        normalized = value.strip().upper()

        # Try exact match first
        for member in cls:
            if member.value.upper() == normalized:
                return member

        # Try with "GRCH" prefix
        if not normalized.startswith("GRCH"):
            normalized = f"GRCH{normalized}"

        for member in cls:
            if member.value.upper() == normalized:
                return member

        valid = ", ".join(m.value for m in cls)
        raise ValueError(f"Invalid GRCh version '{value}'. Valid: {valid} (case-insensitive)")


class Nucleotide(str, Enum):
    A = "A"  # Adenine
    T = "T"  # Thymine (DNA only)
    U = "U"  # Uracil (RNA only)
    C = "C"  # Cytosine
    G = "G"  # Guanine
    MISSING = "-"


class InDel(str, Enum):
    I = "I"  # Insertion
    D = "D"  # Deletion


class MatchType(Enum):
    REFERENCE_CALL = "Reference call"
    VARIANT_CALL = "Variant call"
    NO_CALL = "No call"


class VariantType(str, Enum):
    """Type of genetic variant based on ref/alt sequences."""

    SNV = "snv"  # Single nucleotide variant (ref=1, alt=1)
    INSERTION = "insertion"  # ref=1, alt>1 (e.g., A→AG)
    DELETION = "deletion"  # ref>1, alt=1 (e.g., AG→A)
    INDEL = "indel"  # ref>1, alt>1, same length (e.g., AC→CT)
    MNV = "mnv"  # Multi-nucleotide variant (different lengths)


class AllelesMeta(type):
    @property
    def A(cls):
        return Alleles(Nucleotide.A)

    @property
    def T(cls):
        return Alleles(Nucleotide.T)

    @property
    def C(cls):
        return Alleles(Nucleotide.C)

    @property
    def G(cls):
        return Alleles(Nucleotide.G)

    @property
    def U(cls):
        return Alleles(Nucleotide.U)

    @property
    def I(cls):
        return Alleles(InDel.I)

    @property
    def D(cls):
        return Alleles(InDel.D)

    @property
    def NOT_A(cls):
        return Alleles({n for n in Nucleotide if n != Nucleotide.A and n != Nucleotide.MISSING})

    @property
    def NOT_T(cls):
        return Alleles({n for n in Nucleotide if n != Nucleotide.T and n != Nucleotide.MISSING})

    @property
    def NOT_U(cls):
        return Alleles({n for n in Nucleotide if n != Nucleotide.U and n != Nucleotide.MISSING})

    @property
    def NOT_C(cls):
        return Alleles({n for n in Nucleotide if n != Nucleotide.C and n != Nucleotide.MISSING})

    @property
    def NOT_G(cls):
        return Alleles({n for n in Nucleotide if n != Nucleotide.G and n != Nucleotide.MISSING})


class Alleles(metaclass=AllelesMeta):
    @classmethod
    def from_letter(cls, letter: str):
        """Alternate constructor for a single nucleotide or indel letter."""
        letter = letter.upper()
        mapping = {
            "A": Nucleotide.A,
            "T": Nucleotide.T,
            "C": Nucleotide.C,
            "G": Nucleotide.G,
            "U": Nucleotide.U,
            "I": InDel.I,
            "D": InDel.D,
        }
        if letter not in mapping:
            raise ValueError(f"Invalid allele letter: {letter}")
        return cls(mapping[letter])

    @classmethod
    def from_not_letter(cls, letter: str):
        """Alternate constructor for complement sets, e.g. NOT_A."""
        letter = letter.upper()
        mapping = {
            "A": Nucleotide.A,
            "T": Nucleotide.T,
            "C": Nucleotide.C,
            "G": Nucleotide.G,
            "U": Nucleotide.U,
        }
        if letter not in mapping:
            raise ValueError(f"Invalid NOT_ allele letter: {letter}")
        return cls({n for n in Nucleotide if n != mapping[letter] and n != Nucleotide.MISSING})

    def __init__(
        self,
        nucleotides: Nucleotide | InDel | tuple[Nucleotide | InDel, ...],
    ):
        if isinstance(nucleotides, (Nucleotide, InDel)):
            self.nucleotides = {nucleotides}
        else:
            self.nucleotides = set(nucleotides)

    def __eq__(self, other):
        if isinstance(other, Alleles):
            return self.nucleotides == other.nucleotides
        return False

    def __contains__(self, nucleotide: Nucleotide | InDel) -> bool:
        return nucleotide in self.nucleotides

    def __len__(self) -> int:
        return len(self.nucleotides)

    def __iter__(self):
        return iter(self.nucleotides)

    def add(self, nucleotide: Nucleotide | InDel):
        self.nucleotides.add(nucleotide)

    def remove(self, nucleotide: Nucleotide | InDel):
        self.nucleotides.remove(nucleotide)


@dataclass
class SNP:
    ploidy: str


class DiploidSNP(tuple):
    def __new__(
        cls,
        nucleotide1: Nucleotide | InDel,
        nucleotide2: Nucleotide | InDel,
    ):
        return super().__new__(cls, (nucleotide1, nucleotide2))

    def count(self, nucleotide: Nucleotide | InDel) -> int:
        return super().count(nucleotide)

    def is_homozygous(self) -> bool:
        return self[0] == self[1]

    def is_heterozygous(self) -> bool:
        return self[0] != self[1]


@dataclass
class VariantMatch:
    variant_call: VariantCall
    snp: DiploidSNP
    match_type: MatchType
    source_row: VariantRow | None = None
    context: dict[str, Any] = field(default_factory=dict)

    def __str__(self):
        """
        Provides a pretty-printed string representation of the VariantMatch for easy reading.
        Displays relevant details of the match.
        """
        call = self.variant_call
        genotype = "".join(nuc.value for nuc in self.snp)
        match_type = self.match_type.name

        if call.rsid:
            if isinstance(call.rsid, RSID):
                rsid_repr = "/".join(sorted(call.rsid.aliases))
            else:
                rsid_repr = str(call.rsid)
        elif call.chromosome and call.position is not None:
            rsid_repr = f"chr{call.chromosome}:{call.position}"
        else:
            rsid_repr = "?"

        return f"{match_type}: {rsid_repr} ref={self.snp[0].value} genotype: {genotype}"

    def __repr__(self):
        return self.__str__()

    @property
    def genotype(self) -> DiploidSNP:
        """Return the diploid genotype tuple."""

        return self.snp

    @property
    def genotype_string(self) -> str:
        """Return genotype as joined string (e.g., AG)."""

        return "".join(nuc.value for nuc in self.snp)

    @property
    def genotype_sorted(self) -> str:
        """Return genotype letters sorted alphabetically."""

        return "".join(sorted(nuc.value for nuc in self.snp))

    def count(self, allele: Nucleotide | InDel) -> int:
        """Count specific allele occurrences within the genotype."""

        return self.snp.count(allele)

    @property
    def ref_count(self) -> int:
        """Number of reference alleles present."""

        return sum(1 for allele in self.snp if allele in self.variant_call.ref)

    @property
    def alt_count(self) -> int:
        """Number of alternate alleles present."""

        return sum(1 for allele in self.snp if allele in self.variant_call.alt)

    @property
    def has_variant(self) -> bool:
        """Return True if any alternate allele is present."""

        return self.alt_count > 0

    @property
    def is_homozygous_variant(self) -> bool:
        """Return True if both alleles are alternate alleles."""

        return self.alt_count == 2

    @property
    def is_homozygous_reference(self) -> bool:
        """Return True if both alleles are reference alleles."""

        return self.ref_count == 2

    @property
    def is_heterozygous(self) -> bool:
        """Return True if genotype contains one ref and one alt allele."""

        return self.ref_count == 1 and self.alt_count == 1

    @property
    def has_missing(self) -> bool:
        """Return True if genotype contains any missing allele."""

        return any(allele == Nucleotide.MISSING for allele in self.snp)

    @property
    def raw_line(self) -> str | None:
        """Original TSV line for the matched variant, if available."""

        return self.source_row.raw_line if self.source_row else None

    def as_dict(self) -> dict[str, object]:
        """Return a lightweight summary dictionary for reporting."""

        # Format rsid as a string - use first alias if multiple exist
        rsid_value = None
        if isinstance(self.variant_call.rsid, RSID) and self.variant_call.rsid.aliases:
            rsid_value = sorted(self.variant_call.rsid.aliases)[0]
        elif self.variant_call.rsid:
            rsid_value = str(self.variant_call.rsid)

        return {
            "rsid": rsid_value,
            "chromosome": self.variant_call.chromosome,
            "position": self.variant_call.position,
            "genotype": self.genotype_string,
            "genotype_sorted": self.genotype_sorted,
            "ref_count": self.ref_count,
            "alt_count": self.alt_count,
            "gene": self.variant_call.gene,
            "raw_line": self.raw_line,
        }

    def to_report_row(
        self, participant_id: str | None = None, filename: str | None = None
    ) -> dict[str, object]:
        """Convert VariantMatch to a flat dictionary for TSV reporting.

        Args:
            participant_id: Optional participant/sample ID
            filename: Optional input filename

        Returns:
            Dictionary with all fields formatted for reporting
        """
        vc = self.variant_call

        # Format rsid as string
        rsid = None
        if vc.rsid and hasattr(vc.rsid, "aliases"):
            rsid = sorted(vc.rsid.aliases)[0] if vc.rsid.aliases else None
        elif vc.rsid:
            rsid = str(vc.rsid)

        # Format variant_type
        variant_type = (
            vc.variant_type.value
            if hasattr(vc.variant_type, "value")
            else str(vc.variant_type)
            if vc.variant_type
            else None
        )

        # Format match_type
        match_type = (
            self.match_type.name if hasattr(self.match_type, "name") else str(self.match_type)
        )

        # Format ref and alt
        ref = (
            vc.ref
            if isinstance(vc.ref, str)
            else ",".join(sorted(vc.ref))
            if hasattr(vc.ref, "__iter__")
            else str(vc.ref)
        )
        alt = (
            vc.alt
            if isinstance(vc.alt, str)
            else ",".join(sorted(vc.alt))
            if hasattr(vc.alt, "__iter__")
            else str(vc.alt)
        )

        return {
            "participant_id": participant_id,
            "gene": vc.gene,
            "rsid": rsid,
            "chromosome": vc.chromosome,
            "position": vc.position,
            "genotype": self.genotype_sorted,
            "ref": ref,
            "alt": alt,
            "variant_type": variant_type,
            "match_type": match_type,
            "ref_count": self.ref_count,
            "alt_count": self.alt_count,
            "filename": filename,
        }


@dataclass
class RSID:
    aliases: set

    def __init__(self, *aliases: str | RSID):
        self.aliases = set()
        for alias in aliases:
            if isinstance(alias, RSID):
                self.aliases.update(alias.aliases)
            else:
                self.aliases.add(alias)

    def matches(self, rsid: str | RSID) -> bool:
        if isinstance(rsid, RSID):
            return not self.aliases.isdisjoint(rsid.aliases)
        return rsid in self.aliases


@dataclass
class VariantCall:
    rsid: RSID | str | Iterable | None = None
    ploidy: str = "diploid"
    ref: str | set[str] | Alleles = "-"
    alt: str | set[str] | Alleles = "-"
    chromosome: str | None = None
    position: int | None = None
    assembly: GRCh | str | None = None
    variant_type: VariantType | None = None
    gene: str | None = None

    def __post_init__(self):
        if isinstance(self.rsid, str):
            self.rsid = RSID(self.rsid)
        elif isinstance(self.rsid, Iterable) and not isinstance(self.rsid, RSID):
            self.rsid = RSID(*self.rsid)

        if isinstance(self.assembly, str):
            self.assembly = GRCh.parse(self.assembly)
        elif self.assembly is not None and not isinstance(self.assembly, GRCh):
            raise ValueError(f"assembly must be str, GRCh enum, or None, got {type(self.assembly)}")

        if isinstance(self.chromosome, str):
            self.chromosome = self.chromosome.strip()

        # Auto-detect variant type only for simple string cases
        if self.variant_type is None and isinstance(self.ref, str) and isinstance(self.alt, str):
            ref_len = len(self.ref)
            alt_len = len(self.alt)

            if ref_len == 1 and alt_len == 1:
                self.variant_type = VariantType.SNV
            elif ref_len == 1 and alt_len > 1:
                self.variant_type = VariantType.INSERTION
            elif ref_len > 1 and alt_len == 1:
                self.variant_type = VariantType.DELETION
            elif ref_len == alt_len:
                self.variant_type = VariantType.INDEL
            else:
                self.variant_type = VariantType.MNV

    @staticmethod
    def _normalize_chromosome(value: str | None) -> str | None:
        if value is None:
            return None
        normalized = str(value).strip()
        if not normalized:
            return None

        lower = normalized.lower()
        if lower.startswith("chr"):
            normalized = normalized[3:]

        return normalized

    def _matches_rsid(self, rsid: str | RSID | None) -> bool:
        if self.rsid is None:
            return False
        if rsid is None:
            return False
        if isinstance(rsid, RSID):
            return self.rsid.matches(rsid)

        token = str(rsid).strip()
        if not token:
            return False

        upper_token = token.upper()
        if upper_token in {".", "<NA>", "NA", "N/A"}:
            return False

        return self.rsid.matches(token)

    def _matches_coordinates(
        self,
        chromosome: str | None,
        position: int | None,
        assembly: GRCh | str | None,
    ) -> bool:
        if self.chromosome is None or self.position is None:
            return False
        if chromosome is None or position is None:
            return False

        if isinstance(assembly, str):
            try:
                assembly = GRCh.parse(assembly)
            except ValueError:
                return False

        self_chrom = self._normalize_chromosome(self.chromosome)
        other_chrom = self._normalize_chromosome(chromosome)

        if self_chrom is None or other_chrom is None:
            return False

        same_chromosome = self_chrom.lower() == other_chrom.lower()
        if not same_chromosome:
            return False

        if int(self.position) != int(position):
            return False

        if self.assembly is None or assembly is None:
            return True

        return self.assembly == assembly

    def matches_variant_call(self, other: VariantCall) -> bool:
        if self.rsid and other.rsid and self.rsid.matches(other.rsid):
            return True

        if self._matches_coordinates(other.chromosome, other.position, other.assembly):
            return True

        return other._matches_coordinates(self.chromosome, self.position, self.assembly)

    def _matches_allele(self, nucleotide: str | None, allele: str | set[str] | Alleles) -> bool:
        """
        Check if a nucleotide matches an allele specification.

        Handles three types of allele specifications:
        - str: Single nucleotide ("A") or multi-char sequence ("CTGC")
        - set[str]: Multiple possible alleles ({"C", "G"})
        - Alleles: Set-based allele spec (Alleles.NOT_T, etc.)
        """
        if nucleotide is None:
            return False

        if isinstance(allele, Alleles):
            # Convert string nucleotide to Nucleotide enum for Alleles check
            try:
                nuc_enum = Nucleotide[nucleotide.upper()]
                return nuc_enum in allele
            except (KeyError, AttributeError):
                return False
        elif isinstance(allele, str):
            if len(allele) == 1:
                # Single nucleotide: exact match
                return nucleotide == allele
            else:
                # Multi-char sequence: check if nucleotide appears in sequence
                return nucleotide in allele
        elif isinstance(allele, set):
            # Set of possible alleles
            return nucleotide in allele

        return False

    def _match_snv(self, diploid_snp: DiploidSNP) -> MatchType:
        """Match logic for single nucleotide variants - handles str, set[str], and Alleles."""
        # Convert nucleotides to strings for matching
        nuc0_str = diploid_snp[0].value if isinstance(diploid_snp[0], Nucleotide) else None
        nuc1_str = diploid_snp[1].value if isinstance(diploid_snp[1], Nucleotide) else None

        # Check homozygous reference
        if (
            nuc0_str
            and nuc1_str
            and nuc0_str == nuc1_str
            and self._matches_allele(nuc0_str, self.ref)
        ):
            return MatchType.REFERENCE_CALL

        # Check if either allele matches alt
        if self._matches_allele(nuc0_str, self.alt) or self._matches_allele(nuc1_str, self.alt):
            return MatchType.VARIANT_CALL

        return MatchType.NO_CALL

    def _match_insertion(self, diploid_snp: DiploidSNP) -> MatchType:
        """
        Match logic for insertions (ref=A, alt=AG).
        - I alleles = pathogenic (has insertion)
        - D alleles = reference (no insertion)
        """
        if InDel.I in diploid_snp:
            return MatchType.VARIANT_CALL
        elif InDel.D in diploid_snp:
            return MatchType.REFERENCE_CALL
        else:
            return MatchType.NO_CALL

    def _match_deletion(self, diploid_snp: DiploidSNP) -> MatchType:
        """
        Match logic for deletions (ref=AG, alt=A).
        - D alleles = pathogenic (has deletion)
        - I alleles = reference (no deletion)
        """
        if InDel.D in diploid_snp:
            return MatchType.VARIANT_CALL
        elif InDel.I in diploid_snp:
            return MatchType.REFERENCE_CALL
        else:
            return MatchType.NO_CALL

    def _match_complex(self, diploid_snp: DiploidSNP) -> MatchType:
        """
        Match logic for indels/MNVs (ref=AC, alt=CT or ref=CTGC, alt=TTTA).
        Check if genotype nucleotides appear in ref or alt using _matches_allele().

        Example: rs397509270 ref="AC" alt="CT", genotype=AA
        - A appears in ref (AC) → REFERENCE_CALL

        Example: rs397508883 ref="CTGC", alt="TTTA", genotype=AA
        - A appears only in alt (TTTA) → VARIANT_CALL

        If nucleotide appears in both or neither → NO_CALL
        """
        # Get actual nucleotides (filter out I/D/missing)
        nucs = [a for a in diploid_snp if isinstance(a, Nucleotide) and a != Nucleotide.MISSING]

        if not nucs:
            # Only I/D or missing alleles
            return MatchType.NO_CALL

        # Check if each nucleotide matches ref or alt using the helper
        in_ref_list = [self._matches_allele(nuc.value, self.ref) for nuc in nucs]
        in_alt_list = [self._matches_allele(nuc.value, self.alt) for nuc in nucs]

        all_in_ref = all(in_ref_list)
        all_in_alt = all(in_alt_list)
        any_in_ref = any(in_ref_list)
        any_in_alt = any(in_alt_list)

        if all_in_alt and not any_in_ref:
            # All nucleotides only in alt → pathogenic
            return MatchType.VARIANT_CALL
        elif all_in_ref and not any_in_alt:
            # All nucleotides only in ref → reference
            return MatchType.REFERENCE_CALL
        else:
            # Ambiguous: nucleotides in both, neither, or mixed → NO_CALL
            return MatchType.NO_CALL

    def _determine_match_type(self, diploid_snp: DiploidSNP) -> MatchType:
        """Determine match type based on variant type."""
        if self.variant_type == VariantType.SNV or self.variant_type is None:
            return self._match_snv(diploid_snp)
        elif self.variant_type == VariantType.INSERTION:
            return self._match_insertion(diploid_snp)
        elif self.variant_type == VariantType.DELETION:
            return self._match_deletion(diploid_snp)
        elif self.variant_type in (VariantType.INDEL, VariantType.MNV):
            return self._match_complex(diploid_snp)
        else:
            # Fallback to SNV logic
            return self._match_snv(diploid_snp)

    def filter_variant_row(self, variant_row: VariantRow) -> VariantMatch | None:
        """
        Filters for matching VariantRow based on rsid and/or genomic coordinates.
        """

        if not (
            self._matches_rsid(variant_row.rsid)
            or self._matches_coordinates(
                variant_row.chromosome, variant_row.position, variant_row.assembly
            )
        ):
            return None

        if not (isinstance(variant_row.genotype, str) and len(variant_row.genotype) == 2):
            return None

        def parse_allele(raw: str) -> Nucleotide | InDel:
            """Best-effort parsing of a genotype allele, falling back to missing."""

            normalized = raw.strip()
            if normalized == "-":
                return Nucleotide.MISSING

            try:
                return Nucleotide(normalized)
            except ValueError:
                try:
                    return InDel(normalized)
                except ValueError:
                    print(
                        f"Invalid allele value '{raw}' cast as MISSING '-'",
                        flush=True,
                    )
                    return Nucleotide.MISSING

        allele1 = parse_allele(variant_row.genotype[0])
        allele2 = parse_allele(variant_row.genotype[1])

        diploid_snp = DiploidSNP(allele1, allele2)

        # Use type-specific matching logic
        match_type = self._determine_match_type(diploid_snp)
        return VariantMatch(
            variant_call=self,
            snp=diploid_snp,
            match_type=match_type,
            source_row=variant_row,
        )


@dataclass
class VariantRow:
    rsid: str
    chromosome: str  # str to allow "X","Y","MT"
    position: int
    genotype: str  # Keep as str until needed
    ploidy: str = "diploid"
    assembly: GRCh | str | None = None  # Genome reference build (e.g., GRCh37, GRCh38)
    gs: float | None = None
    baf: float | None = None
    lrr: float | None = None
    raw_line: str | None = None

    def __post_init__(self):
        # Parse assembly to GRCh enum if it's a string
        if isinstance(self.assembly, str):
            self.assembly = GRCh.parse(self.assembly)
        elif self.assembly is not None and not isinstance(self.assembly, GRCh):
            raise ValueError(f"assembly must be str, GRCh enum, or None, got {type(self.assembly)}")


@dataclass
class MatchList:
    variant_calls: Iterable[VariantCall]
    enable_position_clustering: bool = True
    reference_matches: list = field(default_factory=list)
    variant_matches: list = field(default_factory=list)
    no_call_matches: list = field(default_factory=list)
    all_matches: list = field(default_factory=list)
    match_lookup: dict[str, VariantMatch] = field(init=False, default_factory=dict, repr=False)
    _variant_call_index: dict[str, list[VariantCall]] = field(
        init=False, default_factory=dict, repr=False
    )
    _variant_call_order: dict[int, int] = field(init=False, default_factory=dict, repr=False)
    _position_clusters: dict[tuple[str, int], list[dict[str, Any]]] = field(
        init=False,
        default_factory=dict,
        repr=False,
    )

    def __post_init__(self):
        self.variant_calls = list(self.variant_calls)
        self.match_lookup = {}
        self._variant_call_index = {}
        self._variant_call_order = {id(call): idx for idx, call in enumerate(self.variant_calls)}
        self._position_clusters = {}

        for call in self.variant_calls:
            for key in self._keys_for_variant_call(call):
                normalized = key.strip().lower()
                if not normalized:
                    continue
                self._variant_call_index.setdefault(normalized, []).append(call)

    def match_rows(
        self,
        variant_rows: Iterable[VariantRow],
        *,
        enable_multi_variant: bool | None = None,
    ) -> MatchList:
        """Match rows against configured VariantCall objects with position-aware grouping."""

        self.reference_matches = []
        self.variant_matches = []
        self.no_call_matches = []
        self.all_matches = []
        self.match_lookup = {}
        self._position_clusters = {}

        if enable_multi_variant is None:
            clustering_enabled = self.enable_position_clustering
        else:
            clustering_enabled = bool(enable_multi_variant)
            self.enable_position_clustering = clustering_enabled

        # Group variants by position for better reporting
        position_matches = {}  # (chrom, pos) -> list of (variant_call, match)

        for variant_row in variant_rows:
            position_key = (str(variant_row.chromosome), int(variant_row.position))

            # Get all candidate variant calls for this row
            candidate_calls = self._candidate_calls_for_row(variant_row)

            # Try to match against each candidate
            matches_at_position = []
            for variant_call in candidate_calls:
                match = variant_call.filter_variant_row(variant_row)
                if match is not None:
                    matches_at_position.append((variant_call, match))

            if matches_at_position:
                # Store all matches at this position for later prioritization
                if position_key not in position_matches:
                    position_matches[position_key] = []
                position_matches[position_key].extend(matches_at_position)

        # Process matches by position
        for position_key, matches in position_matches.items():
            chrom_key = str(position_key[0])
            pos_key = int(position_key[1])

            # Check if we have any real matches (not NO_CALLs) at this position
            real_matches = [
                (vc, m)
                for vc, m in matches
                if m.match_type in (MatchType.REFERENCE_CALL, MatchType.VARIANT_CALL)
            ]
            no_call_matches = [(vc, m) for vc, m in matches if m.match_type == MatchType.NO_CALL]

            def compute_priority_components(
                item: tuple[VariantCall, VariantMatch],
            ) -> dict[str, bool]:
                variant_call, match = item
                has_rsid = variant_call.rsid is not None
                has_coordinates = (
                    variant_call.chromosome is not None and variant_call.position is not None
                )

                exact_coordinate = False
                exact_rsid = False

                if match.source_row is not None:
                    row = match.source_row
                    if row.chromosome and row.position is not None and has_coordinates:
                        row_chrom = VariantCall._normalize_chromosome(row.chromosome)
                        call_chrom = VariantCall._normalize_chromosome(variant_call.chromosome)
                        if (
                            row_chrom
                            and call_chrom
                            and row_chrom.lower() == call_chrom.lower()
                            and int(row.position) == int(variant_call.position)
                        ):
                            exact_coordinate = True

                    if has_rsid and row.rsid:
                        rsid_token = row.rsid.strip()
                        if (
                            rsid_token
                            and rsid_token.upper() not in {".", "<NA>", "NA", "N/A"}
                            and variant_call._matches_rsid(rsid_token)
                        ):
                            exact_rsid = True

                return {
                    "has_rsid": has_rsid,
                    "has_coordinates": has_coordinates,
                    "exact_rsid": exact_rsid,
                    "exact_coordinate": exact_coordinate,
                }

            if real_matches:
                # Sort matches to prioritize exact rsID + coordinate matches

                def match_priority(item: tuple[VariantCall, VariantMatch]) -> tuple[int, int, int]:
                    components = compute_priority_components(item)
                    exact_rsid = int(components["exact_rsid"])
                    exact_coordinate = int(components["exact_coordinate"])
                    has_rsid = int(components["has_rsid"])
                    has_coordinates = int(components["has_coordinates"])

                    # Priority tuple: higher values should sort first, so negative for descending order
                    return (
                        -(exact_rsid and exact_coordinate),
                        -exact_rsid,
                        -exact_coordinate,
                        -has_rsid,
                        -has_coordinates,
                    )

                real_matches.sort(key=match_priority)

                # We have real matches - only report the best real match
                # Priority: VARIANT_CALL > REFERENCE_CALL
                best_pair = real_matches[0]
                best_match = best_pair[1]
                if best_match.match_type != MatchType.VARIANT_CALL:
                    for pair in real_matches:
                        if pair[1].match_type == MatchType.VARIANT_CALL:
                            best_pair = pair
                            best_match = pair[1]
                            break

                if best_match:
                    if clustering_enabled:
                        # Reorder cluster entries to put higher-priority matches first
                        ordered_matches = real_matches + [
                            item for item in matches if item not in real_matches
                        ]
                        priority_components = {
                            id(match): compute_priority_components((variant_call, match))
                            for variant_call, match in ordered_matches
                        }
                        cluster_entries = self._build_position_cluster(
                            ordered_matches, best_match, priority_components
                        )
                        self._position_clusters[(chrom_key, pos_key)] = cluster_entries

                        alternative_count = sum(
                            1 for entry in cluster_entries if not entry["is_selected"]
                        )
                        same_rsid_alternatives = sum(
                            1
                            for entry in cluster_entries
                            if entry["shares_rsid_with_selected"] and not entry["is_selected"]
                        )

                        best_match.context = {
                            "position_cluster": cluster_entries,
                            "ambiguity": {
                                "position_key": {
                                    "chromosome": chrom_key,
                                    "position": pos_key,
                                },
                                "has_alternatives": alternative_count > 0,
                                "alternative_count": alternative_count,
                                "same_rsid_alternative_count": same_rsid_alternatives,
                            },
                            "priority": priority_components.get(id(best_match)),
                        }
                    else:
                        best_match.context = {
                            "priority": compute_priority_components(best_pair),
                        }

                    self.all_matches.append(best_match)
                    if best_match.match_type == MatchType.REFERENCE_CALL:
                        self.reference_matches.append(best_match)
                    elif best_match.match_type == MatchType.VARIANT_CALL:
                        self.variant_matches.append(best_match)
                    self._register_match(best_match)

            elif no_call_matches:
                # Only NO_CALLs at this position - report all of them to show
                # that multiple variants were tested but none matched
                priority_components_all = {
                    id(match): compute_priority_components((variant_call, match))
                    for variant_call, match in matches
                }
                if clustering_enabled:
                    self._position_clusters[(chrom_key, pos_key)] = self._build_position_cluster(
                        matches,
                        None,
                        priority_components_all,
                    )
                for _variant_call, match in no_call_matches:
                    if clustering_enabled:
                        cluster_entries = self._build_position_cluster(
                            matches, match, priority_components_all
                        )
                        alternative_count = sum(
                            1 for entry in cluster_entries if not entry["is_selected"]
                        )
                        same_rsid_alternatives = sum(
                            1
                            for entry in cluster_entries
                            if entry["shares_rsid_with_selected"] and not entry["is_selected"]
                        )

                        match.context = {
                            "position_cluster": cluster_entries,
                            "ambiguity": {
                                "position_key": {
                                    "chromosome": chrom_key,
                                    "position": pos_key,
                                },
                                "has_alternatives": alternative_count > 0,
                                "alternative_count": alternative_count,
                                "same_rsid_alternative_count": same_rsid_alternatives,
                            },
                            "priority": priority_components_all.get(id(match)),
                        }
                    else:
                        match.context = {
                            "priority": priority_components_all.get(id(match)),
                        }

                    self.all_matches.append(match)
                    self.no_call_matches.append(match)
                    self._register_match(match)

        return self

    def _rsid_alias_set(self, rsid: RSID | str | None) -> set[str]:
        if isinstance(rsid, RSID):
            return {str(alias) for alias in rsid.aliases}
        if rsid:
            return {str(rsid)}
        return set()

    def _format_allele_for_context(self, allele: str | set[str] | Alleles) -> str:
        if isinstance(allele, str):
            return allele
        if isinstance(allele, set):
            return "/".join(sorted(str(a) for a in allele))
        if isinstance(allele, Alleles):
            values: list[str] = []
            for token in allele.nucleotides:
                if hasattr(token, "value"):
                    values.append(str(token.value))
                else:
                    values.append(str(token))
            return "/".join(sorted(values))
        return str(allele)

    def _format_genotype_for_context(self, snp: DiploidSNP) -> str:
        symbols: list[str] = []
        for allele in snp:
            if hasattr(allele, "value"):
                symbols.append(str(allele.value))
            else:
                symbols.append(str(allele))
        return "".join(symbols)

    def _format_variant_type(self, variant_call: VariantCall) -> str | None:
        vtype = getattr(variant_call, "variant_type", None)
        if vtype is None:
            return None
        if hasattr(vtype, "value"):
            return str(vtype.value)
        return str(vtype)

    def _build_position_cluster(
        self,
        matches: list[tuple[VariantCall, VariantMatch]],
        selected_match: VariantMatch | None,
        priority_components: dict[int, dict[str, bool]] | None = None,
    ) -> list[dict[str, Any]]:
        selected_aliases = (
            self._rsid_alias_set(selected_match.variant_call.rsid)
            if selected_match and selected_match.variant_call
            else set()
        )

        cluster: list[dict[str, Any]] = []
        for variant_call, match in matches:
            aliases = self._rsid_alias_set(variant_call.rsid)
            alias_list = sorted(aliases)

            chrom = variant_call.chromosome
            pos = variant_call.position
            if chrom is None and match.source_row is not None:
                chrom = match.source_row.chromosome
            if pos is None and match.source_row is not None:
                pos = match.source_row.position

            entry = {
                "is_selected": match is selected_match,
                "match_type": match.match_type.name,
                "rsids": alias_list,
                "rsid_display": "/".join(alias_list) if alias_list else None,
                "chromosome": str(chrom) if chrom is not None else None,
                "position": int(pos) if pos is not None else None,
                "ref": self._format_allele_for_context(variant_call.ref),
                "alt": self._format_allele_for_context(variant_call.alt),
                "genotype": self._format_genotype_for_context(match.snp),
                "variant_type": self._format_variant_type(variant_call),
                "raw_line": match.source_row.raw_line if match.source_row else None,
                "shares_rsid_with_selected": bool(selected_aliases & aliases)
                if selected_aliases
                else False,
                "priority": priority_components.get(id(match)) if priority_components else None,
            }
            cluster.append(entry)

        return cluster

    def __iter__(self):
        return iter(self.all_matches)

    def __getitem__(self, index):
        return self.all_matches[index]

    def iter_reference(self):
        return iter(self.reference_matches)

    def iter_variant(self):
        return iter(self.variant_matches)

    def iter_no_call(self):
        return iter(self.no_call_matches)

    def categorize_matches(self) -> tuple[list, list, list]:
        """Return matches categorized as (reference_calls, variant_calls, no_calls)."""
        return (self.reference_matches, self.variant_matches, self.no_call_matches)

    def categorize_report_rows(
        self, participant_id: str | None = None, filename: str | None = None
    ) -> tuple[list[dict], list[dict], list[dict]]:
        """Return matches categorized and converted to report rows ready for TSV output.

        Args:
            participant_id: Optional participant/sample ID to include in each row
            filename: Optional filename to include in each row

        Returns:
            Tuple of (reference_rows, variant_rows, no_call_rows) where each is a list of dicts
        """
        ref_rows = [m.to_report_row(participant_id, filename) for m in self.reference_matches]
        var_rows = [m.to_report_row(participant_id, filename) for m in self.variant_matches]
        no_rows = [m.to_report_row(participant_id, filename) for m in self.no_call_matches]
        return (ref_rows, var_rows, no_rows)

    def get_position_summary(self) -> dict:
        """Get a summary of matches grouped by position."""
        from collections import defaultdict

        position_groups = defaultdict(list)
        for match in self.all_matches:
            call = match.variant_call
            if call.chromosome and call.position is not None:
                chrom_key = str(call.chromosome)
                pos_key = int(call.position)
                key = (chrom_key, pos_key)
                position_groups[key].append(match)

        # Create summary
        summary = {
            "total_positions": len(position_groups),
            "positions_with_multiple_variants": 0,
            "details": [],
        }

        for (chrom, pos), matches in position_groups.items():
            if len(matches) > 1:
                summary["positions_with_multiple_variants"] += 1

            position_info = {
                "chromosome": chrom,
                "position": pos,
                "match_count": len(matches),
                "matches": [],
            }

            cluster = self._position_clusters.get((str(chrom), int(pos)))
            if cluster:
                position_info["cluster_size"] = len(cluster)
                position_info["alternative_count"] = sum(
                    1 for entry in cluster if not entry.get("is_selected")
                )
                position_info["same_rsid_alternative_count"] = sum(
                    1
                    for entry in cluster
                    if entry.get("shares_rsid_with_selected") and not entry.get("is_selected")
                )

            for match in matches:
                call = match.variant_call
                match_info = {
                    "rsid": str(call.rsid) if call.rsid else None,
                    "ref": getattr(call, "ref_label", str(call.ref)),
                    "alt": getattr(call, "alt_label", str(call.alt)),
                    "match_type": match.match_type.name,
                    "genotype": "".join(
                        n.value if hasattr(n, "value") else str(n) for n in match.snp
                    ),
                    "variant_type": self._format_variant_type(call),
                }
                position_info["matches"].append(match_info)

        summary["details"].append(position_info)

        return summary

    def __str__(self):
        """
        Provides a pretty-printed string representation of the MatchList for easy reading.
        Each match is displayed on its own line with relevant details.
        """
        return "\n".join(str(match) for match in self.all_matches)

    def __repr__(self):
        return self.__str__()

    def _register_match(self, match: VariantMatch) -> None:
        """Index match by rsID aliases and coordinates for fast lookup."""

        keys: set[str] = set()

        rsid = match.variant_call.rsid
        if isinstance(rsid, RSID):
            keys.update(rsid.aliases)
        elif rsid:
            keys.add(str(rsid))

        chrom = match.variant_call.chromosome
        pos = match.variant_call.position
        if chrom and pos is not None:
            normalized_chrom = match.variant_call._normalize_chromosome(chrom)
            if normalized_chrom:
                keys.add(f"chr{normalized_chrom}:{pos}")
                keys.add(f"{normalized_chrom}:{pos}")

        for key in keys:
            normalized = str(key).strip().lower()
            if normalized and normalized not in self.match_lookup:
                self.match_lookup[normalized] = match

    def find_match(self, variant_call: VariantCall):
        """Return the first match corresponding to the VariantCall, if any."""

        for key in self._keys_for_variant_call(variant_call):
            normalized = key.strip().lower()
            if not normalized:
                continue
            match = self.match_lookup.get(normalized)
            if match is not None:
                return match
        return None

    def get(self, identifier, default=None):
        """Lookup a match by VariantCall, rsID alias, or coordinate string."""

        if isinstance(identifier, VariantCall):
            match = self.find_match(identifier)
            return match if match is not None else default

        if isinstance(identifier, RSID):
            for alias in sorted(identifier.aliases):
                match = self.get(alias)
                if match is not None:
                    return match
            return default

        if identifier is None:
            return default

        key = str(identifier).strip().lower()
        if not key:
            return default

        return self.match_lookup.get(key, default)

    def __getattr__(self, item: str):
        if item.startswith("__"):
            raise AttributeError(f"MatchList has no attribute '{item}'")

        match = self.get(item)
        if match is not None:
            return match
        raise AttributeError(f"MatchList has no attribute '{item}'")

    def _keys_for_variant_call(self, call: VariantCall) -> set[str]:
        keys: set[str] = set()

        rsid = call.rsid
        if isinstance(rsid, RSID):
            keys.update(str(alias) for alias in rsid.aliases)
        elif rsid:
            keys.add(str(rsid))

        chrom = call.chromosome
        pos = call.position
        if chrom and pos is not None:
            normalized = VariantCall._normalize_chromosome(chrom)
            if normalized:
                keys.add(f"{normalized}:{pos}")
                keys.add(f"chr{normalized}:{pos}")

        return keys

    def _keys_for_variant_row(self, variant_row: VariantRow) -> set[str]:
        keys: set[str] = set()

        rsid = (variant_row.rsid or "").strip()
        if rsid:
            upper = rsid.upper()
            if upper not in {".", "<NA>", "NA", "N/A"}:
                keys.add(rsid.lower())

        chrom = VariantCall._normalize_chromosome(variant_row.chromosome)
        if chrom and variant_row.position is not None:
            pos = int(variant_row.position)
            keys.add(f"{chrom}:{pos}".lower())
            keys.add(f"chr{chrom}:{pos}".lower())

        return keys

    def _candidate_calls_for_row(self, variant_row: VariantRow) -> list[VariantCall]:
        keys = self._keys_for_variant_row(variant_row)
        if not keys:
            # No keys means we can't index this row - return empty list
            # Only check if variant calls have no indexable keys either
            return [] if self._variant_call_index else self.variant_calls

        candidate_ids: set[int] = set()
        for key in keys:
            candidate_ids.update(id(call) for call in self._variant_call_index.get(key, []))

        if not candidate_ids:
            # No matches found in the index - return empty list
            return []

        ordered_calls: list[VariantCall] = []
        remaining = set(candidate_ids)
        for call in self.variant_calls:
            ident = id(call)
            if ident in remaining:
                ordered_calls.append(call)
                remaining.remove(ident)
                if not remaining:
                    break

        return ordered_calls
