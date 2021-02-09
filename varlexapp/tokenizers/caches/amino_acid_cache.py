"""A module to cache amino acid codes."""
from typing import Set

# 20 standard amino acids
AMINO_ACID_CODES = {
    'A': 'Ala', 'C': 'Cys', 'D': 'Asp', 'E': 'Glu', 'F': 'Phe', 'G': 'Gly',
    'H': 'His', 'I': 'Ile', 'K': 'Lys', 'L': 'Leu', 'M': 'Met', 'N': 'Asn',
    'P': 'Pro', 'Q': 'Gln', 'R': 'Arg', 'S': 'Ser', 'T': 'Thr', 'V': 'Val',
    'W': 'Trp', 'Y': 'Tyr'
}


class AminoAcidCache:
    """A class to cache amino acid codes."""

    def __init__(self) -> None:
        """Initialize the AminoAcidCache class."""
        self.__amino_acid_codes = self.__load_amino_acid_codes()

    def __contains__(self, item: str) -> bool:
        """Return whether a string in the amino acid cache."""
        return item in self.__amino_acid_codes

    def __load_amino_acid_codes(self) -> Set[str]:
        """Load amino acid cache with amino acid codes.

        :return: A set of valid amino acid codes.
        """
        return (set(AMINO_ACID_CODES.keys()) | set(AMINO_ACID_CODES.values()))
