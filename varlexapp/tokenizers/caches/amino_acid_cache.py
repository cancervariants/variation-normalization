"""A module to cache amino acid codes."""
from typing import Set

# 20 standard amino acids
THREE_LETTER_CODE = {'Ala', 'Cys', 'Asp', 'Glu', 'Phe', 'Gly', 'His', 'Ile',
                     'Lys', 'Leu', 'Met', 'Asn', 'Pro', 'Gln', 'Arg', 'Ser',
                     'Thr', 'Val', 'Trp', 'Tyr'}

ONE_LETTER_CODE = set('ACDEFGHIKLMNPQRSTVWY')


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
        return (THREE_LETTER_CODE | ONE_LETTER_CODE)
