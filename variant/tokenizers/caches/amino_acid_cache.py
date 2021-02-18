"""A module to cache amino acid codes."""
from typing import Set
from variant import AMINO_ACID_PATH
import csv


class AminoAcidCache:
    """A class to cache amino acid codes."""

    def __init__(self,
                 amino_acids_file_path=AMINO_ACID_PATH) -> None:
        """Initialize the AminoAcidCache class."""
        self._amino_acids_file = amino_acids_file_path
        self._amino_acid_code_conversion = dict()
        self.__amino_acid_codes = self.__load_amino_acid_codes()

    def __contains__(self, item: str) -> bool:
        """Return whether a string in the amino acid cache."""
        return item in self.__amino_acid_codes

    def __load_amino_acid_codes(self) -> Set[str]:
        """Load amino acid cache with amino acid codes.

        :return: A set of valid amino acid codes.
        """
        with open(self._amino_acids_file, 'r') as f:
            data = list(csv.reader(f))
            for row in data:
                self._amino_acid_code_conversion[row[2]] = row[1]

        return ({item for sublist in data for item in sublist})
