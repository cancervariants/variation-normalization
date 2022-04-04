"""A module to cache amino acid codes."""
from typing import Set
import csv

from variation import AMINO_ACID_PATH


class AminoAcidCache:
    """A class to cache amino acid codes."""

    def __init__(self, amino_acids_file_path: str = AMINO_ACID_PATH) -> None:
        """Initialize the AminoAcidCache class."""
        self._amino_acids_file = amino_acids_file_path
        self.amino_acid_code_conversion = dict()
        self.__amino_acid_codes = self.__load_amino_acid_codes()

    def __contains__(self, item: str) -> bool:
        """Return whether a string in the amino acid cache."""
        return item.upper() in self.__amino_acid_codes

    def __load_amino_acid_codes(self) -> Set[str]:
        """Load amino acid cache with amino acid codes.

        :return: A set of valid amino acid codes.
        """
        with open(self._amino_acids_file, "r") as f:
            data = list(csv.reader(f))
            for row in data:
                self.amino_acid_code_conversion[row[2]] = row[1]

        return {item.upper() for sublist in data for item in sublist}

    def convert_three_to_one(self, three_letter_amino_acid: str) -> str:
        """Convert a 3 letter amino acid code to a 1 letter amino acid code.

        :param str three_letter_amino_acid: Amino Acid Code to convert
        :return: A str of the one letter protein code
        """
        if three_letter_amino_acid.upper() == "TER":
            return "*"
        for one_letter, three_letter in self.amino_acid_code_conversion.items():  # noqa: E501
            if three_letter.upper() == three_letter_amino_acid.upper():
                return one_letter.upper()
