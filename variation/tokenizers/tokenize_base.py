"""Module for commonly used tokenization methods."""
from typing import Tuple, Optional, Union, List
import re

from variation.tokenizers.caches import NucleotideCache, AminoAcidCache


class TokenizeBase:
    """Class for Tokenize methods."""

    def __init__(self, amino_acid_cache: AminoAcidCache,
                 nucleotide_cache: NucleotideCache) -> None:
        """Initialize Token Base class.

        :param AminoAcidCache amino_acid_cache: Valid amino acid codes
        :param NucleotideCache nucleotide_cache: Valid nucleotides
        """
        self.nucleotide_cache = nucleotide_cache
        self.amino_acid_cache = amino_acid_cache
        self.splitter_char_digit = re.compile("([a-zA-Z]+)([0-9]+)")

    def get_amino_acid_and_pos(
        self, part: List, used_one_letter: bool
    ) -> Optional[Tuple[str, int, bool]]:
        """Return amino acid and position.

        :param List part: Tokenized input string
        :param bool used_one_letter: `True` if used 1 letter AA code.
            `False` if used 3 letter AA code.
        :return: Three letter AA code, position, and whether or not
            one letter AA code was used
        """
        try:
            char_and_digits = self.splitter_char_digit.match(part).groups()
        except AttributeError:
            return None

        if len(char_and_digits) != 2 or len(part) != \
                (len(char_and_digits[0]) + len(char_and_digits[1])):
            return None

        aa = char_and_digits[0]
        if len(aa) == 1:
            if not used_one_letter:
                used_one_letter = True
            tmp_aa = \
                self.amino_acid_cache.amino_acid_code_conversion[aa.upper()]
        else:
            tmp_aa = aa
        pos = char_and_digits[1]

        if not self.amino_acid_cache.__contains__(tmp_aa) \
                or not pos.isdigit():
            return None
        return aa.upper(), pos, used_one_letter

    def get_protein_inserted_sequence(self, parts: List,
                                      used_one_letter: bool) -> Optional[str]:
        """Return inserted sequence for protein reference sequence.

        :param List parts: Tokenized input string
        :param bool used_one_letter: `True` if used 1 letter AA code.
            `False` if used 3 letter AA code.
        :return: Inserted sequence
        """
        # Check inserted sequences
        inserted_sequence = ""
        if used_one_letter:
            for i in range(len(parts[1])):
                aa = parts[1][i:i + 1]
                if len(aa) != 1:
                    return None
                try:
                    self.amino_acid_cache.amino_acid_code_conversion[aa.upper()]  # noqa: E501
                except KeyError:
                    return None
                else:
                    inserted_sequence += aa.upper()
        else:
            for i in range(0, len(parts[1]), 3):
                aa = parts[1][i:i + 3]
                if len(aa) != 3 or not self.amino_acid_cache.__contains__(aa):
                    if aa != "ter":
                        return None
                inserted_sequence += \
                    self.amino_acid_cache.convert_three_to_one(aa)

        if inserted_sequence == "":
            return None
        return inserted_sequence

    def get_aa_pos_range(self,
                         parts: List) -> Optional[Tuple[str, str, str, int, bool]]:
        """Get amino acid(s) and positions(s) for protein reference sequence.

        :param List parts: Tokenized input string
        :return: Beginning AA, End AA,  Beginning position, End position,
            Whether or not one letter code was used
        """
        aa_start = None
        aa_end = None
        pos_start = None
        pos_end = None
        used_one_letter = False

        if "_" in parts[0] and parts[0].count("_") == 1:
            aa_pos_range = parts[0].split("_")
            if len(aa_pos_range) != 2 or \
                    not aa_pos_range[0] or not aa_pos_range[1]:
                return None

            start_aa_pos = \
                self.get_amino_acid_and_pos(
                    aa_pos_range[0], used_one_letter
                )

            if start_aa_pos:
                used_one_letter = start_aa_pos[2]

            end_aa_pos = \
                self.get_amino_acid_and_pos(
                    aa_pos_range[1], used_one_letter
                )

            if start_aa_pos and end_aa_pos:
                aa_start = start_aa_pos[0]
                pos_start = start_aa_pos[1]
                aa_end = end_aa_pos[0]
                pos_end = end_aa_pos[1]
                used_one_letter = end_aa_pos[2]

        else:
            aa_and_pos = \
                self.get_amino_acid_and_pos(
                    parts[0], used_one_letter
                )
            if aa_and_pos:
                aa_start = aa_and_pos[0]
                pos_start = aa_and_pos[1]
                used_one_letter = aa_and_pos[2]

        return aa_start, aa_end, pos_start, pos_end, used_one_letter

    def get_positions_deleted(self, parts: List) -> Optional[Tuple[str, str]]:
        """Return position(s) deleted for transcript and genomic references.

        :param List parts: Tokenized input string
        :return: Start position deleted and end position deleted
        """
        if "_" in parts[0] and parts[0].count("_") == 1:
            positions = self.get_valid_digits(parts[0])
            if not positions:
                return None
            start_pos_del, end_pos_del = positions
            if start_pos_del > end_pos_del:
                return None
        else:
            start_pos_del = parts[0]
            end_pos_del = None
            if not start_pos_del.isdigit():
                return None
        return start_pos_del, end_pos_del

    def get_transcript_genomic_inserted_sequence(
        self, parts: List
    ) -> Optional[Tuple[Union[str, int], Union[str, int]]]:
        """Return inserted sequence for transcript and genomic references.

        :param List parts: Tokenized input string
        :return: Start inserted sequence and end inserted sequence
        """
        # Check inserted sequences
        if "_" in parts[1] and parts[1].count("_") == 1:
            # Replaced by sequence positions
            inserted_sequences = self.get_valid_digits(parts[1])
            if not inserted_sequences:
                return None
            inserted_sequence1, inserted_sequence2 = inserted_sequences
            if inserted_sequence1 > inserted_sequence2:
                return None
        else:
            # Replaced by nucleotides
            inserted_sequence1 = self.get_sequence(parts[1])
            inserted_sequence2 = None
        return inserted_sequence1, inserted_sequence2

    def get_sequence(self, part: str) -> Optional[str]:
        """Return validated sequence for transcript and genomic references.

        :param str part: Sequence to validate
        :return: Sequence of nucleotides
        """
        for char in part:
            if char.upper() not in self.nucleotide_cache.base_nucleotides:
                return None
        return part.upper()

    def get_valid_digits(self, part: str) -> Optional[Tuple[str, str]]:
        """Return valid digits after splitting on `_`.

        :param str part: Range of digits
        :return: Digits represented as strings
        """
        digits = part.split("_")
        digit1 = digits[0]
        digit2 = digits[1]
        if not digit1.isdigit() or not digit2.isdigit():
            return None
        return digit1, digit2
