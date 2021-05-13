"""Module for Amino Acid commonly used methods during tokenization."""
from typing import Optional, Tuple

from .caches import AminoAcidCache
import re


class AminoAcidBase:
    """The Amino Acid Base Class."""

    def __init__(self, amino_acid_cache: AminoAcidCache) -> None:
        """Initialize the DelIns Base Class.

        :param AminoAcidCache amino_acid_cache: Valid amino acid codes.
        """
        self.amino_acid_cache = amino_acid_cache
        self.splitter_char_digit = re.compile("([a-zA-Z]+)([0-9]+)")

    def get_amino_acid_and_pos(self, part, used_one_letter)\
            -> Optional[Tuple[str, int, bool]]:
        """Return amino acid and position.

        :param list part: Tokenized input string
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
            aa = self.amino_acid_cache.amino_acid_code_conversion[aa.upper()]
        pos = char_and_digits[1]

        if not self.amino_acid_cache.__contains__(aa) \
                or not pos.isdigit():
            return None
        return aa.capitalize(), pos, used_one_letter

    def get_inserted_sequence(self, parts, used_one_letter) -> Optional[str]:
        """Return inserted sequence.

        :param list parts: Tokenized input string
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
                    aa = self.amino_acid_cache.amino_acid_code_conversion[aa.upper()]  # noqa: E501
                except KeyError:
                    return None
                else:
                    inserted_sequence += aa

        else:
            for i in range(0, len(parts[1]), 3):
                aa = parts[1][i:i + 3]
                if len(aa) != 3 or not self.amino_acid_cache.__contains__(aa):
                    if aa != 'ter':
                        return None
                inserted_sequence += aa.capitalize()

        if inserted_sequence == '':
            return None
        return inserted_sequence

    def get_possible_range(self, parts)\
            -> Optional[Tuple[str, str, str, int, bool]]:
        """Get amino acid(s) and positions(s).

        :param list parts: Tokenized input string
        :return: Beginning AA, End AA,  Beginning position, End position,
            Whether or not one letter code was used
        """
        aa_start = None
        aa_end = None
        pos_start = None
        pos_end = None
        used_one_letter = False

        if '_' in parts[0] and parts[0].count('_') == 1:
            aa_pos_range = parts[0].split('_')
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
