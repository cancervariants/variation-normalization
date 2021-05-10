"""A module for Amino Acid DelIns Tokenization Class."""
import re
from typing import Optional
from pydantic.error_wrappers import ValidationError
from .caches import AminoAcidCache
from .tokenizer import Tokenizer
from variant.schemas.token_response_schema import AminoAcidDelInsToken, \
    TokenMatchType


class AminoAcidDelIns(Tokenizer):
    """Class for tokenizing DelIns on the protein reference sequence."""

    def __init__(self, amino_acid_cache: AminoAcidCache) -> None:
        """Initialize the DelIns Base Class.

        :param AminoAcidCache amino_acid_cache: Valid amino acid codes.
        """
        self.amino_acid_cache = amino_acid_cache
        self.splitter = re.compile(r'delins')
        self.splitter_char_digit = re.compile("([a-zA-Z]+)([0-9]+)")
        self.parts = None

    def match(self, input_string: str) -> Optional[AminoAcidDelInsToken]:
        """Return token that match the input string."""
        if input_string is None:
            return None

        input_string = str(input_string).lower()

        if 'delins' not in input_string:
            return None

        if input_string.startswith('p.'):
            input_string = input_string[2:]

        if input_string.startswith('(') and input_string.endswith(')'):
            input_string = input_string[1:-1]

        self.parts = {
            'used_one_letter': False,
            'token': input_string,
            'input_string': input_string,
            'match_type': TokenMatchType.UNSPECIFIED.value,
            'start_aa_del': None,
            'start_pos_del': None,
            'end_aa_del': None,
            'end_pos_del': None,
            'inserted_sequence': None
        }

        parts = self.splitter.split(input_string)
        self._get_parts(parts)

        try:
            return AminoAcidDelInsToken(**self.parts)
        except ValidationError:
            return None

    def _get_parts(self, parts):
        """Get parts for DelIns.

        :param list parts: Parts of input string
        """
        if len(parts) != 2:
            return None

        # Get reference sequence
        self._get_positions_deleted(parts)
        self._get_inserted_sequence(parts)

    def _get_positions_deleted(self, parts):
        """Set position(s) deleted.

        :param list parts: Tokenized input string
        """
        # Check positions deleted
        if '_' in parts[0] and parts[0].count('_') == 1:
            aa_pos_range = parts[0].split('_')
            if len(aa_pos_range) != 2 or \
                    not aa_pos_range[0] or not aa_pos_range[1]:
                return None

            start_aa_pos_del = \
                self._get_amino_acid_and_pos(aa_pos_range[0])
            end_aa_pos_del = \
                self._get_amino_acid_and_pos(aa_pos_range[1])

            if start_aa_pos_del and end_aa_pos_del:
                self.parts['start_aa_del'] = start_aa_pos_del[0]
                self.parts['start_pos_del'] = start_aa_pos_del[1]
                self.parts['end_aa_del'] = end_aa_pos_del[0]
                self.parts['end_pos_del'] = end_aa_pos_del[1]

        else:
            aa_and_pos = \
                self._get_amino_acid_and_pos(parts[0])
            if aa_and_pos:
                self.parts['start_aa_del'] = aa_and_pos[0]
                self.parts['start_pos_del'] = aa_and_pos[1]

    def _get_amino_acid_and_pos(self, part):
        """Return amino acid and position deleted.

        :param list part: Tokenized input string
        """
        try:
            char_and_digits = self.splitter_char_digit.match(part).groups()
        except AttributeError:
            return None

        if len(char_and_digits) != 2 or len(part) != \
                (len(char_and_digits[0]) + len(char_and_digits[1])):
            return None

        aa_del = char_and_digits[0]
        if len(aa_del) == 1:
            if not self.parts['used_one_letter']:
                self.parts['used_one_letter'] = True
            aa_del = self.amino_acid_cache.amino_acid_code_conversion[aa_del.upper()]  # noqa: E501
        pos_del = char_and_digits[1]

        if not self.amino_acid_cache.__contains__(aa_del) \
                or not pos_del.isdigit():
            return None
        return aa_del.capitalize(), pos_del

    def _get_inserted_sequence(self, parts):
        """Return inserted sequence.

        :param list parts: Tokenized input string
        """
        # Check inserted sequences
        inserted_sequence = ""
        if self.parts['used_one_letter']:
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
        self.parts['inserted_sequence'] = inserted_sequence
