"""A module for Amino Acid DelIns Tokenization Class."""
import re
from typing import Optional
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
        """Return tokens that match the input string."""
        if input_string is None:
            return None

        input_string = str(input_string).lower()

        if input_string.startswith('(') and input_string.endswith(')'):
            input_string = input_string[1:-1]

        if 'delins' not in input_string or not input_string.startswith('p.'):
            return None

        self.parts = {
            'start_aa_del': None,
            'start_pos_del': None,
            'end_aa_del': None,
            'end_pos_del': None,
            'inserted_sequence': None
        }

        parts = self.splitter.split(input_string)
        self._get_parts(parts)

        params = {
            'token': input_string,
            'input_string': input_string,
            'match_type': TokenMatchType.UNSPECIFIED.value,
            'start_aa_del': self.parts['start_aa_del'],
            'start_pos_del': self.parts['start_pos_del'],
            'end_aa_del': self.parts['end_aa_del'],
            'end_pos_del': self.parts['end_pos_del'],
            'inserted_sequence': self.parts['inserted_sequence']
        }
        return AminoAcidDelInsToken(**params)

    def _get_parts(self, parts):
        """Get parts for DelIns.

        :param list parts: Parts of input string
        """
        if len(parts) != 2:
            return

        # Get reference sequence
        deleted_amino_acids = parts[0][:1]  # noqa: F841
        parts[0] = parts[0][2:]
        self._get_positions_deleted(parts)
        self._get_inserted_sequence(parts)

    def _get_positions_deleted(self, parts):
        """Return position(s) deleted."""
        # Check positions deleted
        if '_' in parts[0] and parts[0].count('_') == 1:
            positions = self._get_valid_digits(parts[0])
            if not positions:
                return
            start_pos_del, end_pos_del = positions
            if start_pos_del > end_pos_del:
                return
        else:
            char_and_digits = self.splitter_char_digit.match(parts[0]).groups()

            if len(char_and_digits) != 2:
                return None

            start_aa_del = char_and_digits[0]
            start_pos_del = char_and_digits[1]

            if not self.amino_acid_cache.__contains__(start_aa_del):
                return None

            if not start_pos_del.isdigit():
                return None

            self.parts['start_aa_del'] = start_aa_del.capitalize()
            self.parts['start_pos_del'] = start_pos_del

    def _get_inserted_sequence(self, parts):
        """Return inserted sequence."""
        # Check inserted sequences
        inserted_sequence = ""
        for i in range(0, len(parts[1]), 3):
            aa = parts[1][i:i + 3]
            if len(aa) != 3:
                return None
            if not self.amino_acid_cache.__contains__(aa):
                return None
            inserted_sequence += aa.capitalize()

        self.parts['inserted_sequence'] = inserted_sequence
