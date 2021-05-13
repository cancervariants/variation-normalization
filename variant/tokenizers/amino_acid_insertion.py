"""A module for Amino Acid Insertion Tokenization Class."""
import re
from typing import Optional
from pydantic.error_wrappers import ValidationError
from .caches import AminoAcidCache
from .tokenizer import Tokenizer
from .amino_acid_base import AminoAcidBase
from variant.schemas.token_response_schema import AminoAcidInsertionToken, \
    TokenMatchType


class AminoAcidInsertion(Tokenizer):
    """Class for tokenizing Insertions on the protein reference sequence."""

    def __init__(self, amino_acid_cache: AminoAcidCache) -> None:
        """Initialize the Amino Acid Deletion Class.

        :param AminoAcidCache amino_acid_cache: Valid amino acid codes.
        """
        self.splitter = re.compile(r'ins')
        self.splitter_char_digit = re.compile("([a-zA-Z]+)([0-9]+)")
        self.parts = None
        self.amino_acid_base = AminoAcidBase(amino_acid_cache)

    def match(self, input_string: str) -> Optional[AminoAcidInsertionToken]:
        """Return token that match the input string."""
        if input_string is None:
            return None

        self.parts = {
            'used_one_letter': False,
            'token': input_string,
            'input_string': input_string,
            'match_type': TokenMatchType.UNSPECIFIED.value,
            'start_aa_flank': None,
            'start_pos_flank': None,
            'end_aa_flank': None,
            'end_pos_flank': None
        }

        input_string = str(input_string).lower()

        if 'c.' in input_string or 'g.' in input_string:
            return None

        if input_string.startswith('p.'):
            input_string = input_string[2:]

        if input_string.startswith('(') and input_string.endswith(')'):
            input_string = input_string[1:-1]

        if 'ins' not in input_string:
            return None

        parts = self.splitter.split(input_string)
        self._get_parts(parts)

        try:
            return AminoAcidInsertionToken(**self.parts)
        except ValidationError:
            return None

    def _get_parts(self, parts):
        """Get parts for Amino Acid Insertion.

        :param list parts: Parts of input string
        """
        if len(parts) != 2:
            return

        range_aa_pos = self.amino_acid_base.get_possible_range(parts)
        if range_aa_pos:
            self.parts['start_aa_flank'] = range_aa_pos[0]
            self.parts['end_aa_flank'] = range_aa_pos[1]
            self.parts['start_pos_flank'] = range_aa_pos[2]
            self.parts['end_pos_flank'] = range_aa_pos[3]
            self.parts['used_one_letter'] = range_aa_pos[4]
        self.parts['inserted_sequence'] = \
            self.amino_acid_base.get_inserted_sequence(
                parts, self.parts['used_one_letter'])
