"""A module for Amino Acid DelIns Tokenization Class."""
import re
from typing import Optional
from pydantic.error_wrappers import ValidationError
from .caches import AminoAcidCache, NucleotideCache
from .tokenizer import Tokenizer
from .tokenize_base import TokenizeBase
from variation.schemas.token_response_schema import AminoAcidDelInsToken, \
    TokenMatchType


class AminoAcidDelIns(Tokenizer):
    """Class for tokenizing DelIns on the protein reference sequence."""

    def __init__(self, amino_acid_cache: AminoAcidCache,
                 nucleotide_cache: NucleotideCache) -> None:
        """Initialize the Amino Acid DelIns Class.

        :param AminoAcidCache amino_acid_cache: Valid amino acid codes
        :param NucleotideCache nucleotide_cache: Valid nucleotides
        """
        self.splitter = re.compile('delins')
        self.parts = None
        self.tokenize_base = TokenizeBase(amino_acid_cache, nucleotide_cache)

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
        range_aa_pos = self.tokenize_base.get_aa_pos_range(parts)
        if range_aa_pos:
            self.parts['start_aa_del'] = range_aa_pos[0]
            self.parts['end_aa_del'] = range_aa_pos[1]
            self.parts['start_pos_del'] = range_aa_pos[2]
            self.parts['end_pos_del'] = range_aa_pos[3]
            self.parts['used_one_letter'] = range_aa_pos[4]
        self.parts['inserted_sequence'] = \
            self.tokenize_base.get_protein_inserted_sequence(
                parts, self.parts['used_one_letter'])
