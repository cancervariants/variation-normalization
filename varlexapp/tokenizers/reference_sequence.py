"""A module for Reference Sequence Tokenization."""
from typing import Optional
from .tokenizer import Tokenizer
from ..models import Token

REFSEQ_PREFIXES = ["NC_", "NT_", "NW_", "NG_", "NM_", "NR_", "NP_", "LRG_"]


class ReferenceSequence(Tokenizer):
    """Class for Reference Sequence."""

    def match(self, input_string: str) -> Optional[Token]:
        """Return a Reference Sequence match if it exists.

        :param str input_string: The input string to match
        """
        if input_string[:3] in REFSEQ_PREFIXES:
            return Token(input_string, 'ReferenceSequence', input_string)
        return None
