"""Module for Protein Termination tokenizer."""
import re
from typing import Optional
from .caches import AminoAcidCache
from .tokenizer import Tokenizer
from variation.schemas.token_response_schema import Token, TokenMatchType


class ProteinTermination(Tokenizer):
    """The protein termination tokenizer class."""

    def __init__(self, amino_acid_cache: AminoAcidCache) -> None:
        """Initialize the Protein Termination tokenizer class."""
        self.__amino_acid_cache = amino_acid_cache
        self.__splitter = re.compile(r'\d+')

    def match(self, input_string: str) -> Optional[Token]:
        """Find Protein Termination token matches."""
        potential_protein = self.__splitter.split(input_string)
        if (len(potential_protein) == 2 and
                potential_protein[0] in self.__amino_acid_cache and
                potential_protein[1].upper() == 'X'):
            return Token(
                token=input_string,
                token_type='ProteinTermination',
                input_string=input_string,
                match_type=TokenMatchType.UNSPECIFIED
            )
        else:
            return None
