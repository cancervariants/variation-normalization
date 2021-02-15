"""Module for tokenizing Protein Frameshift."""
import re
from typing import Optional
from .tokenizer import Tokenizer
from .caches import AminoAcidCache
from varlexapp.schemas.token_response_schema import Token


class ProteinFrameshift(Tokenizer):
    """The Protein Frameshift Tokenizer class."""

    def __init__(self, amino_acid_cache: AminoAcidCache) -> None:
        """Initialize the Protein Frameshift Tokenizer class."""
        self.__amino_acid_cache = amino_acid_cache
        self.__splitter = re.compile(r'\d+fs(\*\d+)?')

    def match(self, input_string: str) -> Optional[Token]:
        """Return token match for Protein Frameshift."""
        potential_protein = self.__splitter.split(input_string)
        if (len(potential_protein) == 2 and
                potential_protein[0] in self.__amino_acid_cache and
                not potential_protein[1]):
            return Token(input_string, 'ProteinFrameshift', input_string)
        else:
            return None
