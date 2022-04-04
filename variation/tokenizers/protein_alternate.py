"""Module for Protein Alternate classification."""
import re
from typing import Optional

from variation.schemas.token_response_schema import Token, TokenMatchType
from .tokenizer import Tokenizer
from .caches import AminoAcidCache


class ProteinAlternate(Tokenizer):
    """The Protein Alternate Tokenization class."""

    def __init__(self, amino_acid_cache: AminoAcidCache) -> None:
        """Initialize the Protein Alternate Tokenizer class."""
        self.__amino_acid_cache = amino_acid_cache
        self.__splitter = re.compile(r"\d+")

    def match(self, input_string: str) -> Optional[Token]:
        """Return Protein Alternate tokens if input string matches."""
        potential_protein = self.__splitter.split(input_string)
        if (len(potential_protein) == 2 and  # noqa: W504
                potential_protein[0] in self.__amino_acid_cache and   # noqa: W504
                not potential_protein[1]):
            return Token(
                token=potential_protein[0],
                token_type="ProteinAlternate",
                input_string=input_string,
                match_type=TokenMatchType.UNSPECIFIED
            )
        else:
            return None
