"""Module for Protein DelIns Tokenizer."""
import re
from typing import Pattern, Optional
from .tokenizer import Tokenizer
from .caches import AminoAcidCache
from variation.schemas.token_response_schema import Token, TokenMatchType


class ProteinDelIns(Tokenizer):
    """The Protein DelIns Tokenizer class."""

    def __init__(self, amino_acid_cache: AminoAcidCache) -> None:
        """Initialize the Protein DelIns Tokenizer."""
        self.__amino_acid_cache = amino_acid_cache
        # TODO: Fix this. Below is protein fs regex
        self.__splitter: Pattern[str] = re.compile(r'\d+fs(\*\d+)?')

    def match(self, input_string: str) -> Optional[Token]:
        """Return token matches for Protein DelIns."""
        potential_protein = self.__splitter.split(input_string)
        conditions = (
            len(potential_protein) == 2,
            potential_protein[0] in self.__amino_acid_cache,
            not potential_protein[1]
        )

        if all(conditions):
            return Token(
                token=input_string,
                token_type='ProteinDelins',
                match_type=TokenMatchType.UNSPECIFIED.value,
                input_string=input_string
            )
        else:
            return None
