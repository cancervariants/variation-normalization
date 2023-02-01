"""Module for Protein Alternate classification."""
import re
from typing import Optional

from bioutils.sequences import aa3_to_aa1_lut

from variation.schemas.token_response_schema import Token, TokenMatchType
from .tokenizer import Tokenizer


class ProteinAlternate(Tokenizer):
    """The Protein Alternate Tokenization class."""

    def __init__(self) -> None:
        """Initialize the Protein Alternate Tokenizer class."""
        self.__splitter = re.compile(r"\d+")

    def match(self, input_string: str) -> Optional[Token]:
        """Return Protein Alternate tokens if input string matches."""
        potential_protein = self.__splitter.split(input_string)
        if all((len(potential_protein) == 2,
                potential_protein[0] in aa3_to_aa1_lut,
                not potential_protein[1])):
            return Token(
                token=potential_protein[0],
                token_type="ProteinAlternate",
                input_string=input_string,
                match_type=TokenMatchType.UNSPECIFIED
            )
        else:
            return None
