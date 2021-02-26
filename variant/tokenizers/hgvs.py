"""Module for HGVS tokenization."""
from typing import Optional
from .tokenizer import Tokenizer
from variant.schemas.token_response_schema import Token, TokenMatchType
from hgvs.parser import Parser
from hgvs.exceptions import HGVSParseError, HGVSInvalidVariantError
from hgvs.validator import IntrinsicValidator


class HGVS(Tokenizer):
    """The HGVS tokenizer class."""

    def __init__(self) -> None:
        """Initialize the HGVS tokenizer class."""
        self.parser = Parser()
        self.validator = IntrinsicValidator()

    def match(self, input_string: str) -> Optional[Token]:
        """Return token matches from input string."""
        try:
            variant = self.parser.parse_hgvs_variant(input_string)
            self.validator.validate(variant)
            return Token(
                token=input_string,
                token_type='HGVS',
                input_string=input_string,
                match_type=TokenMatchType.UNSPECIFIED
            )
        except (HGVSParseError, HGVSInvalidVariantError):
            return None
