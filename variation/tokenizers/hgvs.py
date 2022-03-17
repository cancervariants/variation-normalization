"""Module for HGVS tokenization."""
from typing import Optional

from hgvs.parser import Parser
from hgvs.exceptions import HGVSParseError, HGVSInvalidVariantError
from hgvs.validator import IntrinsicValidator

from variation.tokenizers.reference_sequence import REFSEQ_PREFIXES
from variation.schemas.token_response_schema import Token, TokenMatchType, \
    Nomenclature
from .tokenizer import Tokenizer
from .locus_reference_genomic import LocusReferenceGenomic


class HGVS(Tokenizer):
    """The HGVS tokenizer class."""

    def __init__(self) -> None:
        """Initialize the HGVS tokenizer class."""
        self.parser = Parser()
        self.validator = IntrinsicValidator()

    def match(self, input_string: str) -> Optional[Token]:
        """Return token matches from input string."""
        valid_prefix = False
        for prefix in REFSEQ_PREFIXES:
            if input_string.upper().startswith(prefix):
                valid_prefix = True
                break

        if not valid_prefix:
            return None

        try:
            variation = self.parser.parse_hgvs_variant(input_string)
            self.validator.validate(variation)
            if input_string[:3].upper().startswith("LRG"):
                lrg_token = LocusReferenceGenomic().match(
                    input_string.split(":")[0]
                )
                lrg_token.input_string = input_string
                return lrg_token
            else:
                return Token(
                    token=input_string,
                    token_type="HGVS",
                    input_string=input_string,
                    match_type=TokenMatchType.UNSPECIFIED,
                    nomenclature=Nomenclature.HGVS
                )
        except (HGVSParseError, HGVSInvalidVariantError, AttributeError):
            return None
