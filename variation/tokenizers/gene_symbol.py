"""Module for Gene Symbol tokenization."""
from typing import Optional
from gene.schemas import MatchType
from .tokenizer import Tokenizer
from variation.schemas.token_response_schema import GeneMatchToken, \
    TokenMatchType
from variation import GENE_NORMALIZER


class GeneSymbol(Tokenizer):
    """Class for gene symbol tokenization."""

    def __init__(self) -> None:
        """Initialize the gene symbol tokenizer class."""
        self.gene_normalizer = GENE_NORMALIZER

    def match(self, input_string: str) -> Optional[GeneMatchToken]:
        """Return token matches from input string."""
        norm_resp = self.gene_normalizer.normalize(input_string)

        norm_match_type = norm_resp['match_type']
        if norm_match_type == MatchType.CONCEPT_ID:
            match_type = TokenMatchType.ID
        elif norm_match_type == MatchType.SYMBOL:
            match_type = TokenMatchType.SYMBOL
        elif norm_match_type == MatchType.ALIAS:
            match_type = TokenMatchType.ALIAS
        elif norm_match_type == MatchType.PREV_SYMBOL:
            match_type = TokenMatchType.PREVIOUS
        else:
            match_type = TokenMatchType.UNSPECIFIED

        if match_type != TokenMatchType.UNSPECIFIED:
            label = norm_resp['gene_descriptor']['label']
        else:
            label = input_string.upper()

        if norm_match_type != 0:
            return GeneMatchToken(
                token=label,
                input_string=input_string,
                matched_value=label,
                match_type=match_type
            )
        else:
            return None
