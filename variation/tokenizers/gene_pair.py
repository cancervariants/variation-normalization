"""Module for tokenizing gene pair."""
from typing import Optional

from gene.query import QueryHandler as GeneQueryHandler

from variation.schemas.token_response_schema import GenePairMatchToken, \
    TokenMatchType
from .tokenizer import Tokenizer
from .gene_symbol import GeneSymbol


class GenePair(Tokenizer):
    """The gene pair class."""

    def __init__(self) -> None:
        """Initialize the Gene Pair class."""
        self.__gene_matcher = GeneSymbol(GeneQueryHandler())

    def match(self, input_string: str) -> Optional[GenePairMatchToken]:
        """Return tokens that match the input string."""
        potential_genes = input_string.split("-")
        if len(potential_genes) == 2:
            first_match = self.__gene_matcher.match(potential_genes[0].upper())
            second_match = \
                self.__gene_matcher.match(potential_genes[1].upper())
            if first_match and second_match:
                return GenePairMatchToken(
                    token=f"{first_match.token}-{second_match.token}",
                    match_type=TokenMatchType.UNSPECIFIED,
                    input_string=input_string,
                    left_gene_token=first_match,
                    right_gene_token=second_match
                )
        return None
