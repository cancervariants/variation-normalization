"""Module for Gene Symbol tokenization."""
from typing import Optional

from gene.query import QueryHandler as GeneQueryHandler

from variation.schemas.token_response_schema import GeneMatchToken
from variation.tokenizers.tokenizer import Tokenizer


class GeneSymbol(Tokenizer):
    """Class for gene symbol tokenization."""

    def __init__(self, gene_normalizer: GeneQueryHandler) -> None:
        """Initialize the gene symbol tokenizer class.

        :param QueryHandler gene_normalizer: Instance to gene normalizer
            QueryHandler
        """
        self.gene_normalizer = gene_normalizer
        self._gene_cache = dict()

    def match(self, input_string: str) -> Optional[GeneMatchToken]:
        """Return token matches from input string."""
        upper_input = input_string.upper()
        if upper_input in self._gene_cache:
            return self._gene_cache[upper_input]

        norm_resp = self.gene_normalizer.normalize(upper_input)
        norm_match_type = norm_resp.match_type

        if norm_match_type != 0:
            gene_descriptor = norm_resp.gene_descriptor
            label = gene_descriptor.label
            gene_match_token = GeneMatchToken(
                token=label,
                input_string=input_string,
                matched_value=label,
                gene_descriptor=gene_descriptor
            )
            self._gene_cache[upper_input] = gene_match_token
            return gene_match_token
        else:
            return None
