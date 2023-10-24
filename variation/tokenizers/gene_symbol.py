"""Module for Gene Symbol tokenization."""
from typing import Optional

from gene.query import QueryHandler as GeneQueryHandler

from variation.schemas.token_response_schema import GeneToken
from variation.tokenizers.tokenizer import Tokenizer


class GeneSymbol(Tokenizer):
    """Class for gene symbol tokenization."""

    def __init__(self, gene_normalizer: GeneQueryHandler) -> None:
        """Initialize the gene symbol tokenizer class.

        :param gene_normalizer: Instance to gene normalizer QueryHandler
        """
        self.gene_normalizer = gene_normalizer

    def match(self, input_string: str) -> Optional[GeneToken]:
        """Return tokens that are genes

        :param input_string: Input string
        :return: GeneToken if match was found
        """
        norm_resp = self.gene_normalizer.normalize(input_string)
        norm_match_type = norm_resp.match_type

        if norm_match_type != 0:
            gene = norm_resp.gene
            label = gene.label
            gene_match_token = GeneToken(
                token=label,
                input_string=input_string,
                matched_value=label,
                gene=gene,
            )
            return gene_match_token

        return None
