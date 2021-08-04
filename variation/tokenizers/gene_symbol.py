"""Module for Gene Symbol tokenization."""
from typing import Optional
from .tokenizer import Tokenizer
from .caches import GeneSymbolCache
from variation.schemas.token_response_schema import GeneMatchToken, \
    TokenMatchType


class GeneSymbol(Tokenizer):
    """Class for gene symbol tokenization."""

    def __init__(self, gene_cache: GeneSymbolCache) -> None:
        """Initialize the gene symbol tokenizer class."""
        self.gene_cache = gene_cache

    def match(self, input_string: str) -> Optional[GeneMatchToken]:
        """Return token matches from input string."""
        uppercase_input = input_string.upper()
        if uppercase_input in self.gene_cache.gene_symbols:
            return GeneMatchToken(
                token=uppercase_input,
                input_string=input_string,
                matched_value=uppercase_input,
                match_type=TokenMatchType.SYMBOL
            )
        elif uppercase_input in self.gene_cache.gene_ids:
            return GeneMatchToken(
                token=self.gene_cache.gene_ids[uppercase_input],
                input_string=input_string,
                matched_value=uppercase_input,
                match_type=TokenMatchType.ID
            )
        elif uppercase_input in self.gene_cache.gene_aliases:
            return GeneMatchToken(
                token=self.gene_cache.gene_aliases[uppercase_input],
                input_string=input_string,
                matched_value=uppercase_input,
                match_type=TokenMatchType.ALIAS
            )
        elif uppercase_input in self.gene_cache.previous_identifiers:
            return GeneMatchToken(
                token=self.gene_cache.previous_identifiers[uppercase_input],
                input_string=input_string,
                matched_value=uppercase_input,
                match_type=TokenMatchType.PREVIOUS
            )
        else:
            return None
