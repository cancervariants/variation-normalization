from typing import Optional

from .tokenizer import Tokenizer
from .caches import GeneSymbolCache
from ..models import GeneMatchToken, TokenMatchType

class GeneSymbol(Tokenizer):
    def __init__(self, gene_cache: GeneSymbolCache) -> None:
        self.gene_cache = gene_cache

    def match(self, input_string: str) -> Optional[GeneMatchToken]:
        uppercase_input = input_string.upper()
        if uppercase_input in self.gene_cache.gene_symbols:
            return GeneMatchToken(uppercase_input, input_string, uppercase_input, TokenMatchType.SYMBOL)
        elif uppercase_input in self.gene_cache.gene_ids:
            return GeneMatchToken(self.gene_cache.gene_ids[uppercase_input], input_string, uppercase_input, TokenMatchType.ID)
        elif uppercase_input in self.gene_cache.gene_aliases:
            return GeneMatchToken(self.gene_cache.gene_aliases[uppercase_input], input_string, uppercase_input, TokenMatchType.ALIAS)
        elif uppercase_input in self.gene_cache.previous_identifiers:
            return GeneMatchToken(self.gene_cache.previous_identifiers[uppercase_input], input_string, uppercase_input, TokenMatchType.PREVIOUS)
        else:
            return None
