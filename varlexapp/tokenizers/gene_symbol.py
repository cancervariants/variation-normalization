from .tokenizer import Tokenizer
from ..models import Token, TokenMatchType

class GeneSymbol(Tokenizer):
    def __init__(self, gene_cache):
        self.gene_cache = gene_cache

    def match(self, input_string):
        uppercase_input = input_string.upper()
        if uppercase_input in self.gene_cache.gene_symbols:
            return Token(input_string, 'GeneSymbol', TokenMatchType.SYMBOL)
        elif uppercase_input in self.gene_cache.gene_ids:
            return Token(self.gene_cache.gene_ids[uppercase_input], 'GeneSymbol', TokenMatchType.ID)
        elif uppercase_input in self.gene_cache.gene_aliases:
            return Token(self.gene_cache.gene_aliases[uppercase_input], 'GeneSymbol', TokenMatchType.ALIAS)
        elif uppercase_input in self.gene_cache.previous_identifiers:
            return Token(self.gene_cache.previous_identifiers[uppercase_input], 'GeneSymbol', TokenMatchType.PREVIOUS)
        else:
            return None
