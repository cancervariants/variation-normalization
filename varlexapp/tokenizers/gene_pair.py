from .tokenizer import Tokenizer
from .gene_symbol import GeneSymbol
from ..models import Token

class GenePair(Tokenizer):
    def __init__(self, gene_cache):
        self.__gene_cache = gene_cache
        self.__gene_matcher = GeneSymbol(gene_cache)

    def match(self, input_string):
        potential_genes = input_string.split('-')
        if len(potential_genes) == 2:
            first_match = self.__gene_matcher.match(potential_genes[0].upper())
            second_match = self.__gene_matcher.match(potential_genes[1].upper())
            if (first_match and second_match):
                return Token(f"{first_match.token}-{second_match.token}", 'GenePair', input_string)
        return None
