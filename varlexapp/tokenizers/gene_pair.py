from .tokenizer import Tokenizer
from .gene_symbol import GeneSymbol
from ..models import Token

class GenePair(Tokenizer):
    def __init__(self, gene_cache):
        self.__gene_cache = gene_cache
        self.__gene_matcher = GeneSymbol(gene_cache)

    def match(self, input_string):
        potential_genes = input_string.split('-')
        if (len(potential_genes) == 2 and
                self.__gene_matcher.match(potential_genes[0].upper()) and
                self.__gene_matcher.match(potential_genes[1].upper())):
            return Token(input_string, 'GenePair')
        else:
            return None
