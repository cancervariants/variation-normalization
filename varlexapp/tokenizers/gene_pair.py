from .tokenizer import Tokenizer
from ..models import Token

class GenePair(Tokenizer):
    def __init__(self, gene_cache):
        self.__gene_cache = gene_cache

    def match(self, input_string):
        potential_genes = input_string.split('-')
        if (len(potential_genes) == 2 and
                potential_genes[0].upper() in self.__gene_cache and
                potential_genes[1].upper() in self.__gene_cache):
            return Token(input_string, 'GenePair')
        else:
            return None
