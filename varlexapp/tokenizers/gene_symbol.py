from .tokenizer import Tokenizer
from .token import Token

class GeneSymbol(Tokenizer):
    def __init__(self, gene_cache):
        self.__gene_cache = gene_cache

    def match(self, input_string):
        if input_string.upper() in self.__gene_cache:
            return Token(input_string, 'GeneSymbol')
        else:
            return None
