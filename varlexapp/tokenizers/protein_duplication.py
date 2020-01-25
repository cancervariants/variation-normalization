import re

from .tokenizer import Tokenizer
from ..models import Token

class ProteinDuplication(Tokenizer):
    def __init__(self, amino_acid_cache):
        self.__amino_acid_cache = amino_acid_cache
        self.__splitter = re.compile(r'\d+_?')

    def match(self, input_string):
        potential_proteins = self.__splitter.split(input_string)
        if (len(potential_proteins) == 3 and
                potential_proteins[0] in self.__amino_acid_cache and
                potential_proteins[1] in self.__amino_acid_cache) and 
                potential_proteins[2].upper() == "DUP"
            return Token(input_string, 'ProteinDuplication')  
        else:
            return None

