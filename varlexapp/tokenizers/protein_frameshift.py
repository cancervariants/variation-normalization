import re

from .tokenizer import Tokenizer
from ..models import Token

class ProteinFrameshift(Tokenizer):
    def __init__(self, amino_acid_cache):
        self.__amino_acid_cache = amino_acid_cache
        self.__splitter = re.compile(r'\d+fs(\*\d+)?')
        self.__splitter_2 = re.compile(r'fs(\*\d+)?')

    def match(self, input_string):
        potential_protein = self.__splitter.split(input_string)
        if (len(potential_protein) == 2 and
                potential_protein[0] in self.__amino_acid_cache and
                not potential_protein[1]):  # confirming this is the end of the string
            return Token(input_string, 'ProteinFrameshift')
        else:
            potential_protein = self.__splitter_2.split(input_string)
            potential_amino_acid = re.split(r'\d+', potential_protein[0])
            if (len(potential_amino_acid) == 2) and
                potential_amino_acid[0] in self.__amino_acid_cache and
                potential_amino_acid[1] in self.__amino_acid_cache:
                return Token(input_string, 'ProteinFrameshift')
            else:
                return None