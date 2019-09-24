#\\b{}\\d+(_{}\\d+)?((DEL|del)|(INS|ins)|>)+{}*\\b
#protein digits (protein digits)? (del|ins|>) protein
#why "one or more" on the del/ins/>
import re

from .tokenizer import Tokenizer
from .token import Token

class ProteinFrameshift(Tokenizer):
    def __init__(self, amino_acid_cache):
        self.__amino_acid_cache = amino_acid_cache
        self.__splitter = re.compile(r'\d+fs(\*\d+)?')

    def match(self, input_string):
        potential_protein = self.__splitter.split(input_string)
        conditions = (
                len(potential_protein) == 2,
                potential_protein[0] in self.__amino_acid_cache,
                not potential_protein[1]
        )

        if all(conditions):
            return Token(input_string, 'ProteinFrameshift')

