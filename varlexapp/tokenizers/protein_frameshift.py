import re

from typing import Optional

from .tokenizer import Tokenizer
from .caches import AminoAcidCache
from ..models import Token

class ProteinFrameshift(Tokenizer):
    def __init__(self, amino_acid_cache: AminoAcidCache) -> None:
        self.__amino_acid_cache = amino_acid_cache
        self.__splitter = re.compile(r'\d+fs(\*\d+)?')

    def match(self, input_string: str) -> Optional[Token]:
        potential_protein = self.__splitter.split(input_string)
        if (len(potential_protein) == 2 and
                potential_protein[0] in self.__amino_acid_cache and
                not potential_protein[1]):
            return Token(input_string, 'ProteinFrameshift', input_string)
        else:
            return None

