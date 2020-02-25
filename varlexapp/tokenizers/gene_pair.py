from typing import Optional

from .tokenizer import Tokenizer
from .gene_symbol import GeneSymbol
from ..models import GenePairMatchToken
from .caches import GeneSymbolCache

class GenePair(Tokenizer):
    def __init__(self, gene_cache: GeneSymbolCache) -> None:
        self.__gene_cache = gene_cache
        self.__gene_matcher = GeneSymbol(gene_cache)

    def match(self, input_string: str) -> Optional[GenePairMatchToken]:
        potential_genes = input_string.split('-')
        if len(potential_genes) == 2:
            first_match = self.__gene_matcher.match(potential_genes[0].upper())
            second_match = self.__gene_matcher.match(potential_genes[1].upper())
            if (first_match and second_match):
                return GenePairMatchToken(f"{first_match.token}-{second_match.token}", input_string, first_match, second_match)
        return None
