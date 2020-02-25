from .token import Token
from . import TokenMatchType

class GeneMatchToken(Token):
    def __init__(self, gene_symbol: str, input_string: str, matched_value: str, match_type: TokenMatchType) -> None:
        super().__init__(gene_symbol, 'GeneSymbol', input_string, match_type)
        self.matched_value = matched_value

