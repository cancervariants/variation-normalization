from .token import Token
from . import GeneMatchToken

class GenePairMatchToken(Token):
    def __init__(self, token: str, input_string: str, left_gene_token: GeneMatchToken, right_gene_token: GeneMatchToken) -> None:
        super().__init__(token, 'GenePair', input_string)
        self.left_gene_token = left_gene_token
        self.right_gene_token = right_gene_token
