from .token import Token

class GenePairMatchToken(Token):
    def __init__(self, token, input_string, left_gene_token, right_gene_token):
        super().__init__(token, 'GenePair', input_string)
        self.left_gene_token = left_gene_token
        self.right_gene_token = right_gene_token
