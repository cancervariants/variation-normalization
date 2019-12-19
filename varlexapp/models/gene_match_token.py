from .token import Token

class GeneMatchToken(Token):
    def __init__(self, token, input_string, matched_value, match_type):
        super().__init__(token, 'GeneSymbol', input_string, match_type)
        self.matched_value = matched_value

