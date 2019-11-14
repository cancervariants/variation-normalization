from .token_match_type import TokenMatchType

class Token:
    def __init__(self, token, token_type, input_string, match_type = TokenMatchType.UNSPECIFIED):
        self.token = token
        self.token_type = token_type
        self.match_type = match_type
        self.input_string = input_string
