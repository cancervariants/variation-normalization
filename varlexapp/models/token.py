from .token_match_type import TokenMatchType

class Token:
    def __init__(self, token, token_type, match_type = TokenMatchType.UNSPECIFIED):
        self.token = token
        self.token_type = token_type
        self.match_type = match_type
