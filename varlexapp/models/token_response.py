from typing import List

from . import Token

class TokenResponse:
    def __init__(self, search_term: str, tokens: List[Token]) -> None:
        self.search_term = search_term
        self.tokens = tokens
