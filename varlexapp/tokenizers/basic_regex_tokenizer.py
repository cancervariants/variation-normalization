import re

from typing import Pattern, Optional

from .tokenizer import Tokenizer
from ..models import Token

class BasicRegexTokenizer(Tokenizer):
    def __init__(self) -> None:
        self.matcher: Pattern[str] = re.compile(self.pattern(), re.IGNORECASE)

    def match(self, input_string: str) -> Optional[Token]:
        if self.matcher.match(input_string):
            return Token(input_string, self.token_type(), input_string)
        else:
            return None

    def pattern(self) -> str:
        pass

    def token_type(self) -> str:
        pass
