"""Module for basic regex tokenizer."""
import re
from typing import Pattern, Optional
from .tokenizer import Tokenizer
from variation.schemas.token_response_schema import Token, TokenMatchType


class BasicRegexTokenizer(Tokenizer):
    """The basic regex tokenizer class."""

    def __init__(self) -> None:
        """Initialize the basic regex tokenizer class."""
        self.matcher: Pattern[str] = re.compile(self.pattern(), re.IGNORECASE)

    def match(self, input_string: str) -> Optional[Token]:
        """Return the token the input string matches."""
        if self.matcher.match(input_string):
            return Token(
                token=input_string,
                token_type=self.token_type(),
                input_string=input_string,
                match_type=TokenMatchType.UNSPECIFIED
            )
        else:
            return None

    def pattern(self) -> str:
        """TODO."""
        pass

    def token_type(self) -> str:
        """TODO."""
        pass
