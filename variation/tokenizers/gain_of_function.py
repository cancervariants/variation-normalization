"""Module for gain of function tokenization."""
from .basic_regex_tokenizer import BasicRegexTokenizer


class GainOfFunction(BasicRegexTokenizer):
    """Gain of function tokenizer class."""

    def pattern(self) -> str:
        """Return gain of function regex."""
        return r'\bGAIN[ -]OF[ -]FUNCTION\b'

    def token_type(self) -> str:
        """Return gain of function token type."""
        return 'GainOfFunction'
