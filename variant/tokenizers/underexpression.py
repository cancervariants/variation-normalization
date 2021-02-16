"""Module for under expression tokenization."""
from .basic_regex_tokenizer import BasicRegexTokenizer


class UnderExpression(BasicRegexTokenizer):
    """Under expression tokenizer class."""

    def pattern(self) -> str:
        """Return under expression regex."""
        return r'\bunder ?expression\b'

    def token_type(self) -> str:
        """Return under expression token type."""
        return 'UnderExpression'
