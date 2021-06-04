"""Module for over expression tokenization."""
from .basic_regex_tokenizer import BasicRegexTokenizer


class OverExpression(BasicRegexTokenizer):
    """Over Expression tokenizer class."""

    def pattern(self) -> str:
        """Return over expression regex."""
        return r'\bover ?expression\b'

    def token_type(self) -> str:
        """Return over expression token type."""
        return 'OverExpression'
