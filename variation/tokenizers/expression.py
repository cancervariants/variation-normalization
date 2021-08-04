"""Module for expression tokenization."""
from .basic_regex_tokenizer import BasicRegexTokenizer


class Expression(BasicRegexTokenizer):
    """Expression tokenizer class."""

    def pattern(self) -> str:
        """Return expression regex."""
        return r'\b((?<!over )expression|(?<!under )expression)\b'

    def token_type(self) -> str:
        """Return expression token type."""
        return 'Expression'
