"""Module for wild type tokenization."""
from .basic_regex_tokenizer import BasicRegexTokenizer


class WildType(BasicRegexTokenizer):
    """Wild type tokenizer class."""

    def pattern(self):
        """Return wild type regex."""
        return r'\b(wild type)|wt\b'

    def token_type(self):
        """Return wild type token type."""
        return 'WildType'
