"""Module for Deletion Tokenization."""
from .basic_regex_tokenizer import BasicRegexTokenizer


class Deletion(BasicRegexTokenizer):
    """Deletion Tokenizer class."""

    def pattern(self) -> str:
        """Return regex for Deletion."""
        return r'\b(del|deletion|(copy number)? ?loss)\b'

    def token_type(self) -> str:
        """Return deletion token type."""
        return 'Deletion'
