"""Module for Loss of Function Tokenization."""
from .basic_regex_tokenizer import BasicRegexTokenizer


class LossOfFunction(BasicRegexTokenizer):
    """Loss of Function tokenizer class."""

    def pattern(self) -> str:
        """Return regex for loss of function."""
        return r'\bLOSS[ -]OF[ -]FUNCTION\b'

    def token_type(self) -> str:
        """Return loss of function token type."""
        return 'LossOfFunction'
