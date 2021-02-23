"""Module for Fusion tokenization."""
from .basic_regex_tokenizer import BasicRegexTokenizer


class Fusion(BasicRegexTokenizer):
    """Fusion tokenizer class."""

    def pattern(self) -> str:
        """Return regex for fusion."""
        return r'\bfusion(s)?\b'

    def token_type(self) -> str:
        """Return fusion token type."""
        return 'Fusion'
