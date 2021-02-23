"""Module for amplification tokenization."""
from .basic_regex_tokenizer import BasicRegexTokenizer


class Amplification(BasicRegexTokenizer):
    """Amplification tokenizer class."""

    def pattern(self) -> str:
        """Return amplification regex."""
        return r'\b(amp|amplification|(copy number)? ?gain)\b'

    def token_type(self) -> str:
        """Return amplification token type."""
        return 'Amplification'
