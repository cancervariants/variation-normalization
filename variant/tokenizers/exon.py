"""Module for exon tokenization."""
from .basic_regex_tokenizer import BasicRegexTokenizer


class Exon(BasicRegexTokenizer):
    """Exon tokenizer class."""

    def pattern(self) -> str:
        """Return exon regex."""
        return r'\bexon \d+\b'

    def token_type(self) -> str:
        """Return exon token type."""
        return 'Exon'
