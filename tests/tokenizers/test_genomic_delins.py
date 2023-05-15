"""A module for testing Genomic DelIns tokenization."""
import unittest

from variation.tokenizers import GenomicDelIns
from variation.schemas.token_response_schema import TokenType
from .tokenizer_base import TokenizerBase


class TestGenomicDelInsTokenizer(TokenizerBase, unittest.TestCase):
    """A class for testing Genomic DelIns Tokenization."""

    def tokenizer_instance(self):
        """Return Genomic DelIns instance."""
        return GenomicDelIns()

    def token_type(self):
        """Return genomic delins token type."""
        return TokenType.GENOMIC_DELINS

    def fixture_name(self):
        """Return the fixture name for genomic delins."""
        return "genomic_delins"
