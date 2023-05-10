"""A module for testing Coding DNA DelIns tokenization."""
import unittest

from variation.tokenizers import CodingDNADelIns
from variation.schemas.token_response_schema import TokenType
from .tokenizer_base import TokenizerBase


class TestCodingDNADelInsTokenizer(TokenizerBase, unittest.TestCase):
    """A class for testing Coding DNA DelIns Tokenization."""

    def tokenizer_instance(self):
        """Return Coding DNA DelIns instance."""
        return CodingDNADelIns()

    def token_type(self):
        """Return DNA coding delins token type."""
        return TokenType.CODING_DNA_DELINS

    def fixture_name(self):
        """Return the fixture name for DNA coding delins."""
        return "coding_dna_delins"
