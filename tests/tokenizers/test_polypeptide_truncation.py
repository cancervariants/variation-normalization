"""A module for testing Polypeptide Truncation tokenization."""
import unittest

from variation.tokenizers import PolypeptideTruncation
from variation.schemas.token_response_schema import TokenType
from .tokenizer_base import TokenizerBase


class TestPolypeptideTruncationTokenizer(TokenizerBase, unittest.TestCase):
    """A class for testing Polypeptide Truncation  Tokenization."""

    def tokenizer_instance(self):
        """Return Polypeptide Truncation instance."""
        return PolypeptideTruncation()

    def token_type(self):
        """Return Polypeptide Truncation token type."""
        return TokenType.POLYPEPTIDE_TRUNCATION

    def fixture_name(self):
        """Return the fixture name for Polypeptide Truncation."""
        return "polypeptide_truncation"
