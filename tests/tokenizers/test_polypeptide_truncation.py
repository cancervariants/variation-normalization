"""A module for testing Polypeptide Truncation tokenization."""
import unittest

from variation.tokenizers.caches import AminoAcidCache
from variation.tokenizers import PolypeptideTruncation
from .tokenizer_base import TokenizerBase


class TestPolypeptideTruncationTokenizer(TokenizerBase, unittest.TestCase):
    """A class for testing Polypeptide Truncation  Tokenization."""

    def tokenizer_instance(self):
        """Return Polypeptide Truncation instance."""
        return PolypeptideTruncation(AminoAcidCache())

    def token_type(self):
        """Return Polypeptide Truncation token type."""
        return "PolypeptideTruncation"

    def fixture_name(self):
        """Return the fixture name for Polypeptide Truncation."""
        return "polypeptide_truncation"
