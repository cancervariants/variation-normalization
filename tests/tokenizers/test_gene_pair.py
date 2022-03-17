"""Module for testing gene pair tokenization."""
import unittest

from variation.tokenizers import GenePair
from .tokenizer_base import TokenizerBase


class TestGenePairTokenizer(TokenizerBase, unittest.TestCase):
    """Gene Pair Tokenizer Test class."""

    def tokenizer_instance(self):
        """Return Gene Pair Tokenizer instance."""
        return GenePair()

    def token_type(self):
        """Return Gene Pair token type."""
        return "GenePair"

    def fixture_name(self):
        """Return Gene Pair fixture name."""
        return "gene_pair"
