"""Module for testing gene pair tokenization."""
import unittest
from variant.tokenizers import GenePair
from variant.tokenizers.caches import GeneSymbolCache
from .tokenizer_base import TokenizerBase


class TestGenePairTokenizer(TokenizerBase, unittest.TestCase):
    """Gene Pair Tokenizer Test class."""

    def tokenizer_instance(self):
        """Return Gene Pair Tokenizer instance."""
        return GenePair(GeneSymbolCache())

    def token_type(self):
        """Return Gene Pair token type."""
        return 'GenePair'

    def fixture_name(self):
        """Return Gene Pair fixture name."""
        return 'gene_pair'
