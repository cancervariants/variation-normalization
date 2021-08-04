"""A module for testing the Gene Pair Tokenizer."""
import unittest
from variation.tokenizers import GeneSymbol
from variation.tokenizers.caches import GeneSymbolCache
from .tokenizer_base import TokenizerBase


class TestGenePairTokenizer(TokenizerBase, unittest.TestCase):
    """The Gene Pair Tokenizer class."""

    def tokenizer_instance(self):
        """Return the Gene Pair tokenizer instance."""
        gene_cache = GeneSymbolCache()
        return GeneSymbol(gene_cache)

    def token_type(self):
        """Return the Gene Pair token type."""
        return 'GeneSymbol'

    def fixture_name(self):
        """Return the fixture name for Gene Pair."""
        return 'gene'
