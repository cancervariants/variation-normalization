"""Module for testing gene pair tokenization."""
import unittest
from variant.tokenizers import GenePair
from variant.tokenizers.caches import GeneSymbolCache
from .tokenizer_base import TokenizerBase


class TestGenePairTokenizer(TokenizerBase, unittest.TestCase):
    """Gene Pair Tokenizer Test class."""

    # TODO: don't hardcode this, inject with config
    def tokenizer_instance(self):
        """Return Gene Pair Tokenizer instance."""
        gene_cache = GeneSymbolCache('variant/data/gene_symbols.txt')
        return GenePair(gene_cache)

    def token_type(self):
        """Return Gene Pair token type."""
        return 'GenePair'

    def fixture_name(self):
        """Return Gene Pair fixture name."""
        return 'gene_pair'
