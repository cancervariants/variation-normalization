"""A module for testing the Gene Pair Tokenizer."""
import unittest

from gene.query import QueryHandler as GeneQueryHandler

from variation.tokenizers import GeneSymbol
from .tokenizer_base import TokenizerBase


class TestGenePairTokenizer(TokenizerBase, unittest.TestCase):
    """The Gene Pair Tokenizer class."""

    def tokenizer_instance(self):
        """Return the Gene Pair tokenizer instance."""
        return GeneSymbol(GeneQueryHandler())

    def token_type(self):
        """Return the Gene Pair token type."""
        return "GeneSymbol"

    def fixture_name(self):
        """Return the fixture name for Gene Pair."""
        return "gene"
