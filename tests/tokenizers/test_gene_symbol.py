"""A module for testing the Gene Symbol Tokenizer."""
import unittest

from gene.database.dynamodb import DynamoDbDatabase
from gene.query import QueryHandler as GeneQueryHandler

from variation.tokenizers import GeneSymbol
from .tokenizer_base import TokenizerBase


class TestGeneSymbolTokenizer(TokenizerBase, unittest.TestCase):
    """The Gene Symbol Tokenizer class."""

    def tokenizer_instance(self):
        """Return the Gene Symbol tokenizer instance."""
        return GeneSymbol(GeneQueryHandler(DynamoDbDatabase()))

    def token_type(self):
        """Return the Gene Symbol token type."""
        return "GeneSymbol"

    def fixture_name(self):
        """Return the fixture name for Gene Symbol."""
        return "gene"
