"""A module for testing Protein Insertion Tokenization."""
import unittest

from variation.tokenizers import ProteinInsertion
from variation.schemas.token_response_schema import TokenType
from .tokenizer_base import TokenizerBase


class TestProteinInsertionTokenizer(TokenizerBase, unittest.TestCase):
    """A class for testing Protein Insertion Tokenization."""

    def tokenizer_instance(self):
        """Return Protein Insertion instance."""
        return ProteinInsertion()

    def token_type(self):
        """Return protein insertion token type."""
        return TokenType.PROTEIN_INSERTION

    def fixture_name(self):
        """Return the fixture name for protein insertion."""
        return "protein_insertion"
