"""A module for testing Protein Deletion Tokenization."""
import unittest

from variation.tokenizers import ProteinDeletion
from .tokenizer_base import TokenizerBase


class TestProteinDeletionTokenizer(TokenizerBase, unittest.TestCase):
    """A class for testing Protein Deletion Tokenization."""

    def tokenizer_instance(self):
        """Return Protein Deletion instance."""
        return ProteinDeletion()

    def token_type(self):
        """Return protein deletion token type."""
        return "ProteinDeletion"

    def fixture_name(self):
        """Return the fixture name for protein deletion."""
        return "protein_deletion"
