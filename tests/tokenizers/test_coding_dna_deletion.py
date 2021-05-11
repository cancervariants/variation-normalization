"""A module for testing Coding DNA Deletion Tokenization."""
import unittest
from variant.tokenizers import CodingDNADeletion
from .tokenizer_base import TokenizerBase


class TestCodingDNADeletionTokenizer(TokenizerBase, unittest.TestCase):
    """A class for testing Coding DNA Deletion Tokenization."""

    def tokenizer_instance(self):
        """Return Coding DNA Deletion instance."""
        return CodingDNADeletion()

    def token_type(self):
        """Return Coding DNA deletion token type."""
        return 'CodingDNADeletion'

    def fixture_name(self):
        """Return the fixture name for Coding DNA Deletion."""
        return 'coding_dna_deletion'
