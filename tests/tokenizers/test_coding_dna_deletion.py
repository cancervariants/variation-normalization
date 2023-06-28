"""A module for testing Coding DNA Deletion Tokenization."""
import unittest

from variation.tokenizers import CdnaDeletion
from variation.schemas.token_response_schema import TokenType
from .tokenizer_base import TokenizerBase


class TestCodingDNADeletionTokenizer(TokenizerBase, unittest.TestCase):
    """A class for testing Coding DNA Deletion Tokenization."""

    def tokenizer_instance(self):
        """Return Coding DNA Deletion instance."""
        return CdnaDeletion()

    def token_type(self):
        """Return Coding DNA deletion token type."""
        return TokenType.CODING_DNA_DELETION

    def fixture_name(self):
        """Return the fixture name for Coding DNA Deletion."""
        return "coding_dna_deletion"
