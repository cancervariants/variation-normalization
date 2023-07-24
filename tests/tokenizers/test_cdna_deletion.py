"""A module for testing Cdna Deletion Tokenization."""
import unittest

from variation.tokenizers import CdnaDeletion
from variation.schemas.token_response_schema import TokenType
from .tokenizer_base import TokenizerBase


class TestCdnaDeletionTokenizer(TokenizerBase, unittest.TestCase):
    """A class for testing Cdna Deletion Tokenization."""

    def tokenizer_instance(self):
        """Return Cdna Deletion instance."""
        return CdnaDeletion()

    def token_type(self):
        """Return cdna deletion token type."""
        return TokenType.CDNA_DELETION

    def fixture_name(self):
        """Return the fixture name for Cdna Deletion."""
        return "cdna_deletion"
