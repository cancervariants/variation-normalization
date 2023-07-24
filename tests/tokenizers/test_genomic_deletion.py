"""A module for testing Genomic Deletion Tokenization."""
import unittest

from variation.tokenizers import GenomicDeletion
from variation.schemas.token_response_schema import TokenType
from .tokenizer_base import TokenizerBase


class TestGenomicDeletionTokenizer(TokenizerBase, unittest.TestCase):
    """A class for testing Cdna Deletion Tokenization."""

    def tokenizer_instance(self):
        """Return Genomic Deletion instance."""
        return GenomicDeletion()

    def token_type(self):
        """Return genomic deletion token type."""
        return TokenType.GENOMIC_DELETION

    def fixture_name(self):
        """Return the fixture name for Genomic Deletion."""
        return "genomic_deletion"
