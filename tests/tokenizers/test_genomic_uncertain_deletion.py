"""A module for testing Genomic Uncertain Deletion tokenization."""
import unittest

from variation.tokenizers import GenomicUncertainDeletion
from variation.schemas.token_response_schema import TokenType
from .tokenizer_base import TokenizerBase


class TestGenomicUncertainDeletionTokenizer(TokenizerBase, unittest.TestCase):
    """A class for testing Genomic Uncertain Deletion Tokenization."""

    def tokenizer_instance(self):
        """Return Genomic Uncertain Deletion instance."""
        return GenomicUncertainDeletion()

    def token_type(self):
        """Return genomic uncertain deletion token type."""
        return TokenType.GENOMIC_UNCERTAIN_DELETION

    def fixture_name(self):
        """Return the fixture name for genomic uncertain deletion."""
        return "genomic_uncertain_deletion"
