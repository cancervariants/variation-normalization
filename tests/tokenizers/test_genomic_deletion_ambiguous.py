"""A module for testing Genomic Deletion Ambiguous tokenization."""
import unittest

from variation.tokenizers import GenomicDeletion
from variation.schemas.token_response_schema import TokenType
from .tokenizer_base import TokenizerBase


class TestGenomicDeletionRangeTokenizer(TokenizerBase, unittest.TestCase):
    """A class for testing Genomic Deletion Ambiguous Tokenization."""

    def tokenizer_instance(self):
        """Return Genomic Deletion Ambiguous instance."""
        return GenomicDeletion()

    def token_type(self):
        """Return genomic deletion ambigous token type."""
        return TokenType.GENOMIC_DELETION_AMBIGUOUS

    def fixture_name(self):
        """Return the fixture name for genomic deletion ambiguous."""
        return "genomic_deletion_ambiguous"
