"""A module for testing Genomic Deletion Ranges tokenization."""
import unittest

from variation.tokenizers import GenomicDeletionRange
from .tokenizer_base import TokenizerBase


class TestGenomicDeletionRangeTokenizer(TokenizerBase, unittest.TestCase):
    """A class for testing Genomic Deletion Range Tokenization."""

    def tokenizer_instance(self):
        """Return Genomic Deletion Range instance."""
        return GenomicDeletionRange()

    def token_type(self):
        """Return genomic deletion range token type."""
        return "GenomicDeletionRange"

    def fixture_name(self):
        """Return the fixture name for genomic deletion range."""
        return "genomic_deletion_range"
