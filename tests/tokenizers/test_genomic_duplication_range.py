"""A module for testing Genomic Duplication Range Tokenization."""
import unittest

from variation.tokenizers import GenomicDuplication
from .tokenizer_base import TokenizerBase


class TestGenomicDuplicationRangeTokenizer(TokenizerBase, unittest.TestCase):
    """A class for testing Genomic Duplication Tokenization."""

    def tokenizer_instance(self):
        """Return Genomic Duplication instance."""
        return GenomicDuplication()

    def token_type(self):
        """Return genomic duplication token type."""
        return "GenomicDuplicationRange"

    def fixture_name(self):
        """Return the fixture name for Genomic Duplication."""
        return "genomic_duplication_range"
