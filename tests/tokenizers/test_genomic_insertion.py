"""A module for testing Genomic Insertion Tokenization."""
import unittest

from variation.tokenizers import GenomicInsertion
from variation.schemas.token_response_schema import TokenType
from .tokenizer_base import TokenizerBase


class TestGenomicInsertionTokenizer(TokenizerBase, unittest.TestCase):
    """A class for testing Genomic Insertion Tokenization."""

    def tokenizer_instance(self):
        """Return Genomic Insertion instance."""
        return GenomicInsertion()

    def token_type(self):
        """Return genomic insertion token type."""
        return TokenType.GENOMIC_INSERTION

    def fixture_name(self):
        """Return the fixture name for Genomic Insertion."""
        return "genomic_insertion"
