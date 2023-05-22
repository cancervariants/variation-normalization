"""A module for testing Genomic Duplication Ambiguous Tokenization."""
import unittest

from variation.tokenizers import GenomicDuplication
from variation.schemas.token_response_schema import TokenType
from .tokenizer_base import TokenizerBase


class TestGenomicDuplicationRangeTokenizer(TokenizerBase, unittest.TestCase):
    """A class for testing Genomic Duplication Ambiguous Tokenization."""

    def tokenizer_instance(self):
        """Return Genomic Duplication instance."""
        return GenomicDuplication()

    def token_type(self):
        """Return genomic duplication ambiguous token type."""
        return TokenType.GENOMIC_DUPLICATION_AMBIGUOUS

    def fixture_name(self):
        """Return the fixture name for Genomic Duplication Ambiguous."""
        return "genomic_duplication_ambiguous"
