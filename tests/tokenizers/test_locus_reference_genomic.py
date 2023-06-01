"""A module for testing LRG tokenization."""
import unittest

from variation.tokenizers import LocusReferenceGenomic
from variation.schemas.token_response_schema import TokenType
from .tokenizer_base import TokenizerBase


class TestLocusReferenceGenomicTokenizer(TokenizerBase, unittest.TestCase):
    """A class for testing LRG Tokenization."""

    def tokenizer_instance(self):
        """Return LRG instance."""
        return LocusReferenceGenomic()

    def token_type(self):
        """Return LRG token type."""
        return TokenType.LOCUS_REFERENCE_GENOMIC

    def fixture_name(self):
        """Return fixture name for LRG."""
        return "locus_reference_genomic"
