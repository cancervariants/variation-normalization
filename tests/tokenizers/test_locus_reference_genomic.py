"""A module for testing LRG tokenization."""
import unittest

from variation.tokenizers import LocusReferenceGenomic
from .tokenizer_base import TokenizerBase


class TestLocusReferenceGenomicTokenizer(TokenizerBase, unittest.TestCase):
    """A class for testing LRG Tokenization."""

    def tokenizer_instance(self):
        """Return LRG instance."""
        return LocusReferenceGenomic()

    def token_type(self):
        """Return LRG token type."""
        return "LocusReferenceGenomic"

    def fixture_name(self):
        """Return fixture name for LRG."""
        return "locus_reference_genomic"
