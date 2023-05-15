"""A module for testing Genomic Substitution tokenization."""
import unittest

from variation.tokenizers import GenomicSubstitution
from variation.schemas.token_response_schema import TokenType
from .tokenizer_base import TokenizerBase


class TestGenomicSubstitutionTokenizer(TokenizerBase, unittest.TestCase):
    """A class for testing Genomic Substitution Tokenization."""

    def tokenizer_instance(self):
        """Return Genomic Substitution instance."""
        return GenomicSubstitution()

    def token_type(self):
        """Return genomic substitution token type."""
        return TokenType.GENOMIC_SUBSTITUTION

    def fixture_name(self):
        """Return the fixture name for genomic substitution."""
        return "genomic_substitution"
