"""A module for testing Cdna Substitution tokenization."""
import unittest

from variation.tokenizers import CdnaSubstitution
from variation.schemas.token_response_schema import TokenType
from .tokenizer_base import TokenizerBase


class TestCdnaSubstitutionTokenizer(TokenizerBase, unittest.TestCase):
    """A class for testing Cdna Substitution Tokenization."""

    def tokenizer_instance(self):
        """Return Cdna Substitution instance."""
        return CdnaSubstitution()

    def token_type(self):
        """Return cDNA substitution token type."""
        return TokenType.CDNA_SUBSTITUTION

    def fixture_name(self):
        """Return the fixture name for cDNA substitution."""
        return "cdna_substitution"
