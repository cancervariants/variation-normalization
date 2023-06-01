"""A module for testing Protein Substitution tokenization."""
import unittest

from variation.tokenizers import ProteinSubstitution
from variation.schemas.token_response_schema import TokenType
from .tokenizer_base import TokenizerBase


class TestProteinSubstitutionTokenizer(TokenizerBase, unittest.TestCase):
    """A class for testing Protein Substitution Tokenization."""

    def tokenizer_instance(self):
        """Return protein substitution instance."""
        return ProteinSubstitution()

    def token_type(self):
        """Return protein substitution token type."""
        return TokenType.PROTEIN_SUBSTITUTION

    def fixture_name(self):
        """Return the fixture name for protein substitution."""
        return "protein_substitution"
