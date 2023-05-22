"""A module for testing protein stop gain tokenization."""
import unittest

from variation.tokenizers import ProteinSubstitution
from variation.schemas.token_response_schema import TokenType
from .tokenizer_base import TokenizerBase


class TestProteinStopGainTokenizer(TokenizerBase, unittest.TestCase):
    """A class for testing Protein Stop Gain Tokenization."""

    def tokenizer_instance(self):
        """Return Protein Substitution instance."""
        return ProteinSubstitution()

    def token_type(self):
        """Return Polypeptide Truncation token type."""
        return TokenType.PROTEIN_STOP_GAIN

    def fixture_name(self):
        """Return the fixture name for Polypeptide Truncation."""
        return "protein_stop_gain"
