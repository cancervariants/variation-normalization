"""A module for testing HGVS tokenization."""
import unittest

from variation.tokenizers import HGVS
from variation.schemas.token_response_schema import TokenType
from .tokenizer_base import TokenizerBase


class TestHGVSTokenizer(TokenizerBase, unittest.TestCase):
    """A class for testing HGVS Tokenization."""

    def tokenizer_instance(self):
        """Return HGVS instance."""
        return HGVS()

    def token_type(self):
        """Return HGVS token type."""
        return TokenType.HGVS

    def fixture_name(self):
        """Return fixture name for HGVS."""
        return "hgvs"
