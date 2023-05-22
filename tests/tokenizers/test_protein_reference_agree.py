"""A module for testing Reference Agree tokenization."""
import unittest

from variation.tokenizers import ProteinReferenceAgree
from variation.schemas.token_response_schema import TokenType
from .tokenizer_base import TokenizerBase


class TestReferenceAgreeTokenizer(TokenizerBase, unittest.TestCase):
    """A class for testing Reference Agree Tokenization."""

    def tokenizer_instance(self):
        """Return Reference Agree instance."""
        return ProteinReferenceAgree()

    def token_type(self):
        """Return Reference Agree token type."""
        return TokenType.REFERENCE_AGREE

    def fixture_name(self):
        """Return the fixture name for Reference Agree."""
        return "reference_agree"
