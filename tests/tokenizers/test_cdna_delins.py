"""A module for testing Cdna DelIns tokenization."""
import unittest

from variation.tokenizers import CdnaDelIns
from variation.schemas.token_response_schema import TokenType
from .tokenizer_base import TokenizerBase


class TestCdnaDelInsTokenizer(TokenizerBase, unittest.TestCase):
    """A class for testing Cdna DelIns Tokenization."""

    def tokenizer_instance(self):
        """Return Cdna DelIns instance."""
        return CdnaDelIns()

    def token_type(self):
        """Return cdna delins token type."""
        return TokenType.CDNA_DELINS

    def fixture_name(self):
        """Return the fixture name for cdna delins."""
        return "cdna_delins"
