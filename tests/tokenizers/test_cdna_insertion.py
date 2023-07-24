"""A module for testing Cdna Insertion Tokenization."""
import unittest

from variation.tokenizers import CdnaInsertion
from variation.schemas.token_response_schema import TokenType
from .tokenizer_base import TokenizerBase


class TestCdnaInsertionTokenizer(TokenizerBase, unittest.TestCase):
    """A class for testing Cdna Insertion Tokenization."""

    def tokenizer_instance(self):
        """Return Cdna Insertion instance."""
        return CdnaInsertion()

    def token_type(self):
        """Return Cdna insertion token type."""
        return TokenType.CDNA_INSERTION

    def fixture_name(self):
        """Return the fixture name for Cdna Insertion."""
        return "cdna_insertion"
