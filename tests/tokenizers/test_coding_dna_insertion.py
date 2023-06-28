"""A module for testing Coding DNA Insertion Tokenization."""
import unittest

from variation.tokenizers import CdnaInsertion
from variation.schemas.token_response_schema import TokenType
from .tokenizer_base import TokenizerBase


class TestCodingDNAInsertionTokenizer(TokenizerBase, unittest.TestCase):
    """A class for testing Coding DNA Insertion Tokenization."""

    def tokenizer_instance(self):
        """Return Coding DNA Insertion instance."""
        return CdnaInsertion()

    def token_type(self):
        """Return Coding DNA insertion token type."""
        return TokenType.CODING_DNA_INSERTION

    def fixture_name(self):
        """Return the fixture name for Coding DNA Insertion."""
        return "coding_dna_insertion"
