"""A module for testing DNA Coding Substitution tokenization."""
import unittest
from variant.tokenizers import CodingDNASubstitution
from .tokenizer_base import TokenizerBase


class TestCodingDNASubstitutionTokenizer(TokenizerBase, unittest.TestCase):
    """A class for testing Coding DNA Substitution Tokenization."""

    def tokenizer_instance(self):
        """Return Coding DNA Substitution instance."""
        return CodingDNASubstitution()

    def token_type(self):
        """Return DNA coding substitution token type."""
        return 'CodingDNASubstitution'

    def fixture_name(self):
        """Return the fixture name for DNA coding substitution."""
        return 'coding_dna_substitution'
