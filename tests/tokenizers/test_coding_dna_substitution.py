"""A module for testing DNA Coding Substitution tokenization."""
import unittest
from variant.tokenizers import DNASubstitution
from .tokenizer_base import TokenizerBase


class TestDNACodingSubstitutionTokenizer(TokenizerBase, unittest.TestCase):
    """A class for testing DNA Coding Substitution Tokenization."""

    def tokenizer_instance(self):
        """Return amino acid substitution instance."""
        return DNASubstitution()

    def token_type(self):
        """Return DNA coding substitution token type."""
        return 'CodingDNASubstitution'

    def fixture_name(self):
        """Return the fixture name for DNA coding substitution."""
        return 'coding_dna_substitution'
