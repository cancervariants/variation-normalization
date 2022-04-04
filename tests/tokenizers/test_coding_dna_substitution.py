"""A module for testing Coding DNA Substitution tokenization."""
import unittest

from variation.tokenizers import CodingDNASubstitution
from .tokenizer_base import TokenizerBase


class TestCodingDNASubstitutionTokenizer(TokenizerBase, unittest.TestCase):
    """A class for testing Coding DNA Substitution Tokenization."""

    def tokenizer_instance(self):
        """Return Coding DNA Substitution instance."""
        return CodingDNASubstitution()

    def token_type(self):
        """Return coding DNA substitution token type."""
        return "CodingDNASubstitution"

    def fixture_name(self):
        """Return the fixture name for coding DNA substitution."""
        return "coding_dna_substitution"
