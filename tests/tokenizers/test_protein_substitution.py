"""A module for testing Protein Substitution tokenization."""
import unittest

from variation.tokenizers.caches import AminoAcidCache
from variation.tokenizers import ProteinSubstitution
from .tokenizer_base import TokenizerBase


class TestProteinSubstitutionTokenizer(TokenizerBase, unittest.TestCase):
    """A class for testing Protein Substitution Tokenization."""

    def tokenizer_instance(self):
        """Return protein substitution instance."""
        return ProteinSubstitution(AminoAcidCache())

    def token_type(self):
        """Return protein substitution token type."""
        return "ProteinSubstitution"

    def fixture_name(self):
        """Return the fixture name for protein substitution."""
        return "protein_substitution"
