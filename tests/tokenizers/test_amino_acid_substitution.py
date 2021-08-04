"""A module for testing Amino Acid Substitution tokenization."""
import unittest
from variation.tokenizers import AminoAcidSubstitution
from .tokenizer_base import TokenizerBase
from variation.tokenizers.caches import AminoAcidCache


class TestAminoAcidSubstitutionTokenizer(TokenizerBase, unittest.TestCase):
    """A class for testing Amino Acid Substitution Tokenization."""

    def tokenizer_instance(self):
        """Return amino acid substitution instance."""
        return AminoAcidSubstitution(AminoAcidCache())

    def token_type(self):
        """Return amino acid substitution token type."""
        return 'AminoAcidSubstitution'

    def fixture_name(self):
        """Return the fixture name for amino acid substitution."""
        return 'amino_acid_substitution'
