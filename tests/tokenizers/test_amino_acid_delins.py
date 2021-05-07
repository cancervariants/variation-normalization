"""A module for testing Amino Acid DelIns tokenization."""
import unittest
from variant.tokenizers import AminoAcidDelIns
from .tokenizer_base import TokenizerBase
from variant.tokenizers.caches import AminoAcidCache


class TestAminoAcidDelInsTokenizer(TokenizerBase, unittest.TestCase):
    """A class for testing Amino Acid DelIns Tokenization."""

    def tokenizer_instance(self):
        """Return AminoAcid DelIns instance."""
        return AminoAcidDelIns(AminoAcidCache())

    def token_type(self):
        """Return amino acid delins token type."""
        return 'AminoAcidDelIns'

    def fixture_name(self):
        """Return the fixture name for amino acid delins."""
        return 'amino_acid_delins'
