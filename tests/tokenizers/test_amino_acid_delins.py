"""A module for testing Amino Acid DelIns tokenization."""
import unittest
from variation.tokenizers import AminoAcidDelIns
from .tokenizer_base import TokenizerBase
from variation.tokenizers.caches import AminoAcidCache, NucleotideCache


class TestAminoAcidDelInsTokenizer(TokenizerBase, unittest.TestCase):
    """A class for testing Amino Acid DelIns Tokenization."""

    def tokenizer_instance(self):
        """Return AminoAcid DelIns instance."""
        return AminoAcidDelIns(AminoAcidCache(), NucleotideCache())

    def token_type(self):
        """Return amino acid delins token type."""
        return 'AminoAcidDelIns'

    def fixture_name(self):
        """Return the fixture name for amino acid delins."""
        return 'amino_acid_delins'
