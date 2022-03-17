"""A module for testing Protein DelIns tokenization."""
import unittest

from variation.tokenizers.caches import AminoAcidCache, NucleotideCache
from variation.tokenizers import ProteinDelIns
from .tokenizer_base import TokenizerBase


class TestProteinDelInsTokenizer(TokenizerBase, unittest.TestCase):
    """A class for testing Protein DelIns Tokenization."""

    def tokenizer_instance(self):
        """Return Protein DelIns instance."""
        return ProteinDelIns(AminoAcidCache(), NucleotideCache())

    def token_type(self):
        """Return protein delins token type."""
        return "ProteinDelIns"

    def fixture_name(self):
        """Return the fixture name for protein delins."""
        return "protein_delins"
