"""A module for testing Coding DNA DelIns tokenization."""
import unittest

from variation.tokenizers.caches import AminoAcidCache, NucleotideCache
from variation.tokenizers import CodingDNADelIns
from .tokenizer_base import TokenizerBase


class TestCodingDNADelInsTokenizer(TokenizerBase, unittest.TestCase):
    """A class for testing Coding DNA DelIns Tokenization."""

    def tokenizer_instance(self):
        """Return Coding DNA DelIns instance."""
        return CodingDNADelIns(AminoAcidCache(), NucleotideCache())

    def token_type(self):
        """Return DNA coding delins token type."""
        return "CodingDNADelIns"

    def fixture_name(self):
        """Return the fixture name for DNA coding delins."""
        return "coding_dna_delins"
