"""A module for testing Amino Acid Insertion Tokenization."""
import unittest
from variant.tokenizers import AminoAcidInsertion
from .tokenizer_base import TokenizerBase
from variant.tokenizers.caches import AminoAcidCache


class TestAminoAcidInsertionTokenizer(TokenizerBase, unittest.TestCase):
    """A class for testing Amino Acid Insertion Tokenization."""

    def tokenizer_instance(self):
        """Return Amino Acid Insertion instance."""
        return AminoAcidInsertion(AminoAcidCache())

    def token_type(self):
        """Return amino acid insertion token type."""
        return 'AminoAcidInsertion'

    def fixture_name(self):
        """Return the fixture name for amino acid insertion."""
        return 'amino_acid_insertion'
