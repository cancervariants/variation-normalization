"""A module for testing Amino Acid Deletion Tokenization."""
import unittest
from variation.tokenizers import AminoAcidDeletion
from .tokenizer_base import TokenizerBase
from variation.tokenizers.caches import AminoAcidCache, NucleotideCache


class TestAminoAcidDeletionTokenizer(TokenizerBase, unittest.TestCase):
    """A class for testing Amino Acid Deletion Tokenization."""

    def tokenizer_instance(self):
        """Return AminoAcid Deletion instance."""
        return AminoAcidDeletion(AminoAcidCache(), NucleotideCache())

    def token_type(self):
        """Return amino acid deletion token type."""
        return 'AminoAcidDeletion'

    def fixture_name(self):
        """Return the fixture name for amino acid deletion."""
        return 'amino_acid_deletion'
