"""A module for testing Protein Insertion Tokenization."""
import unittest

from variation.tokenizers import ProteinInsertion
from variation.tokenizers.caches import AminoAcidCache, NucleotideCache
from .tokenizer_base import TokenizerBase


class TestProteinInsertionTokenizer(TokenizerBase, unittest.TestCase):
    """A class for testing Protein Insertion Tokenization."""

    def tokenizer_instance(self):
        """Return Protein Insertion instance."""
        return ProteinInsertion(AminoAcidCache(), NucleotideCache())

    def token_type(self):
        """Return protein insertion token type."""
        return "ProteinInsertion"

    def fixture_name(self):
        """Return the fixture name for protein insertion."""
        return "protein_insertion"
