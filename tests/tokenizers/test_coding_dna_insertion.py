"""A module for testing Coding DNA Insertion Tokenization."""
import unittest

from variation.tokenizers.caches import AminoAcidCache, NucleotideCache
from variation.tokenizers import CodingDNAInsertion
from .tokenizer_base import TokenizerBase


class TestCodingDNAInsertionTokenizer(TokenizerBase, unittest.TestCase):
    """A class for testing Coding DNA Insertion Tokenization."""

    def tokenizer_instance(self):
        """Return Coding DNA Insertion instance."""
        return CodingDNAInsertion(AminoAcidCache(), NucleotideCache())

    def token_type(self):
        """Return Coding DNA insertion token type."""
        return "CodingDNAInsertion"

    def fixture_name(self):
        """Return the fixture name for Coding DNA Insertion."""
        return "coding_dna_insertion"
