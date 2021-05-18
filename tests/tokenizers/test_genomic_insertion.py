"""A module for testing Genomic Insertion Tokenization."""
import unittest
from variant.tokenizers import GenomicInsertion
from .tokenizer_base import TokenizerBase
from variant.tokenizers.caches import AminoAcidCache, NucleotideCache


class TestGenomicInsertionTokenizer(TokenizerBase, unittest.TestCase):
    """A class for testing Genomic Insertion Tokenization."""

    def tokenizer_instance(self):
        """Return Genomic Insertion instance."""
        return GenomicInsertion(AminoAcidCache(), NucleotideCache())

    def token_type(self):
        """Return genomic insertion token type."""
        return 'GenomicInsertion'

    def fixture_name(self):
        """Return the fixture name for Genomic Insertion."""
        return 'genomic_insertion'
