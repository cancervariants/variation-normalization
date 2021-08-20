"""A module for testing Genomic Copy Number Loss tokenization."""
import unittest
from variation.tokenizers import GenomicCopyNumberLoss
from .tokenizer_base import TokenizerBase
from variation.tokenizers.caches import AminoAcidCache, NucleotideCache


class TestGenomicCopyNumberLossTokenizer(TokenizerBase, unittest.TestCase):
    """A class for testing Genomic Copy Number Loss Tokenization."""

    def tokenizer_instance(self):
        """Return Genomic Copy Number Loss instance."""
        return GenomicCopyNumberLoss(AminoAcidCache(), NucleotideCache())

    def token_type(self):
        """Return genomic copy number loss token type."""
        return 'GenomicCopyNumberLoss'

    def fixture_name(self):
        """Return the fixture name for genomic delins."""
        return 'genomic_copy_number_loss'
