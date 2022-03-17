"""A module for testing Genomic Silent Mutation tokenization."""
import unittest

from variation.tokenizers import GenomicSilentMutation
from .tokenizer_base import TokenizerBase


class TestGenomicSilentMutationTokenizer(TokenizerBase, unittest.TestCase):
    """A class for testing Genomic Silent Mutation Tokenization."""

    def tokenizer_instance(self):
        """Return Genomic Silent Mutation instance."""
        return GenomicSilentMutation()

    def token_type(self):
        """Return genomic silent mutation token type."""
        return "GenomicSilentMutation"

    def fixture_name(self):
        """Return the fixture name for coding  DNA silent mutation."""
        return "genomic_silent_mutation"
