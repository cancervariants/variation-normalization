"""A module for testing Coding DNA Silent Mutation tokenization."""
import unittest

from variation.tokenizers import CodingDNASilentMutation
from .tokenizer_base import TokenizerBase


class TestCodingDNASilentMutationTokenizer(TokenizerBase, unittest.TestCase):
    """A class for testing Coding DNA Silent Mutation Tokenization."""

    def tokenizer_instance(self):
        """Return Coding DNA Silent Mutation instance."""
        return CodingDNASilentMutation()

    def token_type(self):
        """Return coding DNA silent mutation token type."""
        return "CodingDNASilentMutation"

    def fixture_name(self):
        """Return the fixture name for coding DNA silent mutation."""
        return "coding_dna_silent_mutation"
