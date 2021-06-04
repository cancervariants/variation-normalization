"""Module for testing fusion tokenization."""
import unittest
from variation.tokenizers import Fusion
from .tokenizer_base import TokenizerBase


class TestFusionTokenizer(TokenizerBase, unittest.TestCase):
    """Fusion tokenizer test class."""

    def tokenizer_instance(self):
        """Return fusion tokenizer instance."""
        return Fusion()

    def token_type(self):
        """Return fusion token type."""
        return 'Fusion'

    def fixture_name(self):
        """Return fusion fixture name."""
        return 'fusion'
