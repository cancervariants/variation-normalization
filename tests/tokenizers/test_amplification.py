"""Module for testing Amplification Tokenization"""
import unittest

from variation.schemas.token_response_schema import TokenType
from variation.tokenizers import FreeTextCategorical
from .tokenizer_base import TokenizerBase


class TestAmplificationTokenizer(TokenizerBase, unittest.TestCase):
    """A class for testing Amplification tokenization"""

    def tokenizer_instance(self):
        """Return the FreeTextCategorical class instance"""
        return FreeTextCategorical()

    def token_type(self):
        """Return the token type for amplification"""
        return TokenType.AMPLIFICATION

    def fixture_name(self):
        """Return the fixture name for Amplification"""
        return "amplification"
