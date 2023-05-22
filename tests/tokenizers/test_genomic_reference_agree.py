"""A module for testing Reference Agree tokenization on genomic reference sequence."""
import unittest

from variation.tokenizers import CdnaGenomicReferenceAgree
from variation.schemas.token_response_schema import TokenType
from .tokenizer_base import TokenizerBase


class TestGenomicReferenceAgreeTokenizer(TokenizerBase, unittest.TestCase):
    """A class for testing Reference Agree Tokenization on genomic reference
    sequence.
    """

    def tokenizer_instance(self):
        """Return CdnaGenomicReferenceAgree instance."""
        return CdnaGenomicReferenceAgree()

    def token_type(self):
        """Return Reference Agree token type."""
        return TokenType.GENOMIC_REFERENCE_AGREE

    def fixture_name(self):
        """Return the fixture name for Reference Agree."""
        return "genomic_reference_agree"
