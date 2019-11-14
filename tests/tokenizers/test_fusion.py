import unittest

from varlexapp.tokenizers import Fusion
from .tokenizer_base import TokenizerBase

class TestFusionTokenizer(TokenizerBase, unittest.TestCase):
    def tokenizer_instance(self):
        return Fusion()

    def token_type(self):
        return 'Fusion'

    def fixture_name(self):
        return 'fusion'
