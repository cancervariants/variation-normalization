import unittest

from varlexapp.tokenizers import ProteinSubstitution
from .tokenizer_base import TokenizerBase
from varlexapp.tokenizers.caches import AminoAcidCache

class TestProteinSubstitutionTokenizer(TokenizerBase, unittest.TestCase):
    def tokenizer_instance(self):
        return ProteinSubstitution(AminoAcidCache())

    def token_type(self):
        return 'ProteinSubstitution'

    def fixture_name(self):
        return 'protein_substitution'
