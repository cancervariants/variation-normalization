from varlexapp.tokenizers import Fusion
from .tokenizer_test import TokenizerTest

class TestFusionTokenizer(TokenizerTest):
    def tokenizer_instance(self):
        return Fusion()

    def token_type(self):
        return 'Fusion'

    def fixture_name(self):
        return 'fusion'
