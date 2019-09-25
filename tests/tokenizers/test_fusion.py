import unittest

from varlexapp.tokenizers import Fusion

class TestFusionTokenizer(unittest.TestCase):
    def setUp(self):
        self.tokenizer = Fusion()
        self.should_match = [
            'Fusion',
            'fusions',
            'FUSION',
            'fuSIons'
        ]

        self.should_not_match = [
            'afusion',
            'fused'
        ]

    def test_matches(self):
        for x in self.should_match:
            res = self.tokenizer.match(x)
            self.assertIsNotNone(res)
            self.assertEqual('Fusion', res.token_type)

        for x in self.should_not_match:
            res = self.tokenizer.match(x)
            self.assertIsNone(res)



