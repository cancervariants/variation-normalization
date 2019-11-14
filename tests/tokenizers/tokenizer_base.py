import unittest
import yaml

class TokenizerBase(object):

    def setUp(self):
        with open('tests/fixtures/tokenizers.yml') as stream:
            self.all_fixtures = yaml.safe_load(stream)
        self.fixtures = self.all_fixtures.get(
                self.fixture_name(),
                {'should_match': [], 'should_not_match': []}
        )
        self.tokenizer = self.tokenizer_instance()

    def tokenizer_instance(self):
        raise NotImplementedError()

    def token_type(self):
        raise NotImplementedError()

    def fixture_name(self):
        raise NotImplementedError()

    def test_matches(self):
        for x in self.fixtures['should_match']:
            res = self.tokenizer.match(x['token'])
            self.assertIsNotNone(res)
            self.assertEqual(self.token_type(), res.token_type)

    def test_not_matches(self):
        for x in self.fixtures['should_not_match']:
            res = self.tokenizer.match(x['token'])
            self.assertIsNone(res)
