import unittest
import yaml

from varlexapp.tokenizers import Tokenize

class ClassifierBase(object):

    def setUp(self):
        with open('tests/fixtures/classifiers.yml') as stream:
            self.all_fixtures = yaml.safe_load(stream)
        self.fixtures = self.all_fixtures.get(
                self.fixture_name(),
                {'should_match': [], 'should_not_match': []}
        )
        self.classifier = self.classifier_instance()
        self.tokenizer = Tokenize('varlexapp/data/gene_symbols.txt')

    def classifier_instance(self):
        raise NotImplementedError()

    def fixture_name(self):
        raise NotImplementedError()

    def test_matches(self):
        for x in self.fixtures['should_match']:
            tokens = self.tokenizer.perform(x['query'])
            classification = self.classifier.match(tokens)
            self.assertIsNotNone(classification)
            self.assertEqual(str(classification.confidence), x['confidence'])

    def test_not_matches(self):
        for x in self.fixtures['should_not_match']:
            tokens = self.tokenizer.perform(x['query'])
            classification = self.classifier.match(tokens)
            self.assertIsNone(classification)
