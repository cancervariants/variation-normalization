"""A module for testing classifier classes."""
import yaml
from variation.tokenizers import Tokenize
from tests import PROJECT_ROOT


class ClassifierBase:
    """The classifier base class."""

    def setUp(self):
        """Set up the test cases."""
        with open(f'{PROJECT_ROOT}/tests/fixtures/classifiers.yml') as stream:
            self.all_fixtures = yaml.safe_load(stream)
        self.fixtures = self.all_fixtures.get(
            self.fixture_name(),
            {'should_match': [], 'should_not_match': []}
        )
        self.classifier = self.classifier_instance()
        self.tokenizer = Tokenize()

    def classifier_instance(self):
        """Check that the classifier_instance method is implemented."""
        raise NotImplementedError()

    def fixture_name(self):
        """Check that the fixture_name method is implemented."""
        raise NotImplementedError()

    def test_matches(self):
        """Test that classifier matches correctly."""
        for x in self.fixtures['should_match']:
            tokens = self.tokenizer.perform(x['query'], [])
            classification = self.classifier.match(tokens)
            self.assertIsNotNone(classification, msg=x)
            self.assertEqual(x['confidence'],
                             str(classification.confidence),
                             msg=x)

    def test_not_matches(self):
        """Test that classifier matches correctly."""
        for x in self.fixtures['should_not_match']:
            tokens = self.tokenizer.perform(x['query'], [])
            classification = self.classifier.match(tokens)
            self.assertIsNone(classification, msg=x)
