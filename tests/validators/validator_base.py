"""A module for testing validator classes."""
import yaml
from varlexapp import PROJECT_ROOT
from varlexapp.tokenizers import Tokenize


class ValidatorBase(object):
    """The validator base class."""

    def setUp(self):
        """Set up the test cases."""
        with open(f'{PROJECT_ROOT}/tests/fixtures/validators.yml') as stream:
            self.all_fixtures = yaml.safe_load(stream)
        self.fixtures = self.all_fixtures.get(
                self.fixture_name(),
                {'should_match': [], 'should_not_match': []}
        )
        self.tokenizer = Tokenize('varlexapp/data/gene_symbols.txt')
        self.classifier = self.classifier_instance()
        self.validator = self.validator_instance()

    def classifier_instance(self):
        """Check that the classifier_instance method is implemented."""
        raise NotImplementedError()

    def validator_instance(self):
        """Check that the validator_instance method is implemented."""
        raise NotImplementedError()

    def fixture_name(self):
        """Check that the fixture_name method is implemented."""
        raise NotImplementedError()

    def test_matches(self):
        """Test that validator matches correctly."""
        for x in self.fixtures['should_match']:
            tokens = self.tokenizer.perform(x['token'])
            classification = self.classifier.match(tokens)
            res = self.validator.validate(classification)
            self.assertIsNotNone(res, msg=x)

    def test_not_matches(self):
        """Test that validator matches correctly."""
        for x in self.fixtures['should_not_match']:
            res = self.validator_instance().validate(x['token'])
            self.assertIsNone(res, msg=x)
