"""A module for testing validator classes."""
import yaml
from varlexapp import PROJECT_ROOT
from varlexapp.tokenizers import Tokenize
from varlexapp.models import ClassificationType


class ValidatorBase(object):
    """The validator base class."""

    def setUp(self):
        """Set up the test cases."""
        with open(f'{PROJECT_ROOT}/tests/fixtures/validators.yml') as stream:
            self.all_fixtures = yaml.safe_load(stream)
        self.fixtures = self.all_fixtures.get(
                self.fixture_name(),
                {'tests': []}
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
        for x in self.fixtures['tests']:
            tokens = self.tokenizer.perform(x['query'])
            classification = self.classifier.match(tokens)
            validation_results = self.validator.validate(classification)
            for vr in validation_results:
                self.assertEqual(ClassificationType.PROTEIN_SUBSTITUTION,
                                 vr.classification.classification_type,
                                 msg=x)
                self.assertEqual(x['is_valid'], vr.is_valid, msg=x)
            self.assertIsNotNone(validation_results, msg=x)
