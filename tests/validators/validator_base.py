"""A module for testing validator classes."""
import yaml
from tests import PROJECT_ROOT
from variation.tokenizers import Tokenize, GeneSymbol
from variation.tokenizers.caches import AminoAcidCache
from gene.query import QueryHandler as GeneQueryHandler


class ValidatorBase:
    """The validator base class."""

    def setUp(self):
        """Set up the test cases."""
        with open(f'{PROJECT_ROOT}/tests/fixtures/validators.yml') as stream:
            self.all_fixtures = yaml.safe_load(stream)
        self.fixtures = self.all_fixtures.get(
            self.fixture_name(),
            {'should_match': [], 'should_not_match': []}
        )
        self.tokenizer = Tokenize(AminoAcidCache(),
                                  GeneSymbol(GeneQueryHandler()))
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
            tokens = self.tokenizer.perform(x['query'], [])
            classification = self.classifier.match(tokens)
            validation_results = self.validator.validate(
                classification, normalize_endpoint=True
            )
            is_valid = False
            for vr in validation_results:
                if vr.is_valid:
                    is_valid = True
                    break
            self.assertTrue(is_valid, msg=x)
            self.assertIsNotNone(validation_results, msg=x)

    def test_not_matches(self):
        """Test that validator matches correctly."""
        for x in self.fixtures['should_not_match']:
            tokens = self.tokenizer.perform(x['query'], [])
            classification = self.classifier.match(tokens)
            validation_results = self.validator.validate(
                classification, normalize_endpoint=True
            )
            is_valid = False
            for vr in validation_results:
                if vr.is_valid:
                    is_valid = True
                    break
            self.assertFalse(is_valid, msg=x)
