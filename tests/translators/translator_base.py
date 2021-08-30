"""A module for testing translator classes."""
import yaml
from tests import PROJECT_ROOT
from variation.tokenizers import Tokenize, GeneSymbol
from variation.tokenizers.caches import AminoAcidCache
from gene.query import QueryHandler as GeneQueryHandler


class TranslatorBase:
    """The translator base class."""

    def setUp(self):
        """Set up the test cases."""
        with open(f'{PROJECT_ROOT}/tests/fixtures/translators.yml') as stream:
            self.all_fixtures = yaml.safe_load(stream)
        self.fixtures = self.all_fixtures.get(
            self.fixture_name(),
            {'tests': []}
        )
        amino_acid_cache = AminoAcidCache()
        gene_symbol = GeneSymbol(GeneQueryHandler())
        self.tokenizer = Tokenize(amino_acid_cache, gene_symbol)
        self.classifier = self.classifier_instance()
        self.validator = self.validator_instance()
        self.translator = self.translator_instance()

    def classifier_instance(self):
        """Check that the classifier_instance method is implemented."""
        raise NotImplementedError()

    def validator_instance(self):
        """Check that the validator_instance method is implemented."""
        raise NotImplementedError()

    def translator_instance(self):
        """Check that the translator_instance method is implemented."""
        raise NotImplementedError()

    def fixture_name(self):
        """Check that the fixture_name method is implemented."""
        raise NotImplementedError()

    def test_translator(self):
        """Test that translator matches correctly."""
        for x in self.fixtures['tests']:
            tokens = self.tokenizer.perform(x['query'], [])
            classification = self.classifier.match(tokens)
            validation_results = self.validator.validate(
                classification, normalize_endpoint=False
            )
            num_valid = 0
            found = list()
            for vr in validation_results:
                if vr.is_valid:
                    variation = (self.translator.translate(vr)).__dict__
                    if variation['type'] == 'Allele':
                        variation['location'] = variation['location'].dict()
                        del variation['location']['id']
                    elif variation['type'] == 'CopyNumber':
                        variation['subject'] = variation['subject'].dict()
                        del variation['subject']['location']['id']

                    if variation not in found:
                        found.append(variation)
                        num_valid += 1
                        self.assertIn(variation, x['variations'],
                                      msg=x['query'])
            self.assertEqual(len(x['variations']), num_valid, msg=x['query'])
