"""A module for testing translator classes."""
import yaml
from varlexapp import PROJECT_ROOT


class TranslatorBase(object):
    """The translator base class."""

    def setUp(self):
        """Set up the test cases."""
        with open(f'{PROJECT_ROOT}/tests/fixtures/translators.yml') as stream:
            self.all_fixtures = yaml.safe_load(stream)
        self.fixtures = self.all_fixtures.get(
                self.fixture_name(),
                {'should_match': [], 'should_not_match': []}
        )

    def translator_instance(self):
        """Check that the translator_instance method is implemented."""
        raise NotImplementedError()

    def fixture_name(self):
        """Check that the fixture_name method is implemented."""
        raise NotImplementedError()

    def test_matches(self):
        """Test that translator matches correctly."""
        for x in self.fixtures['should_match']:
            res = self.translator_instance().perform(x['query'])  # noqa: F841

    def test_not_matches(self):
        """Test that translator matches correctly."""
        for x in self.fixtures['should_not_match']:
            res = self.translator_instance().perform(x['query'])  # noqa: F841
