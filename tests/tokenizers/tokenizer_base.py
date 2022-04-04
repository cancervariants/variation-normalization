"""A module for testing tokenizer classes."""
import yaml

from tests import PROJECT_ROOT


class TokenizerBase:
    """The tokenizer base class."""

    def setUp(self):
        """Set up the test test cases."""
        with open(f"{PROJECT_ROOT}/tests/fixtures/tokenizers.yml") as stream:
            self.all_fixtures = yaml.safe_load(stream)
        self.fixtures = self.all_fixtures.get(
            self.fixture_name(),
            {"should_match": [], "should_not_match": []}
        )

    def tokenizer_instance(self):
        """Check that tokenizer_instance method is implemented."""
        raise NotImplementedError()

    def token_type(self):
        """Check that token_type method is implemented."""
        raise NotImplementedError()

    def fixture_name(self):
        """Check that fixture_name method is implemented."""
        raise NotImplementedError()

    def test_matches(self):
        """Test that tokenizer matches correctly."""
        for x in self.fixtures["should_match"]:
            res = self.tokenizer_instance().match(x["token"])
            self.assertIsNotNone(res, msg=x)
            self.assertEqual(self.token_type(), res.token_type)

    def test_not_matches(self):
        """Test that tokenizer matches correctly."""
        for x in self.fixtures["should_not_match"]:
            res = self.tokenizer_instance().match(x["token"])
            try:
                self.assertIsNone(res, msg=x)
            except AssertionError:
                assert self.token_type() != res.token_type
