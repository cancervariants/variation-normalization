"""Module for testing Silent Mutation Validator."""
import unittest

from variation.validators import SilentMutation
from variation.classifiers import SilentMutationClassifier
from .validator_base import ValidatorBase


class TestSilentMutationValidator(ValidatorBase, unittest.TestCase):
    """A class to test the Silent Mutation Validator."""

    def validator_instance(self):
        """Return Silent Mutation instance."""
        return SilentMutation(*self.aa_params)

    def classifier_instance(self):
        """Return the Silent Mutation classifier instance."""
        return SilentMutationClassifier()

    def fixture_name(self):
        """Return the fixture name for Silent Mutation."""
        return "silent_mutation"
