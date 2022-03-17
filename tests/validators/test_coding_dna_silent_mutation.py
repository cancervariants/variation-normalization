"""Module for testing Coding DNA Silent Mutation Validator."""
import unittest

from variation.validators import CodingDNASilentMutation
from variation.classifiers import CodingDNASilentMutationClassifier
from .validator_base import ValidatorBase


class TestCodingDNASilentMutationValidator(ValidatorBase, unittest.TestCase):
    """A class to test the Coding DNA Silent Mutation Validator."""

    def validator_instance(self):
        """Return coding DNA silent mutation instance."""
        return CodingDNASilentMutation(*self.params)

    def classifier_instance(self):
        """Return the coding DNA silent mutation classifier instance."""
        return CodingDNASilentMutationClassifier()

    def fixture_name(self):
        """Return the fixture name for coding DNA silent mutation."""
        return "coding_dna_silent_mutation"
