"""Module for testing Coding DNA Substitution Validator."""
import unittest

from variation.validators import CodingDNASubstitution
from variation.classifiers import CodingDNASubstitutionClassifier
from .validator_base import ValidatorBase


class TestCodingDNASubstitutionValidator(ValidatorBase, unittest.TestCase):
    """A class to test the Coding DNA Substitution Validator."""

    def validator_instance(self):
        """Return coding DNA substitution instance."""
        return CodingDNASubstitution(*self.params)

    def classifier_instance(self):
        """Return the coding DNA substitution classifier instance."""
        return CodingDNASubstitutionClassifier()

    def fixture_name(self):
        """Return the fixture name for coding DNA substitution."""
        return "coding_dna_substitution"
