"""Module for testing Cdna Substitution Validator."""
import unittest

from variation.validators import CdnaSubstitution
from variation.classifiers import CdnaSubstitutionClassifier
from .validator_base import ValidatorBase


class TestCdnaSubstitutionValidator(ValidatorBase, unittest.TestCase):
    """A class to test the Cdna Substitution Validator."""

    def validator_instance(self):
        """Return cDNA substitution instance."""
        return CdnaSubstitution(*self.params)

    def classifier_instance(self):
        """Return the cDNA substitution classifier instance."""
        return CdnaSubstitutionClassifier()

    def fixture_name(self):
        """Return the fixture name for cDNA substitution."""
        return "cdna_substitution"
