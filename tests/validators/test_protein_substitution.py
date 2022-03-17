"""Module for testing Protein Substitution Validator."""
import unittest

from variation.validators import ProteinSubstitution
from variation.classifiers import ProteinSubstitutionClassifier
from .validator_base import ValidatorBase


class TestProteinSubstitutionValidator(ValidatorBase, unittest.TestCase):
    """A class to test the Protein Substitution Validator."""

    def validator_instance(self):
        """Return protein substitution instance."""
        return ProteinSubstitution(*self.aa_params)

    def classifier_instance(self):
        """Return the protein substitution classifier instance."""
        return ProteinSubstitutionClassifier()

    def fixture_name(self) -> str:
        """Return the fixture name for protein substitution."""
        return "protein_substitution"
