"""Module for testing Protein DelIns Validator."""
import unittest

from variation.validators import ProteinDelIns
from variation.classifiers import ProteinDelInsClassifier
from .validator_base import ValidatorBase


class TestProteinDelInsValidator(ValidatorBase, unittest.TestCase):
    """A class to test the Protein DelIns Validator."""

    def validator_instance(self):
        """Return protein delins instance."""
        return ProteinDelIns(*self.aa_params)

    def classifier_instance(self):
        """Return the protein delins classifier instance."""
        return ProteinDelInsClassifier()

    def fixture_name(self):
        """Return the fixture name for protein delins."""
        return "protein_delins"
