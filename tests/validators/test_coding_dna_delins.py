"""Module for testing Coding DNA DelIns Validator."""
import unittest

from variation.validators import CodingDNADelIns
from variation.classifiers import CodingDNADelInsClassifier
from .validator_base import ValidatorBase


class TestCodingDNADelInsValidator(ValidatorBase, unittest.TestCase):
    """A class to test the Coding DNA DelIns Validator."""

    def validator_instance(self):
        """Return coding DNA delins instance."""
        return CodingDNADelIns(*self.params)

    def classifier_instance(self):
        """Return the coding DNA delins classifier instance."""
        return CodingDNADelInsClassifier()

    def fixture_name(self):
        """Return the fixture name for coding DNA delins."""
        return "coding_dna_delins"
