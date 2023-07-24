"""Module for testing Cdna DelIns Validator."""
import unittest

from variation.validators import CdnaDelIns
from variation.classifiers import CdnaDelInsClassifier
from .validator_base import ValidatorBase


class TestCdnaDelInsValidator(ValidatorBase, unittest.TestCase):
    """A class to test the Cdna DelIns Validator."""

    def validator_instance(self):
        """Return cDNA delins instance."""
        return CdnaDelIns(*self.params)

    def classifier_instance(self):
        """Return the cDNA delins classifier instance."""
        return CdnaDelInsClassifier()

    def fixture_name(self):
        """Return the fixture name for cDNA delins."""
        return "cdna_delins"
