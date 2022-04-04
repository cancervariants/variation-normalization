"""Module for testing Polypeptide Truncation Validator."""
import unittest

from variation.validators import PolypeptideTruncation
from variation.classifiers import PolypeptideTruncationClassifier
from .validator_base import ValidatorBase


class TestPolypeptideTruncationValidator(ValidatorBase, unittest.TestCase):
    """A class to test the Polypeptide Truncation Validator."""

    def validator_instance(self):
        """Return Polypeptide Truncation instance."""
        return PolypeptideTruncation(*self.aa_params)

    def classifier_instance(self):
        """Return the Polypeptide Truncation classifier instance."""
        return PolypeptideTruncationClassifier()

    def fixture_name(self):
        """Return the fixture name for Polypeptide Truncation."""
        return "polypeptide_truncation"
