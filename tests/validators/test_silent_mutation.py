"""Module for testing Reference Agree Validator."""
import unittest

from variation.validators import ProteinReferenceAgree
from variation.classifiers import ProteinReferenceAgreeClassifier
from .validator_base import ValidatorBase


class TestReferenceAgreeValidator(ValidatorBase, unittest.TestCase):
    """A class to test the Reference Agree Validator."""

    def validator_instance(self):
        """Return Reference Agree instance."""
        return ProteinReferenceAgree(*self.aa_params)

    def classifier_instance(self):
        """Return the Reference Agree classifier instance."""
        return ProteinReferenceAgreeClassifier()

    def fixture_name(self):
        """Return the fixture name for Reference Agree."""
        return "protein_reference_agree"
