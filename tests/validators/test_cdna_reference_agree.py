"""Module for testing Cdna Reference Agree Validator."""
import unittest

from variation.validators import CdnaReferenceAgree
from variation.classifiers import CdnaReferenceAgreeClassifier
from .validator_base import ValidatorBase


class TestCdnaReferenceAgreeValidator(ValidatorBase, unittest.TestCase):
    """A class to test the Cdna Reference Agree Validator."""

    def validator_instance(self):
        """Return cDNA reference agree instance."""
        return CdnaReferenceAgree(*self.params)

    def classifier_instance(self):
        """Return the cDNA reference agree classifier instance."""
        return CdnaReferenceAgreeClassifier()

    def fixture_name(self):
        """Return the fixture name for cDNA reference agree."""
        return "cdna_reference_agree"
