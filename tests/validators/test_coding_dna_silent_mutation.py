"""Module for testing Coding DNA Reference Agree Validator."""
import unittest

from variation.validators import CodingDNAReferenceAgree
from variation.classifiers import CodingDNAReferenceAgreeClassifier
from .validator_base import ValidatorBase


class TestCodingDNAReferenceAgreeValidator(ValidatorBase, unittest.TestCase):
    """A class to test the Coding DNA Reference Agree Validator."""

    def validator_instance(self):
        """Return coding DNA reference agree instance."""
        return CodingDNAReferenceAgree(*self.params)

    def classifier_instance(self):
        """Return the coding DNA reference agree classifier instance."""
        return CodingDNAReferenceAgreeClassifier()

    def fixture_name(self):
        """Return the fixture name for coding DNA reference agree."""
        return "coding_dna_reference_agree"
