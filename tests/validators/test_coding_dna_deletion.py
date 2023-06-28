"""Module for testing Coding DNA Deletion Validator."""
import unittest

from variation.validators import CdnaDeletion
from variation.classifiers import CdnaDeletionClassifier
from .validator_base import ValidatorBase


class TestCodingDNADeletionValidator(ValidatorBase, unittest.TestCase):
    """A class to test the CdnaDeletion Validator."""

    def validator_instance(self):
        """Return coding dna deletion instance."""
        return CdnaDeletion(*self.params)

    def classifier_instance(self):
        """Return the coding dna deletion classifier instance."""
        return CdnaDeletionClassifier()

    def fixture_name(self):
        """Return the fixture name for coding dna deletion."""
        return "coding_dna_deletion"
