"""Module for testing Coding DNA Deletion Validator."""
import unittest

from variation.validators import CodingDNADeletion
from variation.classifiers import CodingDNADeletionClassifier
from .validator_base import ValidatorBase


class TestCodingDNADeletionValidator(ValidatorBase, unittest.TestCase):
    """A class to test the CodingDNADeletion Validator."""

    def validator_instance(self):
        """Return coding dna deletion instance."""
        return CodingDNADeletion(*self.params)

    def classifier_instance(self):
        """Return the coding dna deletion classifier instance."""
        return CodingDNADeletionClassifier()

    def fixture_name(self):
        """Return the fixture name for coding dna deletion."""
        return "coding_dna_deletion"
