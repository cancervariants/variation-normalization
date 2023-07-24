"""Module for testing Cdna Deletion Validator."""
import unittest

from variation.validators import CdnaDeletion
from variation.classifiers import CdnaDeletionClassifier
from .validator_base import ValidatorBase


class TestCdnaDeletionValidator(ValidatorBase, unittest.TestCase):
    """A class to test the CdnaDeletion Validator."""

    def validator_instance(self):
        """Return cdna deletion instance."""
        return CdnaDeletion(*self.params)

    def classifier_instance(self):
        """Return the cdna deletion classifier instance."""
        return CdnaDeletionClassifier()

    def fixture_name(self):
        """Return the fixture name for cdna deletion."""
        return "cdna_deletion"
