"""Module for testing Cdna Insertion Validator."""
import unittest

from variation.validators import CdnaInsertion
from variation.classifiers import CdnaInsertionClassifier
from .validator_base import ValidatorBase


class TestCdnaInsertionValidator(ValidatorBase, unittest.TestCase):
    """A class to test the CdnaInsertion Validator."""

    def validator_instance(self):
        """Return cdna insertion instance."""
        return CdnaInsertion(*self.params)

    def classifier_instance(self):
        """Return the cdna insertion classifier instance."""
        return CdnaInsertionClassifier()

    def fixture_name(self):
        """Return the fixture name for cdna insertion."""
        return "cdna_insertion"
