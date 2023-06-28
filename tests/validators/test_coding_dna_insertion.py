"""Module for testing Coding DNA Insertion Validator."""
import unittest

from variation.validators import CdnaInsertion
from variation.classifiers import CodingDNAInsertionClassifier
from .validator_base import ValidatorBase


class TestCodingDNAInsertionValidator(ValidatorBase, unittest.TestCase):
    """A class to test the CdnaInsertion Validator."""

    def validator_instance(self):
        """Return coding dna insertion instance."""
        return CdnaInsertion(*self.params)

    def classifier_instance(self):
        """Return the coding dna insertion classifier instance."""
        return CodingDNAInsertionClassifier()

    def fixture_name(self):
        """Return the fixture name for coding dna insertion."""
        return "coding_dna_insertion"
