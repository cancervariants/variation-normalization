"""Module for testing Protein Insertion Validator."""
import unittest

from variation.validators import ProteinInsertion
from variation.classifiers import ProteinInsertionClassifier
from .validator_base import ValidatorBase


class TestProteinInsertionValidator(ValidatorBase, unittest.TestCase):
    """A class to test the Protein Insertion Validator."""

    def validator_instance(self):
        """Return protein insertion instance."""
        return ProteinInsertion(*self.aa_params)

    def classifier_instance(self):
        """Return the protein insertion classifier instance."""
        return ProteinInsertionClassifier()

    def fixture_name(self):
        """Return the fixture name for protein insertion."""
        return "protein_insertion"
