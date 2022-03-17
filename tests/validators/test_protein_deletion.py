"""Module for testing Protein Deletion Validator."""
import unittest

from variation.validators import ProteinDeletion
from variation.classifiers import ProteinDeletionClassifier
from .validator_base import ValidatorBase


class TestProteinDeletionValidator(ValidatorBase, unittest.TestCase):
    """A class to test the Protein Deletion Validator."""

    def validator_instance(self):
        """Return protein deletion instance."""
        return ProteinDeletion(*self.aa_params)

    def classifier_instance(self):
        """Return the protein deletion classifier instance."""
        return ProteinDeletionClassifier()

    def fixture_name(self):
        """Return the fixture name for protein deletion."""
        return "protein_deletion"
