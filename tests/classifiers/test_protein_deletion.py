"""Module for testing Protein Deletion Classifier."""
import unittest

from variation.classifiers import ProteinDeletionClassifier
from .classifier_base import ClassifierBase


class TestProteinDeletionClassifier(ClassifierBase, unittest.TestCase):
    """A class to test the Protein Deletion Classifier."""

    def classifier_instance(self):
        """Return ProteinDeletionClassifier instance."""
        return ProteinDeletionClassifier()

    def fixture_name(self):
        """Return ProteinDeletionClassifier fixture name."""
        return "protein_deletion"
