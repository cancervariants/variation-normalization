"""Module for testing Protein Insertion Classifier."""
import unittest

from variation.classifiers import ProteinInsertionClassifier
from .classifier_base import ClassifierBase


class TestProteinInsertionClassifier(ClassifierBase, unittest.TestCase):
    """A class to test the Protein Insertion Classifier."""

    def classifier_instance(self):
        """Return ProteinInsertionClassifier instance."""
        return ProteinInsertionClassifier()

    def fixture_name(self):
        """Return ProteinInsertionClassifier fixture name."""
        return "protein_insertion"
