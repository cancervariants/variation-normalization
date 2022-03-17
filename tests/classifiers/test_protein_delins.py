"""Module for testing Protein DelIns Classifier."""
import unittest

from variation.classifiers import ProteinDelInsClassifier
from .classifier_base import ClassifierBase


class TestProteinDelInsClassifier(ClassifierBase, unittest.TestCase):
    """A class to test the Protein DelIns Classifier."""

    def classifier_instance(self):
        """Return ProteinDelInsClassifier instance."""
        return ProteinDelInsClassifier()

    def fixture_name(self):
        """Return ProteinDelInsClassifier fixture name."""
        return "protein_delins"
