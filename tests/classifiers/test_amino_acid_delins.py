"""Module for testing Amino Acid DelIns Classifier."""
import unittest
from variation.classifiers import AminoAcidDelInsClassifier
from .classifier_base import ClassifierBase


class TestAminoAcidDelInsClassifier(ClassifierBase, unittest.TestCase):
    """A class to test the Amino Acid DelIns Classifier."""

    def classifier_instance(self):
        """Return AminoAcidDelInsClassifier instance."""
        return AminoAcidDelInsClassifier()

    def fixture_name(self):
        """Return AminoAcidDelInsClassifier fixture name."""
        return 'amino_acid_delins'
