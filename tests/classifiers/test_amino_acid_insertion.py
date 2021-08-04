"""Module for testing Amino Acid Insertion Classifier."""
import unittest
from variation.classifiers import AminoAcidInsertionClassifier
from .classifier_base import ClassifierBase


class TestAminoAcidInsertionClassifier(ClassifierBase, unittest.TestCase):
    """A class to test the Amino Acid Insertion Classifier."""

    def classifier_instance(self):
        """Return AminoAcidInsertionClassifier instance."""
        return AminoAcidInsertionClassifier()

    def fixture_name(self):
        """Return AminoAcidInsertionClassifier fixture name."""
        return 'amino_acid_insertion'
