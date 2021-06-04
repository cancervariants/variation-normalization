"""Module for testing Amino Acid Deletion Classifier."""
import unittest
from variation.classifiers import AminoAcidDeletionClassifier
from .classifier_base import ClassifierBase


class TestAminoAcidDeletionClassifier(ClassifierBase, unittest.TestCase):
    """A class to test the Amino Acid Deletion Classifier."""

    def classifier_instance(self):
        """Return AminoAcidDeletionClassifier instance."""
        return AminoAcidDeletionClassifier()

    def fixture_name(self):
        """Return AminoAcidDeletionClassifier fixture name."""
        return 'amino_acid_deletion'
