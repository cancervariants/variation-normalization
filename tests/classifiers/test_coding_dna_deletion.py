"""Module for testing Coding DNA Deletion Classifier."""
import unittest

from variation.classifiers import CodingDNADeletionClassifier
from .classifier_base import ClassifierBase


class TestCodingDNADeletionClassifier(ClassifierBase, unittest.TestCase):
    """A class to test the Coding DNA Deletion Classifier."""

    def classifier_instance(self):
        """Return CodingDNADeletionClassifier instance."""
        return CodingDNADeletionClassifier()

    def fixture_name(self):
        """Return CodingDNADeletionClassifier fixture name."""
        return "coding_dna_deletion"
