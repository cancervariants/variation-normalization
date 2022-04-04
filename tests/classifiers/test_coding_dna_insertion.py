"""Module for testing Coding DNA Insertion Classifier."""
import unittest

from variation.classifiers import CodingDNAInsertionClassifier
from .classifier_base import ClassifierBase


class TestCodingDNAInsertionClassifier(ClassifierBase, unittest.TestCase):
    """A class to test the Coding DNA Insertion Classifier."""

    def classifier_instance(self):
        """Return CodingDNAInsertionClassifier instance."""
        return CodingDNAInsertionClassifier()

    def fixture_name(self):
        """Return CodingDNAInsertionClassifier fixture name."""
        return "coding_dna_insertion"
