"""Module for testing Coding DNA DelIns Classifier."""
import unittest

from variation.classifiers import CodingDNADelInsClassifier
from .classifier_base import ClassifierBase


class TestCodingDNADelInsClassifier(ClassifierBase, unittest.TestCase):
    """A class to test the Coding DNA DelIns Classifier."""

    def classifier_instance(self):
        """Return CodingDNADelInsClassifier instance."""
        return CodingDNADelInsClassifier()

    def fixture_name(self):
        """Return CodingDNADelInsClassifier fixture name."""
        return "coding_dna_delins"
