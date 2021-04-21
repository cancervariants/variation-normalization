"""Module for testing DNA Coding DelIns Classifier."""
import unittest
from variant.classifiers import CodingDNADelInsClassifier
from .classifier_base import ClassifierBase


class TestCodingDNADelInsClassifier(ClassifierBase, unittest.TestCase):
    """A class to test the Coding DNA DelIns Classifier."""

    def classifier_instance(self):
        """Return CodingDNADelInsClassifier instance."""
        return CodingDNADelInsClassifier()

    def fixture_name(self):
        """Return DNACodingDelInsClassifier fixture name."""
        return 'coding_dna_delins'
