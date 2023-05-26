"""Module for testing Coding DNA Reference Agree Classifier."""
import unittest

from variation.classifiers import CodingDNAReferenceAgreeClassifier
from .classifier_base import ClassifierBase


class TestCodingDNAReferenceAgreeClassifier(ClassifierBase, unittest.TestCase):
    """A class to test the Coding DNA Reference Agree Classifier."""

    def classifier_instance(self):
        """Return CodingDNAReferenceAgreeClassifier instance."""
        return CodingDNAReferenceAgreeClassifier()

    def fixture_name(self):
        """Return CodingDNAReferenceAgreeClassifier fixture name."""
        return "coding_dna_reference_agree"
