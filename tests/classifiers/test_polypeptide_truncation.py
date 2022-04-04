"""Module for testing Polypeptide Truncation Classifier."""
import unittest

from variation.classifiers import PolypeptideTruncationClassifier
from .classifier_base import ClassifierBase


class TestPolypeptideTruncationClassifier(ClassifierBase, unittest.TestCase):
    """A class to test the Polypeptide Truncation Classifier."""

    def classifier_instance(self):
        """Return PolypeptideTruncationClassifier instance."""
        return PolypeptideTruncationClassifier()

    def fixture_name(self):
        """Return PolypeptideTruncationClassifier fixture name."""
        return "polypeptide_truncation"
