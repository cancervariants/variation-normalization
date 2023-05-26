"""Module for testing Reference Agree Classifier."""
import unittest

from variation.classifiers import ProteinReferenceAgreeClassifier
from .classifier_base import ClassifierBase


class TestReferenceAgreeClassifier(ClassifierBase, unittest.TestCase):
    """A class to test the Reference Agree Classifier."""

    def classifier_instance(self):
        """Return ProteinReferenceAgreeClassifier instance."""
        return ProteinReferenceAgreeClassifier()

    def fixture_name(self):
        """Return ProteinReferenceAgreeClassifier fixture name."""
        return "reference_agree"
