"""Module for testing Genomic DelIns Classifier."""
import unittest

from variation.classifiers import GenomicDelInsClassifier
from .classifier_base import ClassifierBase


class TestGenomicDelInsClassifier(ClassifierBase, unittest.TestCase):
    """A class to test the Genomic DelIns Classifier."""

    def classifier_instance(self):
        """Return GenomicDelInsClassifier instance."""
        return GenomicDelInsClassifier()

    def fixture_name(self):
        """Return GenomicDelInsClassifier fixture name."""
        return "genomic_delins"
