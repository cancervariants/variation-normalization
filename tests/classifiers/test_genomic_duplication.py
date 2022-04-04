"""Module for testing Genomic Duplication Classifier."""
import unittest

from variation.classifiers import GenomicDuplicationClassifier
from .classifier_base import ClassifierBase


class TestGenomicDuplicationClassifier(ClassifierBase, unittest.TestCase):
    """A class to test the Genomic Duplication Classifier."""

    def classifier_instance(self):
        """Return GenomicDuplicationClassifier instance."""
        return GenomicDuplicationClassifier()

    def fixture_name(self):
        """Return GenomicDuplicationClassifier fixture name."""
        return "genomic_duplication"
