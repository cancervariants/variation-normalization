"""Module for testing Genomic Insertion Classifier."""
import unittest

from variation.classifiers import GenomicInsertionClassifier
from .classifier_base import ClassifierBase


class TestGenomicInsertionClassifier(ClassifierBase, unittest.TestCase):
    """A class to test the Genomic Insertion Classifier."""

    def classifier_instance(self):
        """Return GenomicInsertionClassifier instance."""
        return GenomicInsertionClassifier()

    def fixture_name(self):
        """Return GenomicInsertionClassifier fixture name."""
        return "genomic_insertion"
