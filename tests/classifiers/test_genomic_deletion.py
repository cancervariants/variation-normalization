"""Module for testing Genomic Deletion Classifier."""
import unittest

from variation.classifiers import GenomicDeletionClassifier
from .classifier_base import ClassifierBase


class TestGenomicDeletionClassifier(ClassifierBase, unittest.TestCase):
    """A class to test the Genomic Deletion Classifier."""

    def classifier_instance(self):
        """Return GenomicDeletionClassifier instance."""
        return GenomicDeletionClassifier()

    def fixture_name(self):
        """Return GenomicDeletionClassifier fixture name."""
        return "genomic_deletion"
