"""Module for testing Genomic Deletion Range Classifier."""
import unittest

from variation.classifiers import GenomicDeletionRangeClassifier
from .classifier_base import ClassifierBase


class TestGenomicDeletionRangeClassifier(ClassifierBase,
                                         unittest.TestCase):
    """A class to test the Genomic Deletion Range Classifier."""

    def classifier_instance(self):
        """Return Genomic Deletion Range Classifier instance."""
        return GenomicDeletionRangeClassifier()

    def fixture_name(self):
        """Return Genomic Deletion Range fixture name."""
        return "genomic_deletion_range"
