"""Module for testing Genomic Uncertain Deletion Classifier."""
import unittest

from variation.classifiers import GenomicUncertainDeletionClassifier
from .classifier_base import ClassifierBase


class TestGenomicUncertainDeletionClassifier(ClassifierBase,
                                             unittest.TestCase):
    """A class to test the Genomic Uncertain Deletion Classifier."""

    def classifier_instance(self):
        """Return Genomic Uncertain Deletion Classifier instance."""
        return GenomicUncertainDeletionClassifier()

    def fixture_name(self):
        """Return Genomic Uncertain Deletion fixture name."""
        return "genomic_uncertain_deletion"
