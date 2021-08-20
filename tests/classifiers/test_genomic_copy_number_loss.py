"""Module for testing Genomic Copy Number Loss Classifier."""
import unittest
from variation.classifiers import GenomicCopyNumberLossClassifier
from .classifier_base import ClassifierBase


class TestGenomicCopyNumberLossClassifier(ClassifierBase, unittest.TestCase):
    """A class to test the Genomic Copy Number Loss Classifier."""

    def classifier_instance(self):
        """Return Genomic Copy Number Loss Classifier instance."""
        return GenomicCopyNumberLossClassifier()

    def fixture_name(self):
        """Return Genomic Copy Number Loss fixture name."""
        return 'genomic_copy_number_loss'
