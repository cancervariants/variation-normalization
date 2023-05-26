"""Module for testing Genomic Reference Agree Classifier."""
import unittest

from variation.classifiers import GenomicReferenceAgreeClassifier
from .classifier_base import ClassifierBase


class TestGenomicReferenceAgreeClassifier(ClassifierBase, unittest.TestCase):
    """A class to test the Genomic Reference Agree Classifier."""

    def classifier_instance(self):
        """Return GenomicReferenceAgreeClassifier instance."""
        return GenomicReferenceAgreeClassifier()

    def fixture_name(self):
        """Return GenomicReferenceAgreeClassifier fixture name."""
        return "genomic_reference_agree"
