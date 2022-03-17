"""Module for testing Genomic Silent Mutation Classifier."""
import unittest

from variation.classifiers import GenomicSilentMutationClassifier
from .classifier_base import ClassifierBase


class TestGenomicSilentMutationClassifier(ClassifierBase, unittest.TestCase):
    """A class to test the Genomic Silent Mutation Classifier."""

    def classifier_instance(self):
        """Return GenomicSilentMutationClassifier instance."""
        return GenomicSilentMutationClassifier()

    def fixture_name(self):
        """Return GenomicSilentMutationClassifier fixture name."""
        return "genomic_silent_mutation"
