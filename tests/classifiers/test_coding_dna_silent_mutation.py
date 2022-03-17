"""Module for testing Coding DNA Silent Mutation Classifier."""
import unittest

from variation.classifiers import CodingDNASilentMutationClassifier
from .classifier_base import ClassifierBase


class TestCodingDNASilentMutationClassifier(ClassifierBase, unittest.TestCase):
    """A class to test the Coding DNA Silent Mutation Classifier."""

    def classifier_instance(self):
        """Return CodingDNASilentMutationClassifier instance."""
        return CodingDNASilentMutationClassifier()

    def fixture_name(self):
        """Return CodingDNASilentMutationClassifier fixture name."""
        return "coding_dna_silent_mutation"
