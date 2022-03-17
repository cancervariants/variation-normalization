"""Module for testing Silent Mutation Classifier."""
import unittest

from variation.classifiers import SilentMutationClassifier
from .classifier_base import ClassifierBase


class TestSilentMutationClassifier(ClassifierBase, unittest.TestCase):
    """A class to test the Silent Mutation Classifier."""

    def classifier_instance(self):
        """Return SilentMutationClassifier instance."""
        return SilentMutationClassifier()

    def fixture_name(self):
        """Return SilentMutationClassifier fixture name."""
        return "silent_mutation"
