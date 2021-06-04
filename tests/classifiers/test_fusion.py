"""Module for testing Fusion classification."""
import unittest
from variation.classifiers import FusionClassifier
from .classifier_base import ClassifierBase


class TestGenePairClassifier(ClassifierBase, unittest.TestCase):
    """The gene pair classifier test class."""

    def classifier_instance(self):
        """Return fusion classifier instance."""
        return FusionClassifier()

    def fixture_name(self):
        """Return fusion fixture name."""
        return 'fusion'
