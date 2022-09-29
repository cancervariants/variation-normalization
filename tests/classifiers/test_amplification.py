"""Module for testing Amplification classifier"""
import unittest

from variation.classifiers import AmplificationClassifier
from .classifier_base import ClassifierBase


class TestAmplificationClassifier(ClassifierBase, unittest.TestCase):
    """A class to test the Amplification Classifier"""

    def classifier_instance(self):
        """Return AmplificationClassifier instance"""
        return AmplificationClassifier()

    def fixture_name(self):
        """Return AmplificationClassifier fixture name"""
        return "amplification"
