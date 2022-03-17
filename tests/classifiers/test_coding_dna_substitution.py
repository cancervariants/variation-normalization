"""Module for testing Coding DNA Substitution Classifier."""
import unittest

from variation.classifiers import CodingDNASubstitutionClassifier
from .classifier_base import ClassifierBase


class TestCodingDNASubstitutionClassifier(ClassifierBase, unittest.TestCase):
    """A class to test the Coding DNA Substitution Classifier."""

    def classifier_instance(self):
        """Return CodingDNASubstitutionClassifier instance."""
        return CodingDNASubstitutionClassifier()

    def fixture_name(self):
        """Return CodingDNASubstitutionClassifier fixture name."""
        return "coding_dna_substitution"
