"""Module for testing DNA Coding Substitution Classifier."""
import unittest
from variant.classifiers import CodingDNASubstitutionClassifier
from .classifier_base import ClassifierBase


class TestCodingDNASubstitutionClassifier(ClassifierBase, unittest.TestCase):
    """A class to test the Coding DNA Substitution Classifier."""

    def classifier_instance(self):
        """Return CodingDNASubstitutionClassifier instance."""
        return CodingDNASubstitutionClassifier()

    def fixture_name(self):
        """Return DNACodingSubstitutionClassifier fixture name."""
        return 'coding_dna_substitution'
