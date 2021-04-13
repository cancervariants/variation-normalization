"""Module for testing DNA Coding Substitution Classifier."""
import unittest
from variant.classifiers import DNACodingSubstitutionClassifier
from .classifier_base import ClassifierBase


class TestDNACodingSubstitutionClassifier(ClassifierBase, unittest.TestCase):
    """A class to test the DNA Coding Substitution Classifier."""

    def classifier_instance(self):
        """Return DNACodingSubstitutionClassifier instance."""
        return DNACodingSubstitutionClassifier()

    def fixture_name(self):
        """Return DNACodingSubstitutionClassifier fixture name."""
        return 'dna_coding_substitution'
