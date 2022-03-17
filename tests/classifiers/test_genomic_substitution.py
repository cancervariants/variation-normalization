"""Module for testing Genomic Substitution Classifier."""
import unittest

from variation.classifiers import GenomicSubstitutionClassifier
from .classifier_base import ClassifierBase


class TestGenomicSubstitutionClassifier(ClassifierBase, unittest.TestCase):
    """A class to test the Coding DNA Substitution Classifier."""

    def classifier_instance(self):
        """Return GenomicSubstitutionClassifier instance."""
        return GenomicSubstitutionClassifier()

    def fixture_name(self):
        """Return GenomicSubstitutionClassifier fixture name."""
        return "genomic_substitution"
