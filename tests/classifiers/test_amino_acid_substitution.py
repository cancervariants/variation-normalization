"""Module for testing Amino Acid Substitution Classifier."""
import unittest
from variation.classifiers import AminoAcidSubstitutionClassifier
from .classifier_base import ClassifierBase


class TestAminoAcidSubstitutionClassifier(ClassifierBase, unittest.TestCase):
    """A class to test the Amino Acid Substitution Classifier."""

    def classifier_instance(self):
        """Return AminoAcidSubstitutionClassifier instance."""
        return AminoAcidSubstitutionClassifier()

    def fixture_name(self):
        """Return AminoAcidSubstitutionClassifier fixture name."""
        return 'amino_acid_substitution'
