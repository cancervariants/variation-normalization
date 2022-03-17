"""Module for testing protein Substitution Translator."""
import unittest

from variation.classifiers import ProteinSubstitutionClassifier
from variation.translators import ProteinSubstitution
from variation.validators import ProteinSubstitution as AASUB_V
from .translator_base import TranslatorBase


class TestProteinSubstitutionTranslator(TranslatorBase, unittest.TestCase):
    """A class to test the protein Substitution Translator."""

    def classifier_instance(self):
        """Return protein substitution instance."""
        return ProteinSubstitutionClassifier()

    def validator_instance(self):
        """Return protein substitution instance."""
        return AASUB_V(*self.aa_params)

    def translator_instance(self):
        """Return protein substitution instance."""
        return ProteinSubstitution()

    def fixture_name(self):
        """Return the fixture name for protein substitution."""
        return "protein_substitution"
