"""Module for testing amino acid Substitution Translator."""
import unittest
from variation.classifiers import AminoAcidSubstitutionClassifier
from variation.translators import AminoAcidSubstitution
from variation.validators import AminoAcidSubstitution as AASUB_V
from .translator_base import TranslatorBase


class TestAminoAcidSubstitutionTranslator(TranslatorBase, unittest.TestCase):
    """A class to test the amino acid Substitution Translator."""

    def classifier_instance(self):
        """Return amino acid substitution instance."""
        return AminoAcidSubstitutionClassifier()

    def validator_instance(self):
        """Return amino acid substitution instance."""
        return AASUB_V(*self.aa_params)

    def translator_instance(self):
        """Return amino acid substitution instance."""
        return AminoAcidSubstitution()

    def fixture_name(self):
        """Return the fixture name for amino acid substitution."""
        return 'amino_acid_substitution'
