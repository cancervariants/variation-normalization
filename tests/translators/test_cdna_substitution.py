"""Module for testing Cdna Substitution Translator."""
import unittest

from variation.classifiers import CdnaSubstitutionClassifier
from variation.translators import CdnaSubstitution
from variation.validators import CdnaSubstitution as CDNASUB_V
from .translator_base import TranslatorBase


class TestCdnaSubstitutionTranslator(TranslatorBase, unittest.TestCase):
    """A class to test the Cdna Substitution Translator."""

    def classifier_instance(self):
        """Return cDNA substitution instance."""
        return CdnaSubstitutionClassifier()

    def validator_instance(self):
        """Return cDNA substitution instance."""
        return CDNASUB_V(*self.params)

    def translator_instance(self):
        """Return cDNA substitution instance."""
        return CdnaSubstitution()

    def fixture_name(self):
        """Return the fixture name for cDNA substitution."""
        return "cdna_substitution"
