"""Module for testing Coding DNA Substitution Translator."""
import unittest

from variation.classifiers import CodingDNASubstitutionClassifier
from variation.translators import CodingDNASubstitution
from variation.validators import CodingDNASubstitution as CDNASUB_V
from .translator_base import TranslatorBase


class TestCodingDNASubstitutionTranslator(TranslatorBase, unittest.TestCase):
    """A class to test the Coding DNA Substitution Translator."""

    def classifier_instance(self):
        """Return coding DNA substitution instance."""
        return CodingDNASubstitutionClassifier()

    def validator_instance(self):
        """Return coding DNA substitution instance."""
        return CDNASUB_V(*self.params)

    def translator_instance(self):
        """Return coding DNA substitution instance."""
        return CodingDNASubstitution()

    def fixture_name(self):
        """Return the fixture name for coding DNA substitution."""
        return "coding_dna_substitution"
