"""Module for testing Coding DNA DelIns Translator."""
import unittest

from variation.classifiers import CodingDNADelInsClassifier
from variation.translators import CodingDNADelIns
from variation.validators import CodingDNADelIns as CDNADELINS_V
from .translator_base import TranslatorBase


class TestCodingDNADelInsTranslator(TranslatorBase, unittest.TestCase):
    """A class to test the Coding DNA DelIns Translator."""

    def classifier_instance(self):
        """Return coding DNA delins instance."""
        return CodingDNADelInsClassifier()

    def validator_instance(self):
        """Return coding DNA delins instance."""
        return CDNADELINS_V(*self.params)

    def translator_instance(self):
        """Return coding DNA delins instance."""
        return CodingDNADelIns()

    def fixture_name(self):
        """Return the fixture name for coding DNA delins."""
        return "coding_dna_delins"
