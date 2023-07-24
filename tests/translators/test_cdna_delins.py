"""Module for testing Cdna DelIns Translator."""
import unittest

from variation.classifiers import CdnaDelInsClassifier
from variation.translators import CdnaDelIns
from variation.validators import CdnaDelIns as CDNADELINS_V
from .translator_base import TranslatorBase


class TestCdnaDelInsTranslator(TranslatorBase, unittest.TestCase):
    """A class to test the Cdna DelIns Translator."""

    def classifier_instance(self):
        """Return cDNA delins instance."""
        return CdnaDelInsClassifier()

    def validator_instance(self):
        """Return cDNA delins instance."""
        return CDNADELINS_V(*self.params)

    def translator_instance(self):
        """Return cDNA delins instance."""
        return CdnaDelIns()

    def fixture_name(self):
        """Return the fixture name for cDNA delins."""
        return "cdna_delins"
