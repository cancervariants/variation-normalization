"""Module for testing protein delins Translator."""
import unittest

from variation.classifiers import ProteinDelInsClassifier
from variation.translators import ProteinDelIns
from variation.validators import ProteinDelIns as AAD_V
from .translator_base import TranslatorBase


class TestProteinDelInsTranslator(TranslatorBase, unittest.TestCase):
    """A class to test the protein delins Translator."""

    def classifier_instance(self):
        """Return protein delins instance."""
        return ProteinDelInsClassifier()

    def validator_instance(self):
        """Return protein delins instance."""
        return AAD_V(*self.aa_params)

    def translator_instance(self):
        """Return protein delins instance."""
        return ProteinDelIns()

    def fixture_name(self):
        """Return the fixture name for protein delins."""
        return "protein_delins"
