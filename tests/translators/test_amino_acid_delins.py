"""Module for testing amino acid delins Translator."""
import unittest
from variation.classifiers import AminoAcidDelInsClassifier
from variation.translators import AminoAcidDelIns
from variation.validators import AminoAcidDelIns as AAD_V
from .translator_base import TranslatorBase


class TestAminoAcidDelInsTranslator(TranslatorBase, unittest.TestCase):
    """A class to test the amino acid delins Translator."""

    def classifier_instance(self):
        """Return amino acid delins instance."""
        return AminoAcidDelInsClassifier()

    def validator_instance(self):
        """Return amino acid delins instance."""
        return AAD_V(*self.aa_params)

    def translator_instance(self):
        """Return amino acid delins instance."""
        return AminoAcidDelIns()

    def fixture_name(self):
        """Return the fixture name for amino acid delins."""
        return 'amino_acid_delins'
