"""Module for testing amino acid insertion Translator."""
import unittest
from variation.classifiers import AminoAcidInsertionClassifier
from variation.translators import AminoAcidInsertion
from variation.validators import AminoAcidInsertion as AAI_V
from .translator_base import TranslatorBase


class TestAminoAcidInsertionTranslator(TranslatorBase, unittest.TestCase):
    """A class to test the amino acid insertion Translator."""

    def classifier_instance(self):
        """Return amino acid insertion instance."""
        return AminoAcidInsertionClassifier()

    def validator_instance(self):
        """Return amino acid insertion instance."""
        return AAI_V(*self.aa_params)

    def translator_instance(self):
        """Return amino acid insertion instance."""
        return AminoAcidInsertion()

    def fixture_name(self):
        """Return the fixture name for amino acid insertion."""
        return 'amino_acid_insertion'
