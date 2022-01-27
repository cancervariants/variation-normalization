"""Module for testing amino acid deletion Translator."""
import unittest
from variation.classifiers import AminoAcidDeletionClassifier
from variation.translators import AminoAcidDeletion
from variation.validators import AminoAcidDeletion as AAD_V
from .translator_base import TranslatorBase


class TestAminoAcidDeletionTranslator(TranslatorBase, unittest.TestCase):
    """A class to test the amino acid deletion Translator."""

    def classifier_instance(self):
        """Return amino acid deletion instance."""
        return AminoAcidDeletionClassifier()

    def validator_instance(self):
        """Return amino acid deletion instance."""
        return AAD_V(*self.aa_params)

    def translator_instance(self):
        """Return amino acid deletion instance."""
        return AminoAcidDeletion()

    def fixture_name(self):
        """Return the fixture name for amino acid deletion."""
        return 'amino_acid_deletion'
