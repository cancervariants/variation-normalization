"""Module for testing protein deletion Translator."""
import unittest

from variation.classifiers import ProteinDeletionClassifier
from variation.translators import ProteinDeletion
from variation.validators import ProteinDeletion as AAD_V
from .translator_base import TranslatorBase


class TestProteinDeletionTranslator(TranslatorBase, unittest.TestCase):
    """A class to test the protein deletion Translator."""

    def classifier_instance(self):
        """Return protein deletion instance."""
        return ProteinDeletionClassifier()

    def validator_instance(self):
        """Return protein deletion instance."""
        return AAD_V(*self.aa_params)

    def translator_instance(self):
        """Return protein deletion instance."""
        return ProteinDeletion()

    def fixture_name(self):
        """Return the fixture name for protein deletion."""
        return "protein_deletion"
