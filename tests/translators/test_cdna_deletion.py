"""Module for testing Cdna Deletion Translator."""
import unittest

from variation.classifiers import CdnaDeletionClassifier
from variation.translators import CdnaDeletion
from variation.validators import CdnaDeletion as CDNAD_V
from .translator_base import TranslatorBase


class TestCdnaDeletionTranslator(TranslatorBase, unittest.TestCase):
    """A class to test the Cdna Deletion Translator."""

    def classifier_instance(self):
        """Return cDNA deletion instance."""
        return CdnaDeletionClassifier()

    def validator_instance(self):
        """Return cDNA delins instance."""
        return CDNAD_V(*self.params)

    def translator_instance(self):
        """Return cDNA deletion instance."""
        return CdnaDeletion()

    def fixture_name(self):
        """Return the fixture name for cDNA deletion."""
        return "cdna_deletion"
