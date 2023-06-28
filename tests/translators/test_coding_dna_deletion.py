"""Module for testing Coding DNA Deletion Translator."""
import unittest

from variation.classifiers import CdnaDeletionClassifier
from variation.translators import CdnaDeletion
from variation.validators import CdnaDeletion as CDNAD_V
from .translator_base import TranslatorBase


class TestCodingDNADeletionTranslator(TranslatorBase, unittest.TestCase):
    """A class to test the Coding DNA Deletion Translator."""

    def classifier_instance(self):
        """Return coding DNA deletion instance."""
        return CdnaDeletionClassifier()

    def validator_instance(self):
        """Return coding DNA delins instance."""
        return CDNAD_V(*self.params)

    def translator_instance(self):
        """Return coding DNA deletion instance."""
        return CdnaDeletion()

    def fixture_name(self):
        """Return the fixture name for coding DNA deletion."""
        return "coding_dna_deletion"
