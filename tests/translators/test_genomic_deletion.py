"""Module for testing Genomic Deletion Translator."""
import unittest

from variation.classifiers import GenomicDeletionClassifier
from variation.translators import GenomicDeletion
from variation.validators import GenomicDeletion as GD_V
from .translator_base import TranslatorBase


class TestGenomicDeletionTranslator(TranslatorBase, unittest.TestCase):
    """A class to test the Genomic Deletion Translator."""

    def classifier_instance(self):
        """Return genomic deletion instance."""
        return GenomicDeletionClassifier()

    def validator_instance(self):
        """Return coding DNA deletion instance."""
        return GD_V(*self.params)

    def translator_instance(self):
        """Return genomic deletion instance."""
        return GenomicDeletion()

    def fixture_name(self):
        """Return the fixture name for genomic deletion."""
        return "genomic_deletion"
