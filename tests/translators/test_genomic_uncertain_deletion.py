"""Module for testing Genomic Uncertain Deletion Translator."""
import unittest

from variation.classifiers import GenomicUncertainDeletionClassifier
from variation.translators import GenomicUncertainDeletion
from variation.validators import GenomicUncertainDeletion as GUD_V
from .translator_base import TranslatorBase


class TestGenomicUncertainDeletionTranslator(TranslatorBase,
                                             unittest.TestCase):
    """A class to test the Genomic Uncertain Deletion Translator."""

    def classifier_instance(self):
        """Return genomic uncertain deletion instance."""
        return GenomicUncertainDeletionClassifier()

    def validator_instance(self):
        """Return genomic uncertain deletion instance."""
        return GUD_V(*self.params)

    def translator_instance(self):
        """Return genomic uncertain deletion instance."""
        return GenomicUncertainDeletion()

    def fixture_name(self):
        """Return the fixture name for genomic uncertain deletion."""
        return "genomic_uncertain_deletion"
