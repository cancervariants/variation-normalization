"""Module for testing Genomic Duplication Translator."""
import unittest

from variation.classifiers import GenomicDuplicationClassifier
from variation.translators import GenomicDuplication
from variation.validators import GenomicDuplication as GD_V
from .translator_base import TranslatorBase


class TestGenomicDuplicationTranslator(TranslatorBase, unittest.TestCase):
    """A class to test the Genomic Duplication Translator."""

    def classifier_instance(self):
        """Return genomic duplication instance."""
        return GenomicDuplicationClassifier()

    def validator_instance(self):
        """Return genomic duplication instance."""
        return GD_V(*self.params)

    def translator_instance(self):
        """Return genomic duplication instance."""
        return GenomicDuplication()

    def fixture_name(self):
        """Return the fixture name for genomic duplication."""
        return "genomic_duplication"
