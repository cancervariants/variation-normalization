"""Module for testing Genomic DelIns Translator."""
import unittest

from variation.classifiers import GenomicDelInsClassifier
from variation.translators import GenomicDelIns
from variation.validators import GenomicDelIns as GENOMICDELINS_V
from .translator_base import TranslatorBase


class TestGenomicDelInsTranslator(TranslatorBase, unittest.TestCase):
    """A class to test the Genomic DelIns Translator."""

    def classifier_instance(self):
        """Return Genomic delins instance."""
        return GenomicDelInsClassifier()

    def validator_instance(self):
        """Return genomic delins instance."""
        return GENOMICDELINS_V(*self.params)

    def translator_instance(self):
        """Return genomic delins instance."""
        return GenomicDelIns()

    def fixture_name(self):
        """Return the fixture name for genomic delins."""
        return "genomic_delins"
