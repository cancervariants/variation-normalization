"""Module for testing Genomic Insertion Translator."""
import unittest

from variation.classifiers import GenomicInsertionClassifier
from variation.translators import GenomicInsertion
from variation.validators import GenomicInsertion as GD_V
from .translator_base import TranslatorBase


class TestGenomicInsertionTranslator(TranslatorBase, unittest.TestCase):
    """A class to test the Genomic Insertion Translator."""

    def classifier_instance(self):
        """Return genomic insertion instance."""
        return GenomicInsertionClassifier()

    def validator_instance(self):
        """Return coding DNA insertion instance."""
        return GD_V(*self.params)

    def translator_instance(self):
        """Return genomic insertion instance."""
        return GenomicInsertion()

    def fixture_name(self):
        """Return the fixture name for genomic insertion."""
        return "genomic_insertion"
