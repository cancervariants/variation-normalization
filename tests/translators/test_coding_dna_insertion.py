"""Module for testing Coding DNA Insertion Translator."""
import unittest

from variation.classifiers import CodingDNAInsertionClassifier
from variation.translators import CdnaInsertion
from variation.validators import CdnaInsertion as CDNAD_V
from .translator_base import TranslatorBase


class TestCodingDNAInsertionTranslator(TranslatorBase, unittest.TestCase):
    """A class to test the Coding DNA Insertion Translator."""

    def classifier_instance(self):
        """Return coding DNA insertion instance."""
        return CodingDNAInsertionClassifier()

    def validator_instance(self):
        """Return coding DNA delins instance."""
        return CDNAD_V(*self.params)

    def translator_instance(self):
        """Return coding DNA insertion instance."""
        return CdnaInsertion()

    def fixture_name(self):
        """Return the fixture name for coding DNA insertion."""
        return "coding_dna_insertion"
