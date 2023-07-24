"""Module for testing Cdna Insertion Translator."""
import unittest

from variation.classifiers import CdnaInsertionClassifier
from variation.translators import CdnaInsertion
from variation.validators import CdnaInsertion as CDNAD_V
from .translator_base import TranslatorBase


class TestCdnaInsertionTranslator(TranslatorBase, unittest.TestCase):
    """A class to test the Cdna Insertion Translator."""

    def classifier_instance(self):
        """Return cDNA insertion instance."""
        return CdnaInsertionClassifier()

    def validator_instance(self):
        """Return cDNA delins instance."""
        return CDNAD_V(*self.params)

    def translator_instance(self):
        """Return cDNA insertion instance."""
        return CdnaInsertion()

    def fixture_name(self):
        """Return the fixture name for cDNA insertion."""
        return "cdna_insertion"
