"""Module for testing protein insertion Translator."""
import unittest

from variation.classifiers import ProteinInsertionClassifier
from variation.translators import ProteinInsertion
from variation.validators import ProteinInsertion as AAI_V
from .translator_base import TranslatorBase


class TestProteinInsertionTranslator(TranslatorBase, unittest.TestCase):
    """A class to test the protein insertion Translator."""

    def classifier_instance(self):
        """Return protein insertion instance."""
        return ProteinInsertionClassifier()

    def validator_instance(self):
        """Return protein insertion instance."""
        return AAI_V(*self.aa_params)

    def translator_instance(self):
        """Return protein insertion instance."""
        return ProteinInsertion()

    def fixture_name(self):
        """Return the fixture name for protein insertion."""
        return "protein_insertion"
