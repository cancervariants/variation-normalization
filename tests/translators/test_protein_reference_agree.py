"""Module for testing Reference Agree Translator."""
import unittest

from variation.classifiers import ProteinReferenceAgreeClassifier
from variation.translators import ProteinReferenceAgree
from variation.validators import ProteinReferenceAgree as PRA
from .translator_base import TranslatorBase


class TestProteinReferenceAgreeTranslator(TranslatorBase, unittest.TestCase):
    """A class to test the reference agree Translator."""

    def classifier_instance(self):
        """Return reference agree instance."""
        return ProteinReferenceAgreeClassifier()

    def validator_instance(self):
        """Return reference agree instance."""
        return PRA(*self.aa_params)

    def translator_instance(self):
        """Return reference agree instance."""
        return ProteinReferenceAgree()

    def fixture_name(self):
        """Return the fixture name for reference agree."""
        return "protein_reference_agree"
