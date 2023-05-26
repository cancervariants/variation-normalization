"""Module for testing Coding DNA Reference Agree Translator."""
import unittest

from variation.classifiers import CodingDNAReferenceAgreeClassifier
from variation.translators import CodingDNAReferenceAgree
from variation.validators import CodingDNAReferenceAgree as CDNASM_V
from .translator_base import TranslatorBase


class TestCodingDNAReferenceAgreeTranslator(TranslatorBase, unittest.TestCase):
    """A class to test the Coding DNA Reference Agree Translator."""

    def classifier_instance(self):
        """Return coding DNA reference agree instance."""
        return CodingDNAReferenceAgreeClassifier()

    def validator_instance(self):
        """Return coding DNA reference agree instance."""
        return CDNASM_V(*self.params)

    def translator_instance(self):
        """Return coding DNA reference agree instance."""
        return CodingDNAReferenceAgree()

    def fixture_name(self):
        """Return the fixture name for coding DNA reference agree."""
        return "coding_dna_reference_agree"
