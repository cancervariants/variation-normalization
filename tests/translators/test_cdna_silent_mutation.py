"""Module for testing Cdna Reference Agree Translator."""
import unittest

from variation.classifiers import CdnaReferenceAgreeClassifier
from variation.translators import CdnaReferenceAgree
from variation.validators import CdnaReferenceAgree as CDNASM_V
from .translator_base import TranslatorBase


class TestCdnaReferenceAgreeTranslator(TranslatorBase, unittest.TestCase):
    """A class to test the Cdna Reference Agree Translator."""

    def classifier_instance(self):
        """Return cDNA reference agree instance."""
        return CdnaReferenceAgreeClassifier()

    def validator_instance(self):
        """Return cDNA reference agree instance."""
        return CDNASM_V(*self.params)

    def translator_instance(self):
        """Return cDNA reference agree instance."""
        return CdnaReferenceAgree()

    def fixture_name(self):
        """Return the fixture name for cDNA reference agree."""
        return "cdna_reference_agree"
