"""Module for testing Genomic Reference Agree Translator."""
import unittest

from variation.classifiers import GenomicReferenceAgreeClassifier
from variation.translators import GenomicReferenceAgree
from variation.validators import GenomicReferenceAgree as GENOMICSM_V
from .translator_base import TranslatorBase


class TestGenomicReferenceAgreeTranslator(TranslatorBase, unittest.TestCase):
    """A class to test the Genomic Reference Agree Translator."""

    def classifier_instance(self):
        """Return genomic reference agree instance."""
        return GenomicReferenceAgreeClassifier()

    def validator_instance(self):
        """Return genomic reference agree instance."""
        return GENOMICSM_V(*self.params)

    def translator_instance(self):
        """Return genomic reference agree instance."""
        return GenomicReferenceAgree()

    def fixture_name(self):
        """Return the fixture name for genomic reference agree."""
        return "genomic_reference_agree"
