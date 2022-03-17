"""Module for testing polypeptide truncation Translator."""
import unittest

from variation.classifiers import PolypeptideTruncationClassifier
from variation.translators import PolypeptideTruncation
from variation.validators import PolypeptideTruncation as PT_V
from .translator_base import TranslatorBase


class TestPolypeptideTruncationTranslator(TranslatorBase, unittest.TestCase):
    """A class to test the polypeptide truncation Translator."""

    def classifier_instance(self):
        """Return polypeptide truncation instance."""
        return PolypeptideTruncationClassifier()

    def validator_instance(self):
        """Return polypeptide truncation instance."""
        return PT_V(*self.aa_params)

    def translator_instance(self):
        """Return polypeptide truncation instance."""
        return PolypeptideTruncation()

    def fixture_name(self):
        """Return the fixture name for polypeptide truncation."""
        return "polypeptide_truncation"
