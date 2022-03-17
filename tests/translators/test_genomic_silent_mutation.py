"""Module for testing Genomic Silent Mutation Translator."""
import unittest

from variation.classifiers import GenomicSilentMutationClassifier
from variation.translators import GenomicSilentMutation
from variation.validators import GenomicSilentMutation as GENOMICSM_V
from .translator_base import TranslatorBase


class TestGenomicSilentMutationTranslator(TranslatorBase, unittest.TestCase):
    """A class to test the Genomic Silent Mutation Translator."""

    def classifier_instance(self):
        """Return genomic silent mutation instance."""
        return GenomicSilentMutationClassifier()

    def validator_instance(self):
        """Return genomic silent mutation instance."""
        return GENOMICSM_V(*self.params)

    def translator_instance(self):
        """Return genomic silent mutation instance."""
        return GenomicSilentMutation()

    def fixture_name(self):
        """Return the fixture name for genomic silent mutation."""
        return "genomic_silent_mutation"
