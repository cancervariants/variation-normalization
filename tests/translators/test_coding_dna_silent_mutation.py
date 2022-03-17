"""Module for testing Coding DNA Silent Mutation Translator."""
import unittest

from variation.classifiers import CodingDNASilentMutationClassifier
from variation.translators import CodingDNASilentMutation
from variation.validators import CodingDNASilentMutation as CDNASM_V
from .translator_base import TranslatorBase


class TestCodingDNASilentMutationTranslator(TranslatorBase, unittest.TestCase):
    """A class to test the Coding DNA Silent Mutation Translator."""

    def classifier_instance(self):
        """Return coding DNA silent mutation instance."""
        return CodingDNASilentMutationClassifier()

    def validator_instance(self):
        """Return coding DNA silent mutation instance."""
        return CDNASM_V(*self.params)

    def translator_instance(self):
        """Return coding DNA silent mutation instance."""
        return CodingDNASilentMutation()

    def fixture_name(self):
        """Return the fixture name for coding DNA silent mutation."""
        return "coding_dna_silent_mutation"
