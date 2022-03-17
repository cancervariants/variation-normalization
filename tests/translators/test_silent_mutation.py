"""Module for testing Silent Mutation Translator."""
import unittest

from variation.classifiers import SilentMutationClassifier
from variation.translators import SilentMutation
from variation.validators import SilentMutation as SM_V
from .translator_base import TranslatorBase


class TestSilentMutationTranslator(TranslatorBase, unittest.TestCase):
    """A class to test the silent mutation Translator."""

    def classifier_instance(self):
        """Return silent mutation instance."""
        return SilentMutationClassifier()

    def validator_instance(self):
        """Return silent mutation instance."""
        return SM_V(*self.aa_params)

    def translator_instance(self):
        """Return silent mutation instance."""
        return SilentMutation()

    def fixture_name(self):
        """Return the fixture name for silent mutation."""
        return "silent_mutation"
