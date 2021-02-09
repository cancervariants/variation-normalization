"""Module for testing Protein Substitution Translator."""
import unittest
from varlexapp.translators import ProteinSubstitution
from .translator_base import TranslatorBase


class TestProteinSubstitutionTranslator(TranslatorBase, unittest.TestCase):
    """A class to test the Protein Substitution Translator."""

    def translator_instance(self):
        """Return protein substitution instance."""
        return ProteinSubstitution()

    def fixture_name(self):
        """Return the fixture name for protein substitution."""
        return 'protein_substitution'
