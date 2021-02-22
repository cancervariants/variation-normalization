"""Module for testing Protein Substitution Translator."""
import unittest
from variant.classifiers import SilentMutationClassifier
from variant.translators import SilentMutation
from variant.validators import SilentMutation as SM_V
from .translator_base import TranslatorBase
from variant.data_sources import SeqRepoAccess, TranscriptMappings
from variant import SEQREPO_DATA_PATH, TRANSCRIPT_MAPPINGS_PATH


class TestAminoAcidSubstitutionTranslator(TranslatorBase, unittest.TestCase):
    """A class to test the Protein Substitution Translator."""

    def classifier_instance(self):
        """Return protein substitution instance."""
        return SilentMutationClassifier()

    def validator_instance(self):
        """Return protein substitution instance."""
        return SM_V(SeqRepoAccess(SEQREPO_DATA_PATH),
                    TranscriptMappings(TRANSCRIPT_MAPPINGS_PATH))

    def translator_instance(self):
        """Return protein substitution instance."""
        return SilentMutation()

    def fixture_name(self):
        """Return the fixture name for protein substitution."""
        return 'amino_acid_substitution'
