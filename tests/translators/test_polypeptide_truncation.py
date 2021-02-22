"""Module for testing Protein Substitution Translator."""
import unittest
from variant.classifiers import PolypeptideTruncationClassifier
from variant.translators import PolypeptideTruncation
from variant.validators import PolypeptideTruncation as PT_V
from .translator_base import TranslatorBase
from variant.data_sources import SeqRepoAccess, TranscriptMappings
from variant import SEQREPO_DATA_PATH, TRANSCRIPT_MAPPINGS_PATH


class TestPolypeptideTruncationTranslator(TranslatorBase, unittest.TestCase):
    """A class to test the Protein Substitution Translator."""

    def classifier_instance(self):
        """Return protein substitution instance."""
        return PolypeptideTruncationClassifier()

    def validator_instance(self):
        """Return protein substitution instance."""
        return PT_V(SeqRepoAccess(SEQREPO_DATA_PATH),
                    TranscriptMappings(TRANSCRIPT_MAPPINGS_PATH))

    def translator_instance(self):
        """Return protein substitution instance."""
        return PolypeptideTruncation()

    def fixture_name(self):
        """Return the fixture name for protein substitution."""
        return 'amino_acid_substitution'
