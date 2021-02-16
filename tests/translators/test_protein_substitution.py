"""Module for testing Protein Substitution Translator."""
import unittest
from variant.classifiers import ProteinSubstitutionClassifier
from variant.translators import ProteinSubstitution
from variant.validators import ProteinSubstitution as PSUB_V
from .translator_base import TranslatorBase
from variant.data_sources import SeqRepoAccess, TranscriptMappings
from variant import PROJECT_ROOT


class TestProteinSubstitutionTranslator(TranslatorBase, unittest.TestCase):
    """A class to test the Protein Substitution Translator."""

    def classifier_instance(self):
        """Return protein substitution instance."""
        return ProteinSubstitutionClassifier()

    def validator_instance(self):
        """Return protein substitution instance."""
        return PSUB_V(SeqRepoAccess(
            f"{PROJECT_ROOT}/variant/data/seqrepo/latest"),
            TranscriptMappings(f"{PROJECT_ROOT}/variant/data"
                               f"/transcript_mapping.tsv"))

    def translator_instance(self):
        """Return protein substitution instance."""
        return ProteinSubstitution()

    def fixture_name(self):
        """Return the fixture name for protein substitution."""
        return 'protein_substitution'
