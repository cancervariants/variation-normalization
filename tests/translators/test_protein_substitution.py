"""Module for testing Protein Substitution Translator."""
import unittest
from varlexapp.classifiers import ProteinSubstitutionClassifier
from varlexapp.translators import ProteinSubstitution
from varlexapp.validators import ProteinSubstitution as PSUB_V
from .translator_base import TranslatorBase
from varlexapp.data_sources import SeqRepoAccess, TranscriptMappings
from varlexapp import PROJECT_ROOT


class TestProteinSubstitutionTranslator(TranslatorBase, unittest.TestCase):
    """A class to test the Protein Substitution Translator."""

    def classifier_instance(self):
        """Return protein substitution instance."""
        return ProteinSubstitutionClassifier()

    def validator_instance(self):
        """Return protein substitution instance."""
        return PSUB_V(SeqRepoAccess(
            f"{PROJECT_ROOT}/varlexapp/data/seqrepo/latest"),
            TranscriptMappings(f"{PROJECT_ROOT}/varlexapp/data"
                               f"/transcript_mapping.tsv"))

    def translator_instance(self):
        """Return protein substitution instance."""
        return ProteinSubstitution()

    def fixture_name(self):
        """Return the fixture name for protein substitution."""
        return 'protein_substitution'
