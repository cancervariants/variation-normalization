"""Module for testing Silent Mutation Translator."""
import unittest
from variant.classifiers import SilentMutationClassifier
from variant.translators import SilentMutation
from variant.validators import SilentMutation as SM_V
from .translator_base import TranslatorBase
from variant.tokenizers import GeneSymbol
from variant.tokenizers.caches import GeneSymbolCache, AminoAcidCache
from variant.data_sources import SeqRepoAccess, TranscriptMappings
from variant import SEQREPO_DATA_PATH, TRANSCRIPT_MAPPINGS_PATH


class TestSilentMutationTranslator(TranslatorBase, unittest.TestCase):
    """A class to test the silent mutation Translator."""

    def classifier_instance(self):
        """Return silent mutation instance."""
        return SilentMutationClassifier()

    def validator_instance(self):
        """Return silent mutation instance."""
        return SM_V(SeqRepoAccess(SEQREPO_DATA_PATH),
                    TranscriptMappings(TRANSCRIPT_MAPPINGS_PATH),
                    GeneSymbol(GeneSymbolCache()),
                    AminoAcidCache()
                    )

    def translator_instance(self):
        """Return silent mutation instance."""
        return SilentMutation()

    def fixture_name(self):
        """Return the fixture name for silent mutation."""
        return 'silent_mutation'
