"""Module for testing Silent Mutation Translator."""
import unittest
from variation.classifiers import SilentMutationClassifier
from variation.translators import SilentMutation
from variation.validators import SilentMutation as SM_V
from .translator_base import TranslatorBase
from variation.tokenizers import GeneSymbol
from variation.tokenizers.caches import GeneSymbolCache, AminoAcidCache
from variation.data_sources import TranscriptMappings, SeqRepoAccess, \
    MANETranscriptMappings, UTA
from variation.mane_transcript import MANETranscript


class TestSilentMutationTranslator(TranslatorBase, unittest.TestCase):
    """A class to test the silent mutation Translator."""

    def classifier_instance(self):
        """Return silent mutation instance."""
        return SilentMutationClassifier()

    def validator_instance(self):
        """Return silent mutation instance."""
        seqrepo_access = SeqRepoAccess()
        transcript_mappings = TranscriptMappings()
        uta = UTA()
        return SM_V(
            seqrepo_access, transcript_mappings, GeneSymbol(GeneSymbolCache()),
            MANETranscript(seqrepo_access, transcript_mappings,
                           MANETranscriptMappings(), uta),
            uta,
            AminoAcidCache())

    def translator_instance(self):
        """Return silent mutation instance."""
        return SilentMutation()

    def fixture_name(self):
        """Return the fixture name for silent mutation."""
        return 'silent_mutation'
