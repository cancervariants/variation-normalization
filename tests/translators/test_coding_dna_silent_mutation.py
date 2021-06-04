"""Module for testing Coding DNA Silent Mutation Translator."""
import unittest
from variation.classifiers import CodingDNASilentMutationClassifier
from variation.translators import CodingDNASilentMutation
from variation.validators import CodingDNASilentMutation as CDNASM_V
from .translator_base import TranslatorBase
from variation.tokenizers import GeneSymbol
from variation.tokenizers.caches import GeneSymbolCache
from variation.data_sources import SeqRepoAccess, TranscriptMappings
from variation import SEQREPO_DATA_PATH, TRANSCRIPT_MAPPINGS_PATH


class TestCodingDNASilentMutationTranslator(TranslatorBase, unittest.TestCase):
    """A class to test the Coding DNA Silent Mutation Translator."""

    def classifier_instance(self):
        """Return coding DNA silent mutation instance."""
        return CodingDNASilentMutationClassifier()

    def validator_instance(self):
        """Return coding DNA silent mutation instance."""
        return CDNASM_V(SeqRepoAccess(SEQREPO_DATA_PATH),
                        TranscriptMappings(TRANSCRIPT_MAPPINGS_PATH),
                        GeneSymbol(GeneSymbolCache())
                        )

    def translator_instance(self):
        """Return coding DNA silent mutation instance."""
        return CodingDNASilentMutation()

    def fixture_name(self):
        """Return the fixture name for coding DNA silent mutation."""
        return 'coding_dna_silent_mutation'
