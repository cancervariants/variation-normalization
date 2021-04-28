"""Module for testing Genomic Silent Mutation Translator."""
import unittest
from variant.classifiers import GenomicSilentMutationClassifier
from variant.translators import GenomicSilentMutation
from variant.validators import GenomicSilentMutation as GENOMICSM_V
from .translator_base import TranslatorBase
from variant.tokenizers import GeneSymbol
from variant.tokenizers.caches import GeneSymbolCache
from variant.data_sources import SeqRepoAccess, TranscriptMappings
from variant import SEQREPO_DATA_PATH, TRANSCRIPT_MAPPINGS_PATH


class TestGenomicSilentMutationTranslator(TranslatorBase, unittest.TestCase):
    """A class to test the Genomic Silent Mutation Translator."""

    def classifier_instance(self):
        """Return genomic silent mutation instance."""
        return GenomicSilentMutationClassifier()

    def validator_instance(self):
        """Return genomic silent mutation instance."""
        return GENOMICSM_V(SeqRepoAccess(SEQREPO_DATA_PATH),
                           TranscriptMappings(TRANSCRIPT_MAPPINGS_PATH),
                           GeneSymbol(GeneSymbolCache())
                           )

    def translator_instance(self):
        """Return genomic silent mutation instance."""
        return GenomicSilentMutation()

    def fixture_name(self):
        """Return the fixture name for genomic silent mutation."""
        return 'genomic_silent_mutation'
