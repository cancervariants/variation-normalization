"""Module for testing Coding DNA DelIns Translator."""
import unittest
from variant.classifiers import CodingDNADelInsClassifier
from variant.translators import CodingDNADelIns
from variant.validators import CodingDNADelIns as CDNADELINS_V
from .translator_base import TranslatorBase
from variant.tokenizers import GeneSymbol
from variant.tokenizers.caches import GeneSymbolCache
from variant.data_sources import SeqRepoAccess, TranscriptMappings
from variant import SEQREPO_DATA_PATH, TRANSCRIPT_MAPPINGS_PATH


class TestCodingDNADelInsTranslator(TranslatorBase, unittest.TestCase):
    """A class to test the Coding DNA DelIns Translator."""

    def classifier_instance(self):
        """Return coding DNA delins instance."""
        return CodingDNADelInsClassifier()

    def validator_instance(self):
        """Return coding DNA delins instance."""
        return CDNADELINS_V(SeqRepoAccess(SEQREPO_DATA_PATH),
                            TranscriptMappings(TRANSCRIPT_MAPPINGS_PATH),
                            GeneSymbol(GeneSymbolCache())
                            )

    def translator_instance(self):
        """Return coding DNA delins instance."""
        return CodingDNADelIns()

    def fixture_name(self):
        """Return the fixture name for coding DNA delins."""
        return 'coding_dna_delins'
