"""Module for testing Coding DNA Deletion Translator."""
import unittest
from variant.classifiers import CodingDNADeletionClassifier
from variant.translators import CodingDNADeletion
from variant.validators import CodingDNADeletion as CDNAD_V
from .translator_base import TranslatorBase
from variant.tokenizers import GeneSymbol
from variant.tokenizers.caches import GeneSymbolCache
from variant.data_sources import SeqRepoAccess, TranscriptMappings
from variant import SEQREPO_DATA_PATH, TRANSCRIPT_MAPPINGS_PATH


class TestCodingDNADeletionTranslator(TranslatorBase, unittest.TestCase):
    """A class to test the Coding DNA Deletion Translator."""

    def classifier_instance(self):
        """Return coding DNA deletion instance."""
        return CodingDNADeletionClassifier()

    def validator_instance(self):
        """Return coding DNA delins instance."""
        return CDNAD_V(SeqRepoAccess(SEQREPO_DATA_PATH),
                       TranscriptMappings(TRANSCRIPT_MAPPINGS_PATH),
                       GeneSymbol(GeneSymbolCache())
                       )

    def translator_instance(self):
        """Return coding DNA deletion instance."""
        return CodingDNADeletion()

    def fixture_name(self):
        """Return the fixture name for coding DNA deletion."""
        return 'coding_dna_deletion'
