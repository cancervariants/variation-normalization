"""Module for testing Coding DNA Insertion Translator."""
import unittest
from variant.classifiers import CodingDNAInsertionClassifier
from variant.translators import CodingDNAInsertion
from variant.validators import CodingDNAInsertion as CDNAD_V
from .translator_base import TranslatorBase
from variant.tokenizers import GeneSymbol
from variant.tokenizers.caches import GeneSymbolCache
from variant.data_sources import SeqRepoAccess, TranscriptMappings
from variant import SEQREPO_DATA_PATH, TRANSCRIPT_MAPPINGS_PATH


class TestCodingDNAInsertionTranslator(TranslatorBase, unittest.TestCase):
    """A class to test the Coding DNA Insertion Translator."""

    def classifier_instance(self):
        """Return coding DNA insertion instance."""
        return CodingDNAInsertionClassifier()

    def validator_instance(self):
        """Return coding DNA delins instance."""
        return CDNAD_V(SeqRepoAccess(SEQREPO_DATA_PATH),
                       TranscriptMappings(TRANSCRIPT_MAPPINGS_PATH),
                       GeneSymbol(GeneSymbolCache())
                       )

    def translator_instance(self):
        """Return coding DNA insertion instance."""
        return CodingDNAInsertion()

    def fixture_name(self):
        """Return the fixture name for coding DNA insertion."""
        return 'coding_dna_insertion'
