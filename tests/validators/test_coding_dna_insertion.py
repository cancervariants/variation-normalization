"""Module for testing Coding DNA Insertion Validator."""
import unittest
from variation.validators import CodingDNAInsertion
from variation.classifiers import CodingDNAInsertionClassifier
from .validator_base import ValidatorBase
from variation.tokenizers import GeneSymbol
from variation.tokenizers.caches import GeneSymbolCache
from variation.data_sources import TranscriptMappings, SeqRepoAccess
from variation import SEQREPO_DATA_PATH, TRANSCRIPT_MAPPINGS_PATH


class TestCodingDNAInsertionValidator(ValidatorBase, unittest.TestCase):
    """A class to test the CodingDNAInsertion Validator."""

    def validator_instance(self):
        """Return coding dna insertion instance."""
        return CodingDNAInsertion(SeqRepoAccess(SEQREPO_DATA_PATH),
                                  TranscriptMappings(TRANSCRIPT_MAPPINGS_PATH),
                                  GeneSymbol(GeneSymbolCache()))

    def classifier_instance(self):
        """Return the coding dna insertion classifier instance."""
        return CodingDNAInsertionClassifier()

    def fixture_name(self):
        """Return the fixture name for coding dna insertion."""
        return 'coding_dna_insertion'
