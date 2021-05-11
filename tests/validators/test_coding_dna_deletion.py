"""Module for testing Coding DNA Deletion Validator."""
import unittest
from variant.validators import CodingDNADeletion
from variant.classifiers import CodingDNADeletionClassifier
from .validator_base import ValidatorBase
from variant.tokenizers import GeneSymbol
from variant.tokenizers.caches import GeneSymbolCache
from variant.data_sources import TranscriptMappings, SeqRepoAccess
from variant import SEQREPO_DATA_PATH, TRANSCRIPT_MAPPINGS_PATH


class TestCodingDNADeletionValidator(ValidatorBase, unittest.TestCase):
    """A class to test the CodingDNADeletion Validator."""

    def validator_instance(self):
        """Return coding dna deletion instance."""
        return CodingDNADeletion(SeqRepoAccess(SEQREPO_DATA_PATH),
                                 TranscriptMappings(TRANSCRIPT_MAPPINGS_PATH),  # noqa: E501
                                 GeneSymbol(GeneSymbolCache()))

    def classifier_instance(self):
        """Return the coding dna deletion classifier instance."""
        return CodingDNADeletionClassifier()

    def fixture_name(self):
        """Return the fixture name for coding dna deletion."""
        return 'coding_dna_deletion'
